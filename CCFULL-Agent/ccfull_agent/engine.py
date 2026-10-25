"""Pinned builds, review receipts and isolated execution. No network in this module."""
import copy
import datetime
import hashlib
import json
import math
import re
import shutil
import subprocess
import uuid
from pathlib import Path
from .model import ValidationError, require, validate, render, digest, read_json, energy_grid

ROOT = Path(__file__).resolve().parents[1]
VENDOR = ROOT/"vendor"/"ccfull"
SOURCE_SHA = "548144eac5b8d8fb80b474dc0c686f745a4a54385df80702fdba5c2cdd925823"
INPUT_SHA = "cd6da9eae4372eb08f6cfcde4d775597a994d3b92bb5062af40259ff67de857a"
OUTPUT_SHA = "cc4038926fe349ccb6e60fba2572f731bbd2e3ca95abee321b99fe4f21ce367d"
FLAGS = ["-std=legacy", "-O2"]
NORMAL = "The program terminated normally"
FLOAT = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[EeDd][+-]?\d+)?"

def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+"\n", encoding="utf-8")

def verify_vendor():
    for name, expected in (("ccfull.f", SOURCE_SHA), ("ccfull.inp", INPUT_SHA), ("OUTPUT", OUTPUT_SHA)):
        require(sha(VENDOR/name) == expected, f"Pinned upstream {name} checksum mismatch")
    return SOURCE_SHA

def build():
    verify_vendor()
    compiler = shutil.which("gfortran")
    require(compiler is not None, "gfortran missing: install a Fortran compiler first")
    version = subprocess.run([compiler, "--version"], capture_output=True, text=True, check=True).stdout.splitlines()[0]
    identity = {"source_sha256": SOURCE_SHA, "compiler": compiler, "compiler_version": version, "flags": FLAGS}
    directory = ROOT/"build"/digest(identity)[:16]
    binary = directory/"ccfull"
    if (directory/"build.json").exists():
        record = read_json(directory/"build.json")
        require(record["identity"] == identity and sha(binary) == record["binary_sha256"], "Build cache changed; inspect it before rebuilding")
        return binary, record
    directory.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run([compiler, *FLAGS, "-o", str(binary), str(VENDOR/"ccfull.f")],
                          capture_output=True, text=True, timeout=120)
    (directory/"compiler.stdout.log").write_text(proc.stdout)
    (directory/"compiler.stderr.log").write_text(proc.stderr)
    require(proc.returncode == 0, f"Compilation failed; see {directory}")
    record = {"identity": identity, "binary_sha256": sha(binary), "built_at": now()}
    save(directory/"build.json", record)
    return binary, record

def make_plan(configs, purpose):
    require(isinstance(configs, list) and 1 <= len(configs) <= 30, "Plan must contain 1..30 configurations")
    require(isinstance(purpose, str) and bool(purpose.strip()), "Plan purpose required")
    ids = []
    for c in configs:
        validate(c)
        ids.append(c["id"])
    require(len(ids) == len(set(ids)), "Duplicate configuration ids")
    return {"schema_version": 1, "source_sha256": SOURCE_SHA, "purpose": purpose, "configs": copy.deepcopy(configs)}

def validate_plan(plan):
    require(set(plan) == {"schema_version","source_sha256","purpose","configs"}, "Unexpected plan fields")
    require(plan["schema_version"] == 1 and plan["source_sha256"] == SOURCE_SHA, "Plan/source version mismatch")
    make_plan(plan["configs"], plan["purpose"])
    return digest(plan)

def prepare(plan, directory):
    token = validate_plan(plan)
    d = Path(directory)
    d.mkdir(parents=True, exist_ok=False)
    save(d/"plan.json", plan)
    review = ["# CCFULL calculation review", "", f"Plan digest: `{token}`", "", plan["purpose"], "",
              "**Status: proposed; no research calculation has run.**", "",
              "Confirm the exact plan before creating an approval receipt. Any edit invalidates that receipt.", ""]
    for i,c in enumerate(plan["configs"]):
        info = validate(c)
        deck, answers, prompts = render(c)
        run = d/f"{i:02d}"
        run.mkdir()
        save(run/"config.json", c)
        (run/"ccfull.inp").write_text(deck)
        (run/"stdin.txt").write_text(answers)
        save(run/"prompts.json", prompts)
        review += [f"## {c['id']}", "", f"Potential: {c['potential']['id']}; coupling: {c['couplings']['id']}",
                   f"Channels: {info['channels']}; points: {info['points']}", "",
                   "```json", json.dumps(c,ensure_ascii=False,indent=2), "```", ""]
    (d/"review.md").write_text("\n".join(review), encoding="utf-8")
    return {"directory": str(d), "plan_sha256": token, "status": "proposed"}

def approve(plan, statement, actor):
    token = validate_plan(plan)
    require(isinstance(statement,str) and len(statement.strip()) >= 4, "Record the user's actual confirmation text")
    require(isinstance(actor,str) and bool(actor.strip()), "Confirmation actor required")
    return {"schema_version":1, "plan_sha256":token, "actor":actor, "statement":statement,
            "confirmed_at":now(), "scope":"exact-plan",
            "note":"Workflow receipt, not identity authentication. Agent must not fabricate human approval."}

def check_approval(plan, receipt):
    token = validate_plan(plan)
    require(receipt.get("schema_version") == 1 and receipt.get("scope") == "exact-plan", "Approval receipt required")
    require(receipt.get("plan_sha256") == token, "Missing or stale approval: discuss the changed plan first")
    require(bool(receipt.get("actor")) and bool(receipt.get("statement")), "Incomplete approval receipt")

def parse_output(path):
    lines = Path(path).read_text(errors="replace").splitlines()
    barrier = None
    rows = []
    in_table = False
    for line in lines:
        match = re.search(r"Rb=\s*("+FLOAT+r").*Vb=\s*("+FLOAT+r").*Curv=\s*("+FLOAT+r")",line)
        if match:
            barrier = dict(zip(("Rb_fm","Vb_MeV","curvature_MeV"), [float(s.replace("D","E").replace("d","e")) for s in match.groups()]))
        if "Ecm (MeV)" in line and "sigma" in line:
            in_table = True
            continue
        if in_table and line.strip() and re.match(r"^[\s+-.0-9]",line) and any(ch.isdigit() for ch in line):
            parts = line.split()
            if len(parts)==3:
                try:
                    values = [float(x.replace("D","E").replace("d","e")) for x in parts]
                except ValueError as exc:
                    raise ValidationError("Invalid OUTPUT numeric row: "+line) from exc
                require(all(math.isfinite(x) for x in values), "Non-finite OUTPUT")
                rows.append(dict(zip(("E_cm_MeV","sigma_mb","mean_l_hbar"),values)))
    require(barrier is not None and rows, "Missing barrier or excitation table")
    require(all(math.isfinite(v) and v>0 for v in barrier.values()), "Nonphysical barrier")
    return {"barrier":barrier,"table":rows}

def parse_columns(path):
    rows = []
    for line in Path(path).read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split()
        require(len(parts)==2, f"{path}: expected two columns")
        row = [float(s.replace("D","E").replace("d","e")) for s in parts]
        require(all(math.isfinite(x) for x in row), f"{path}: non-finite output")
        rows.append(row)
    return rows

def collect_run(directory, config, *, probability_upper_tolerance=1e-6):
    # Normal execution remains strict. Only explicit retrospective C0 review uses the bounded override.
    require(1e-6<=probability_upper_tolerance<=2e-4, "Unsupported probability review tolerance")
    d = Path(directory)
    stdout = (d/"stdout.log").read_text(errors="replace")
    require(NORMAL in stdout, "Missing normal termination marker (Fortran STOP may exit zero)")
    result = parse_output(d/"OUTPUT")
    grid = energy_grid(config)
    cross, spin = parse_columns(d/"cross.dat"), parse_columns(d/"spin.dat")
    require(len(result["table"]) == len(grid) == len(cross) == len(spin), "Incomplete/extra energy rows")
    for expected, row, cr, sp in zip(grid,result["table"],cross,spin):
        for actual in (row["E_cm_MeV"],cr[0],sp[0]):
            require(math.isclose(expected,actual,rel_tol=0,abs_tol=5.1e-6), "Energy grid mismatch")
        require(cr[1]>=0 and 0<=sp[1]<=200, "Invalid cross section/mean angular momentum")
        require(math.isclose(row["sigma_mb"],cr[1],rel_tol=6e-6,abs_tol=5.1e-6), "OUTPUT/cross.dat mismatch")
        require(math.isclose(row["mean_l_hbar"],sp[1],rel_tol=0,abs_tol=5.1e-6), "OUTPUT/spin.dat mismatch")
        row.update(E_cm_MeV=expected,sigma_mb=cr[1],mean_l_hbar=sp[1])
    # This pinned source writes one line per partial wave, despite the filename.
    penetration = parse_columns(d/"s-wave.dat")
    require(bool(penetration), "Missing penetration output")
    require(all(-1e-8<=p<=1+probability_upper_tolerance for _,p in penetration), "Penetrability outside [0,1]: investigate numerical accuracy")
    result["warnings"] = sorted(set(line.strip() for line in stdout.splitlines()
                                   if any(s in line.lower() for s in ("warning","shallow potential","oscillatio"))))
    result["warnings"] += validate(config)["warnings"]
    if any(row["mean_l_hbar"]>150 for row in result["table"]):
        result["warnings"].append("Large mean l: fixed l<=200 cutoff requires independent checking")
    result["notes"] = ["spin.dat contains mean l, not a spin distribution.",
                       "s-wave.dat contains partial-wave penetrabilities in this pinned source; it is not an isolated l=0 curve.",
                       "Source truncates l at <=200 with relative incremental sigma stopping criterion 1e-4."]
    return result

def _execute(config, directory, authorization, timeout=60):
    """Internal primitive; CLI research entry is execute_plan, benchmark is fixed below."""
    validate(config)
    binary, build_record = build()
    d = Path(directory)
    d.mkdir(parents=True,exist_ok=False)
    deck, stdin, prompts = render(config)
    save(d/"config.json", config)
    save(d/"authorization.json", authorization)
    save(d/"build.json", build_record)
    (d/"ccfull.inp").write_text(deck)
    (d/"stdin.txt").write_text(stdin)
    save(d/"prompts.json",prompts)
    status = {"status":"running","started_at":now(),"config_sha256":digest(config),"source_sha256":SOURCE_SHA}
    save(d/"status.json",status)
    try:
        with (d/"stdout.log").open("w") as out, (d/"stderr.log").open("w") as err:
            proc = subprocess.run([str(binary)],cwd=d,input=stdin,text=True,stdout=out,stderr=err,timeout=timeout)
        status["returncode"] = proc.returncode
        require(proc.returncode==0, f"CCFULL exited {proc.returncode}; inspect logs")
        result = collect_run(d,config)
        result.update(config=config, config_sha256=digest(config), source_sha256=SOURCE_SHA,
                      binary_sha256=build_record["binary_sha256"], status="succeeded")
        save(d/"result.json",result)
        status["status"]="succeeded"
        status["warnings"]=result["warnings"]
    except (OSError,ValueError,subprocess.TimeoutExpired) as exc:
        status.update(status="failed",error=str(exc))
        raise ValidationError(f"Run failed at {d}: {exc}") from exc
    finally:
        status["finished_at"]=now()
        status["artifacts"]={p.name:sha(p) for p in d.iterdir() if p.is_file() and p.name!="status.json"}
        save(d/"status.json",status)
    return d

def execute_plan(plan, receipt, directory, timeout=60, keep_going=False):
    check_approval(plan,receipt)
    root = Path(directory)
    root.mkdir(parents=True,exist_ok=False)
    save(root/"plan.json",plan)
    save(root/"approval.json",receipt)
    summary={"status":"running","runs":[],"plan_sha256":digest(plan)}
    save(root/"batch.json",summary)
    try:
        for i,c in enumerate(plan["configs"]):
            try:
                d=_execute(c,root/f"{i:02d}",receipt,timeout)
                summary["runs"].append(str(d))
            except Exception as exc:
                summary.setdefault("failures",[]).append({"index":i,"id":c["id"],"directory":str(root/f"{i:02d}"),"error":str(exc)})
                if not keep_going:
                    raise
        if summary.get("failures"):
            raise ValidationError(f"{len(summary['failures'])} of {len(plan['configs'])} runs failed; all requested configurations were attempted; see {root/'batch.json'}")
        summary["status"]="succeeded"
    except Exception as exc:
        summary.update(status="failed",error=str(exc))
        raise
    finally:
        save(root/"batch.json",summary)
    return summary

def official_config():
    """Independent fixed configuration of the upstream nine-line example."""
    source="Hagino official ccfull.inp, pinned SHA256 "+INPUT_SHA
    return {"schema_version":1,"id":"official-16O-144Sm",
      "projectile":{"A":16,"Z":8,"ground_state":"0+","coupling_radius_fm":1.2,"source":source},
      "target":{"A":144,"Z":62,"ground_state":"0+","coupling_radius_fm":1.06,"source":source},
      "potential":{"id":"official","V0_MeV":105.1,"r0_fm":1.1,"a_fm":0.75,"radius_convention":"sum_A_cuberoot","source":source},
      "couplings":{"id":"target-3minus","projectile":{"type":"inert","modes":[]},
        "target":{"type":"vibrational","modes":[{"energy_MeV":1.81,"lambda":3,"beta_C":0.205,"beta_N":0.205,"phonons":1,"source":source}]},
        "target_mutual":"none"},
      "energy":{"frame":"cm","min_MeV":55.0,"max_MeV":70.0,"step_MeV":1.0},
      "radial":{"Rmax_fm":30.0,"dr_fm":0.05},"transfer":{"enabled":False},
      "assumptions":["Official software benchmark, not independently adopted current nuclear-structure data.",
                     "Projectile inert; harmonic target octupole mode; beta_N=beta_C; incoming-wave boundary condition."]}

def benchmark(directory):
    from .analysis import compare_numeric, plot_results
    verify_vendor()
    root=Path(directory)
    root.mkdir(parents=True,exist_ok=False)
    base=official_config()
    authorization={"scope":"software-validation-only","basis":"Fixed official benchmark requested through the benchmark command"}
    paths=[]
    for label,radial in (("official",{"Rmax_fm":30.0,"dr_fm":0.05}),
                         ("half-step",{"Rmax_fm":30.0,"dr_fm":0.025}),
                         ("larger-radius",{"Rmax_fm":40.0,"dr_fm":0.05})):
        c=copy.deepcopy(base)
        c["radial"]=radial
        paths.append(_execute(c,root/label,authorization))
    results=[read_json(p/"result.json") for p in paths]
    reference=parse_output(VENDOR/"OUTPUT")
    regression=compare_numeric(reference,results[0],rtol=1e-4,atol=5e-6,barrier_atol=0.01000001)
    convergence=[compare_numeric(results[0],r,rtol=0.01,atol=0,fields=("sigma_mb",),barrier_atol=None) for r in results[1:]]
    report={"passed":regression["passed"] and all(x["passed"] for x in convergence),
            "reference":{"path":str(VENDOR/"OUTPUT"),"sha256":OUTPUT_SHA},
            "regression":regression,"convergence":dict(zip(("dr_half","Rmax_plus10"),convergence)),
            "runs":[str(p) for p in paths],
            "limits":["Numerical verification only; does not select a physical model.",
                      "Channel-space convergence needs approved additional physical configurations.",
                      "Barrier-distribution energy-step convergence has not been established."]}
    save(root/"verification.json",report)
    plot_results(results,root/"figures",labels=["Official configuration","dr = 0.025 fm","Rmax = 40 fm"])
    lines=["# Official CCFULL validation","","PASS" if report["passed"] else "FAIL","",
           "Three runs share the official physics configuration; only numerical settings vary.","",
           "```json",json.dumps(report,indent=2),"```"]
    (root/"verification.md").write_text("\n".join(lines))
    require(report["passed"],f"Benchmark or convergence failed: {root/'verification.json'}")
    return report
