"""Command line interface. Research run always requires an exact-plan receipt."""
import argparse
import json
import sys
from pathlib import Path
from .model import ValidationError, read_json, validate, beta_from_transition
from .engine import ROOT, build, make_plan, prepare, approve, execute_plan, benchmark, save
from .analysis import compare_runs

def main(argv=None):
    p=argparse.ArgumentParser(prog="ccfull-agent")
    subs=p.add_subparsers(dest="command",required=True)
    subs.add_parser("build")
    b=subs.add_parser("benchmark")
    b.add_argument("--out",required=True)
    v=subs.add_parser("validate")
    v.add_argument("config")
    q=subs.add_parser("prepare")
    q.add_argument("configs",nargs="+")
    q.add_argument("--purpose",required=True)
    q.add_argument("--out",required=True)
    a=subs.add_parser("approve",help="Only after actual human confirmation of this exact plan")
    a.add_argument("plan")
    a.add_argument("--statement",required=True)
    a.add_argument("--actor",required=True)
    a.add_argument("--out",required=True)
    r=subs.add_parser("run")
    r.add_argument("plan")
    r.add_argument("--approval",required=True)
    r.add_argument("--out",required=True)
    r.add_argument("--timeout",type=float,default=60)
    r.add_argument("--keep-going",action="store_true",help="Attempt remaining approved runs after a failure; batch still fails")
    c=subs.add_parser("compare")
    c.add_argument("runs",nargs="+")
    c.add_argument("--dimension",choices=["couplings","potential"],required=True)
    c.add_argument("--out",required=True)
    review=subs.add_parser("review-c0",help="Record actual user acceptance of small C0 probability overshoot; no solver")
    review.add_argument("runs",nargs="+")
    review.add_argument("--statement",required=True)
    review.add_argument("--actor",required=True)
    review.add_argument("--out",required=True)
    e=subs.add_parser("propose")
    e.add_argument("--catalog",required=True,help="Path to a user-supplied evidence catalog")
    e.add_argument("--out",required=True)
    n=subs.add_parser("sources")
    n.add_argument("--nucleus",required=True,help="e.g. 144Sm, 16O; creates URLs, not fetched facts")
    n.add_argument("--out",required=True)
    t=subs.add_parser("beta")
    t.add_argument("--B",type=float,required=True)
    t.add_argument("--unit",choices=["e2_fm2lambda","e2_b_lambda"],required=True)
    t.add_argument("--A",type=int,required=True)
    t.add_argument("--Z",type=int,required=True)
    t.add_argument("--lambda",dest="lam",type=int,required=True)
    t.add_argument("--radius",type=float,required=True)
    t.add_argument("--direction",choices=["up","down"],required=True)
    t.add_argument("--uncertainty",type=float)
    w=subs.add_parser("beta-wu",help="Explicit single-particle Weisskopf convention; magnitude only")
    for name,typ in [("B",float),("A",int),("Z",int),("radius",float)]:
        w.add_argument("--"+name,type=typ,required=True)
    w.add_argument("--lambda",dest="lam",type=int,required=True)
    w.add_argument("--direction",choices=["up","down"],required=True)
    w.add_argument("--uncertainty",type=float)
    f=subs.add_parser("fetch-structure",help="Retrieve NuDat snapshot and parse; no solver")
    f.add_argument("--nucleus",required=True)
    f.add_argument("--out",required=True)
    i=subs.add_parser("import-structure",help="Parse a hashed local NuDat snapshot")
    i.add_argument("snapshot")
    i.add_argument("--nucleus",required=True)
    i.add_argument("--manifest",required=True)
    i.add_argument("--out",required=True)
    s=subs.add_parser("assess-structure",help="Explainable coupling triage, not an approved model")
    s.add_argument("structure")
    s.add_argument("--radius",type=float,required=True)
    s.add_argument("--max-energy",type=float,default=7.0)
    s.add_argument("--min-wu",type=float,default=5.0)
    s.add_argument("--out",required=True)
    i=subs.add_parser("import-experiment")
    i.add_argument("csv")
    i.add_argument("--metadata",required=True)
    i.add_argument("--out",required=True)
    for name in ("compare-experiment","fit-grid"):
        c=subs.add_parser(name)
        c.add_argument("runs",nargs="+")
        c.add_argument("--data",required=True)
        c.add_argument("--interpolation",choices=["none","linear"],default="none")
        c.add_argument("--out",required=True)
        if name=="fit-grid":
            c.add_argument("--scan",required=True)
            c.add_argument("--delta-chi2",type=float,default=1.0)
    s=subs.add_parser("prepare-scan",help="Prepare a bounded parameter/radial scan; never runs it")
    s.add_argument("config")
    s.add_argument("--axes",required=True,help="JSON object: parameter path to explicit values")
    s.add_argument("--purpose",required=True)
    s.add_argument("--out",required=True)
    d=subs.add_parser("diagnose")
    d.add_argument("runs",nargs="+")
    d.add_argument("--out",required=True)
    c=subs.add_parser("convergence")
    c.add_argument("runs",nargs="+")
    c.add_argument("--rtol",type=float,default=.01)
    c.add_argument("--out",required=True)
    args=p.parse_args(argv)
    try:
        if args.command=="build":
            binary,record=build()
            result={"binary":str(binary),**record}
        elif args.command=="benchmark":
            result=benchmark(args.out)
        elif args.command=="validate":
            result=validate(read_json(args.config))
        elif args.command=="prepare":
            result=prepare(make_plan([read_json(x) for x in args.configs],args.purpose),args.out)
        elif args.command=="approve":
            if Path(args.out).exists():
                raise ValidationError("Approval file already exists; use a new destination")
            result=approve(read_json(args.plan),args.statement,args.actor)
            save(args.out,result)
        elif args.command=="run":
            if not 0<args.timeout<=600:
                raise ValidationError("Timeout must be in (0,600] seconds")
            result=execute_plan(read_json(args.plan),read_json(args.approval),args.out,args.timeout,keep_going=args.keep_going)
        elif args.command=="compare":
            result=compare_runs(args.runs,args.out,args.dimension)
        elif args.command=="review-c0":
            from .review import review_c0
            result=review_c0(args.runs,args.statement,args.actor,args.out)
        elif args.command=="propose":
            from .evidence import propose
            result=propose(read_json(args.catalog),args.out)
        elif args.command=="sources":
            from .evidence import source_queries
            result=source_queries(args.nucleus)
            if Path(args.out).exists():
                raise ValidationError("Source query file exists")
            save(args.out,result)
        elif args.command=="beta":
            result=beta_from_transition(args.B,args.unit,args.A,args.Z,args.lam,args.radius,args.direction,args.uncertainty)
        elif args.command=="beta-wu":
            from .structure import beta_from_wu
            result=beta_from_wu(args.B,args.A,args.Z,args.lam,args.radius,args.direction,args.uncertainty)
        elif args.command=="fetch-structure":
            from .structure import fetch_structure
            result=fetch_structure(args.nucleus,args.out)
        elif args.command=="import-structure":
            from .structure import import_structure
            result=import_structure(args.snapshot,args.nucleus,read_json(args.manifest),args.out)
        elif args.command=="assess-structure":
            from .structure import assess_structure
            from .experiment import new_json
            result=new_json(args.out,assess_structure(read_json(args.structure),args.radius,args.max_energy,args.min_wu))
        elif args.command=="import-experiment":
            from .experiment import import_experiment
            result=import_experiment(args.csv,read_json(args.metadata),args.out)
        elif args.command=="compare-experiment":
            from .experiment import compare_experiment
            result=compare_experiment(args.runs,read_json(args.data),args.out,args.interpolation)
        elif args.command=="prepare-scan":
            from .planning import prepare_scan
            result=prepare_scan(read_json(args.config),read_json(args.axes),args.out,args.purpose)
        elif args.command=="fit-grid":
            from .planning import fit_grid
            result=fit_grid(args.runs,read_json(args.data),read_json(args.scan),args.out,args.interpolation,args.delta_chi2)
        elif args.command=="diagnose":
            from .planning import diagnose_run
            from .experiment import new_json
            result=new_json(args.out,{'diagnostics':[diagnose_run(p) for p in args.runs]})
        elif args.command=="convergence":
            from .planning import convergence_report
            result=convergence_report(args.runs,args.out,args.rtol)
        print(json.dumps(result,ensure_ascii=False,indent=2,allow_nan=False))
        return 0
    except (ValidationError,OSError,KeyError,TypeError,ValueError,IndexError) as exc:
        print(json.dumps({"status":"failed","error":str(exc)},ensure_ascii=False),file=sys.stderr)
        return 2

if __name__=="__main__":
    sys.exit(main())
