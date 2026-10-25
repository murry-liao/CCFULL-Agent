"""Strict, unit-explicit configurations for the pinned 2015 CCFULL source."""
import copy
import hashlib
import json
import math

class ValidationError(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise ValidationError(message)

def keys(obj, expected, path):
    require(isinstance(obj, dict), f"{path}: expected object")
    require(set(obj) == set(expected), f"{path}: expected keys {sorted(expected)}; got {sorted(obj)}")

def number(value, path, low=None, high=None, integer=False):
    require(type(value) in (int, float) and math.isfinite(value), f"{path}: finite number required")
    if integer:
        require(type(value) is int, f"{path}: integer required")
    if low is not None:
        require(value >= low, f"{path}: must be >= {low}")
    if high is not None:
        require(value <= high, f"{path}: must be <= {high}")
    return value

def text(value, path):
    require(isinstance(value, str) and bool(value.strip()), f"{path}: nonempty text required")

def digest(obj):
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()).hexdigest()

def read_json(path):
    def invalid(s):
        raise ValidationError(f"Non-finite JSON constant: {s}")
    def pairs(items):
        result = {}
        for k, v in items:
            require(k not in result, f"Duplicate JSON key: {k}")
            result[k] = v
        return result
    with open(path, encoding="utf-8") as stream:
        return json.load(stream, parse_constant=invalid, object_pairs_hook=pairs)

def validate(config):
    keys(config, ["schema_version", "id", "projectile", "target", "potential", "couplings",
                  "energy", "radial", "transfer", "assumptions"], "config")
    require(config["schema_version"] == 1, "Only schema_version 1 is supported")
    text(config["id"], "id")
    require(isinstance(config["assumptions"], list) and config["assumptions"], "Explicit assumptions required")
    for item in config["assumptions"]:
        text(item, "assumptions")
    for side in ("projectile", "target"):
        n = config[side]
        keys(n, ["A", "Z", "ground_state", "coupling_radius_fm", "source"], side)
        number(n["A"], side+".A", 1, 300, True)
        number(n["Z"], side+".Z", 1, min(103, n["A"]), True)
        require(n["ground_state"] == "0+", f"{side}: v1 supports 0+ ground states only")
        number(n["coupling_radius_fm"], side+".coupling_radius_fm", 0.1, 3)
        text(n["source"], side+".source")
    p = config["potential"]
    keys(p, ["id", "V0_MeV", "r0_fm", "a_fm", "radius_convention", "source"], "potential")
    text(p["id"], "potential.id")
    text(p["source"], "potential.source")
    require(p["radius_convention"] == "sum_A_cuberoot", "Potential radius must be r0*(Ap^(1/3)+At^(1/3))")
    for field, low, high in (("V0_MeV", 0.001, 3000), ("r0_fm", 0.1, 3), ("a_fm", 0.01, 5)):
        number(p[field], "potential."+field, low, high)
    couplings = config["couplings"]
    keys(couplings, ["id", "projectile", "target", "target_mutual"], "couplings")
    text(couplings["id"], "couplings.id")
    require(couplings["target_mutual"] in ("none", "all"), "target_mutual must be none or all")
    counts = {}
    warnings = []
    for side, max_modes in (("projectile", 1), ("target", 2)):
        s = couplings[side]
        keys(s, ["type", "modes"], "couplings."+side)
        require(s["type"] in ("inert", "vibrational", "rotational"), "Unsupported intrinsic motion")
        require(isinstance(s["modes"], list), "modes must be an array")
        if s["type"] == "inert":
            require(not s["modes"], "Inert nucleus cannot contain active modes")
            counts[side] = [0]
            continue
        require(1 <= len(s["modes"]) <= (1 if s["type"] == "rotational" else max_modes),
                f"{side}: unsupported number/mix of modes")
        counts[side] = []
        multipoles = []
        for mode in s["modes"]:
            if s["type"] == "vibrational":
                keys(mode, ["energy_MeV", "lambda", "beta_C", "beta_N", "phonons", "source"], side+".mode")
                number(mode["energy_MeV"], "energy_MeV", 0.000001, 100)
                number(mode["lambda"], "lambda", 2, 6, True)
                for k in ("beta_C", "beta_N"):
                    number(mode[k], k, -2, 2)
                number(mode["phonons"], "phonons", 1, 29, True)
                counts[side].append(mode["phonons"])
                multipoles.append(mode["lambda"])
                if mode["phonons"] > 1:
                    warnings.append(f"{side}: harmonic multi-phonon truncation is an assumption, not measured spectrum")
            else:
                keys(mode, ["E2_MeV", "beta2", "beta4", "levels", "source"], side+".rotor")
                number(mode["E2_MeV"], "E2_MeV", 0.000001, 100)
                number(mode["beta2"], "beta2", -2, 2)
                number(mode["beta4"], "beta4", -2, 2)
                number(mode["levels"], "levels", 1, 29, True)
                counts[side].append(mode["levels"])
                warnings.append(f"{side}: rigid K=0 rotor; beta signs must be supported, not inferred from B(E2)")
            text(mode["source"], side+".mode.source")
        require(len(set(multipoles)) == len(multipoles), "Two modes of the same multipolarity need explicit model support; v1 rejects them")
    t = counts["target"]
    if len(t) == 1:
        require(couplings["target_mutual"] == "none", "target_mutual only applies to two target modes")
        nt = t[0] + 1
    elif couplings["target_mutual"] == "all":
        nt = (t[0]+1)*(t[1]+1)
    else:
        nt = 1 + sum(t)
    channels = nt*(counts["projectile"][0]+1)
    require(channels <= 30, f"{channels} channels exceeds pinned source limit 30")
    require(config["transfer"] == {"enabled": False}, "v1 only flags transfer evidence; execution is disabled")
    e = config["energy"]
    keys(e, ["frame", "min_MeV", "max_MeV", "step_MeV"], "energy")
    require(e["frame"] == "cm", "Only explicit center-of-mass energy accepted; convert lab energy first")
    for k in ("min_MeV", "max_MeV", "step_MeV"):
        number(e[k], "energy."+k, 0.000001, 10000)
    require(e["max_MeV"] >= e["min_MeV"], "Reversed energy interval")
    steps = (e["max_MeV"]-e["min_MeV"])/e["step_MeV"]
    require(abs(steps-round(steps)) < 1e-7 and steps <= 10000, "Energy grid must hit its endpoint and have <=10001 points")
    r = config["radial"]
    keys(r, ["Rmax_fm", "dr_fm"], "radial")
    number(r["Rmax_fm"], "Rmax_fm", 5, 100)
    number(r["dr_fm"], "dr_fm", 0.0001, 0.5)
    require(math.ceil(r["Rmax_fm"]/r["dr_fm"])+2 <= 5000, "Conservative radial grid exceeds source limit 5000")
    return {"channels": channels, "points": round(steps)+1, "warnings": warnings}

def energy_grid(c):
    n = validate(c)["points"]
    e = c["energy"]
    return [e["min_MeV"]+i*e["step_MeV"] for i in range(n)]

def render(c):
    """Return deck, exact finite stdin sequence and its human-readable rationale."""
    validate(c)
    co = c["couplings"]
    p, t = c["projectile"], c["target"]
    flags = {"inert": -1, "vibrational": 0, "rotational": 1}
    def mode_line(s, index=0):
        if index >= len(s["modes"]):
            return [0, 0, 2, 0]
        m = s["modes"][index]
        if s["type"] == "vibrational":
            return [m["energy_MeV"], m["beta_C"], m["lambda"], m["phonons"]]
        return [m["E2_MeV"], m["beta2"], m["beta4"], m["levels"]]
    v, e, r = c["potential"], c["energy"], c["radial"]
    lines = [[p["A"], p["Z"], t["A"], t["Z"]],
             [p["coupling_radius_fm"], flags[co["projectile"]["type"]], t["coupling_radius_fm"], flags[co["target"]["type"]]],
             mode_line(co["target"]), mode_line(co["target"], 1), mode_line(co["projectile"]),
             [0, 0, 0], [v["V0_MeV"], v["r0_fm"], v["a_fm"]],
             [e["min_MeV"], e["max_MeV"], e["step_MeV"]], [r["Rmax_fm"], r["dr_fm"]]]
    answers = [("Modified Woods-Saxon?", "n")]
    def beta_answer(m, label):
        different = m["beta_N"] != m["beta_C"]
        answers.append((label+": separate beta_N?", "y" if different else "n"))
        if different:
            answers.append((label+": beta_N", format(m["beta_N"], ".15g")))
    for i, m in enumerate(co["target"]["modes"]):
        if co["target"]["type"] == "vibrational":
            beta_answer(m, f"target mode {i+1}")
        if i == 1:
            include = co["target_mutual"] == "all"
            answers.append(("Include target mutual excitations?", "y" if include else "n"))
            if include:
                answers.append(("Include all target mutual excitations?", "y"))
    for m in co["projectile"]["modes"]:
        if co["projectile"]["type"] == "vibrational":
            beta_answer(m, "projectile")
    for side in ("target", "projectile"):
        s = co[side]
        if s["type"] == "vibrational":
            for i, m in enumerate(s["modes"]):
                if m["phonons"] > 1:
                    answers.append((f"Anharmonic {side} mode {i+1}?", "n"))
    for side in ("target", "projectile"):
        if co[side]["type"] == "rotational":
            answers.append((f"Generalised {side} E2?", "n"))
    deck = "\n".join(",".join(format(x, ".15g") for x in row) for row in lines)+"\n"
    stdin = "\n".join(a for _, a in answers)+"\n"
    return deck, stdin, [{"prompt": q, "answer": a} for q, a in answers]

def beta_from_transition(B, unit, A, Z, multipolarity, radius_fm, direction, uncertainty=None):
    """Magnitude only: 0+ <-> lambda; no nuclear beta or sign inference."""
    number(B, "B", 0.000000000001)
    number(A, "A", 1, 300, True)
    number(Z, "Z", 1, A, True)
    number(multipolarity, "lambda", 2, 6, True)
    number(radius_fm, "radius_fm", 0.01)
    require(unit in ("e2_fm2lambda", "e2_b_lambda"), "Convert W.u. with a separately documented convention")
    require(direction in ("up", "down"), "Specify 0+->lambda (up) or lambda->0+ (down)")
    factor = (100**multipolarity if unit == "e2_b_lambda" else 1)*(2*multipolarity+1 if direction == "down" else 1)
    R = radius_fm*A**(1/3)
    beta = 4*math.pi/(3*Z*R**multipolarity)*math.sqrt(B*factor)
    if uncertainty is not None:
        number(uncertainty, "uncertainty", 0)
    return {"abs_beta_C": beta, "uncertainty": None if uncertainty is None else beta*uncertainty/(2*B),
            "B_up_e2_fm2lambda": B*factor, "R_fm": R,
            "assumption": "0+ ground state, collective surface model; radius fixed; sign and beta_N undetermined"}
