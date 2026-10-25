"""Numerical comparisons and transparent diagnostic figures."""
import copy
import csv
import math
from pathlib import Path
from .model import require, digest, read_json

def compare_numeric(reference, actual, rtol, atol, fields=("sigma_mb","mean_l_hbar"), barrier_atol=0.01000001):
    a,b=reference["table"],actual["table"]
    require(len(a)==len(b),"Cannot compare different energy grids")
    failures=[]
    maximum={field:0.0 for field in fields}
    for x,y in zip(a,b):
        require(math.isclose(x["E_cm_MeV"],y["E_cm_MeV"],abs_tol=5.1e-6,rel_tol=0),"Energy grids differ; no implicit interpolation")
        for field in fields:
            ref,val=x[field],y[field]
            require(math.isfinite(ref) and math.isfinite(val),"Non-finite comparison")
            delta=abs(val-ref)
            relative=delta/abs(ref) if ref else (0.0 if val==0 else None)
            if relative is not None:
                maximum[field]=max(maximum[field],relative)
            if not math.isclose(ref,val,rel_tol=rtol,abs_tol=atol):
                failures.append({"energy_MeV":x["E_cm_MeV"],"field":field,"reference":ref,"actual":val,"relative_error":relative})
    if barrier_atol is not None:
        for k,v in reference["barrier"].items():
            if abs(v-actual["barrier"][k])>barrier_atol:
                failures.append({"field":k,"reference":v,"actual":actual["barrier"][k]})
    return {"passed":not failures,"rtol":rtol,"atol":atol,"max_relative_error":maximum,"failures":failures}

def barrier_distribution(table):
    require(len(table)>=3,"At least three energy points required")
    out=[]
    for left,center,right in zip(table,table[1:],table[2:]):
        h1=center["E_cm_MeV"]-left["E_cm_MeV"]
        h2=right["E_cm_MeV"]-center["E_cm_MeV"]
        require(h1>0 and math.isclose(h1,h2,rel_tol=1e-8,abs_tol=1e-10),"Equal positive energy spacing required")
        value=(right["E_cm_MeV"]*right["sigma_mb"]-2*center["E_cm_MeV"]*center["sigma_mb"]+left["E_cm_MeV"]*left["sigma_mb"])/h1**2
        out.append({"E_cm_MeV":center["E_cm_MeV"],"D_mb_per_MeV":value})
    return out

def plot_results(results, directory, labels=None, *, allow_reviewed=False):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from .engine import save
    out=Path(directory)
    out.mkdir(parents=True,exist_ok=False)
    require(1<=len(results)<=30,"Need 1..30 results")
    for r in results:
        require(r.get("status")=="succeeded" or (allow_reviewed and r.get("status")=="accepted_with_warning"
                and r.get("review",{}).get("scope")=="cross_section_comparison_only"),
                "Only successful or explicitly reviewed results may be plotted")
    labels=[r["config"]["id"] for r in results] if labels is None else labels
    require(len(labels)==len(results),"Each result needs exactly one label")
    colors=["#0072B2","#D55E00","#009E73","#CC79A7","#000000","#6B5B00"]
    styles=["-","--",":","-."]
    markers=["o","s","^","D"]
    metadata={"transformation":"D(E): equal-spacing three-point second derivative of E*sigma; interior points only; no smoothing or interpolation",
              "log_policy":"Log10 sigma only if all points positive; otherwise linear, no points omitted",
              "data":[],"matplotlib_version":matplotlib.__version__}
    with plt.rc_context({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"savefig.facecolor":"white"}):
        reviewed=any(r.get("status")=="accepted_with_warning" for r in results)
        fig,axes=plt.subplots(1,1 if reviewed else 3,figsize=(7,5) if reviewed else (13.5,4.5),layout="constrained",squeeze=False)
        axes=axes.ravel()
        if reviewed:
            fig.suptitle("Cross-section comparison; C0 accepted with numerical warning",fontsize=10)
            metadata['review_scope']='Cross sections only; reviewed C0 mean-l and barrier distributions not assessed.'
            metadata['transformation']='Cross sections as output; no probability clipping and no D(E) calculation.'
        all_positive=all(row["sigma_mb"]>0 for r in results for row in r["table"])
        for i,(r,label) in enumerate(zip(results,labels)):
            table=r["table"]
            x=[z["E_cm_MeV"] for z in table]
            style={"color":colors[i%len(colors)],"linestyle":styles[i%len(styles)],"marker":markers[i%len(markers)],"markersize":3}
            axes[0].plot(x,[z["sigma_mb"] for z in table],label=label,**style)
            if not reviewed:
                axes[1].plot(x,[z["mean_l_hbar"] for z in table],label=label,**style)
            dist=barrier_distribution(table) if len(table)>=3 and not reviewed else []
            if dist:
                axes[2].plot([z["E_cm_MeV"] for z in dist],[z["D_mb_per_MeV"] for z in dist],label=label,**style)
            with (out/f"series-{i:02d}.csv").open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=["E_cm_MeV","sigma_mb"] if reviewed else ["E_cm_MeV","sigma_mb","mean_l_hbar"],extrasaction="ignore")
                writer.writeheader()
                writer.writerows(table)
            with (out/f"barrier-{i:02d}.csv").open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=["E_cm_MeV","D_mb_per_MeV"])
                writer.writeheader()
                writer.writerows(dist)
            metadata["data"].append({"label":label,"config_sha256":r.get("config_sha256"),"energy_step_MeV":r["config"]["energy"]["step_MeV"],"warnings":r.get("warnings",[])})
        if all_positive:
            axes[0].set_yscale("log",base=10)
        axes[0].set_ylabel("Fusion cross section (mb)" + ("; log10 scale" if all_positive else ""))
        if not reviewed:
            axes[1].set_ylabel("Mean angular momentum / hbar")
            axes[2].set_ylabel("D(E) (mb / MeV)")
            axes[2].axhline(0,color="#777777",lw=0.6)
        for ax,title in zip(axes,["Fusion excitation function","Mean angular momentum","Barrier distribution (diagnostic)"]):
            ax.set_xlabel("Center-of-mass energy (MeV)")
            ax.set_title(title,fontsize=11)
            ax.grid(alpha=0.2)
        axes[0].legend(fontsize=7)
        fig.savefig(out/"comparison.png",dpi=180)
        fig.savefig(out/"comparison.svg")
        plt.close(fig)
    save(out/"figure.json",metadata)
    (out/"description.txt").write_text(("Fusion cross sections only; C0 is explicitly accepted with a numerical warning. Original strict failure remains recorded. " if reviewed else "Three panels show calculated fusion cross sections, mean angular momenta and the unsmoothed finite-difference barrier distribution. ")+"Data and transformation metadata are provided alongside the figure. No experimental points or experimental uncertainties are displayed.\n")
    return {"figure":str(out/"comparison.png"),"metadata":str(out/"figure.json")}

def compare_runs(paths, directory, dimension):
    from .engine import save, sha
    require(dimension in ("couplings","potential"),"Choose couplings or potential comparison")
    results=[]
    for path in paths:
        d=Path(path)
        if (d/"review.json").is_file():
            from .review import load_reviewed
            results.append(load_reviewed(d))
            continue
        status=read_json(d/"status.json")
        require(status["status"]=="succeeded","Failed run cannot enter comparison")
        require(sha(d/"result.json")==status["artifacts"]["result.json"],"Result changed after run")
        results.append(read_json(d/"result.json"))
    require(len(results)>=2,"Comparison requires >=2 runs")
    def fixed(c):
        c=copy.deepcopy(c)
        c.pop("id")
        c.pop("assumptions")
        c.pop(dimension)
        return c
    fixed_hash=digest(fixed(results[0]["config"]))
    require(all(digest(fixed(r["config"]))==fixed_hash for r in results),"Comparison changes more than the declared factor")
    require(len({r["source_sha256"] for r in results})==1 and len({r["binary_sha256"] for r in results})==1,"Mixed solver builds")
    for r in results[1:]:
        compare_numeric(results[0],r,rtol=1,atol=0,barrier_atol=None)
    out=Path(directory)
    out.mkdir(parents=True,exist_ok=False)
    differences=[]
    for r in results[1:]:
        rows=[]
        for x,y in zip(results[0]["table"],r["table"]):
            rows.append({"E_cm_MeV":x["E_cm_MeV"],"delta_sigma_mb":y["sigma_mb"]-x["sigma_mb"],
                         "relative_delta":(y["sigma_mb"]-x["sigma_mb"])/x["sigma_mb"] if x["sigma_mb"] else None})
        differences.append({"id":r["config"]["id"],"rows":rows})
    save(out/"comparison.json",{"dimension":dimension,"reference":results[0]["config"]["id"],
                              "runs":[str(p) for p in paths],"differences":differences,
                              "reviews":[r.get("review") for r in results if r.get("review")],
                              "note":"Sensitivity comparison; no goodness-of-fit or physical mechanism selection."})
    return plot_results(results,out/"figures",allow_reviewed=True)
