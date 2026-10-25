"""Evidence-assisted proposal preparation; never starts the solver."""
import copy
import re
from pathlib import Path
from urllib.parse import urlencode
from .model import require, validate, ValidationError
from .engine import save, make_plan, prepare

def source_queries(nucleus):
    m=re.fullmatch(r"(\d{1,3})([A-Z][a-z]?)",nucleus)
    require(m is not None,"Nucleus format must be e.g. 144Sm or 16O")
    mass,symbol=m.groups()
    return {"nucleus":nucleus,"status":"queries_only_not_retrieved_data",
      "nudat":"https://www.nndc.bnl.gov/nudat3/getdataset.jsp?"+urlencode({"nucleus":mass+symbol.upper(),"unc":"nds"}),
      "queries":[f'{nucleus} adopted levels B(E2) B(E3) ENSDF',
                 f'{nucleus} collective excitation deformation coupled channels fusion',
                 f'{nucleus} transfer fusion Q value'],
      "guidance":"Use connected browser/research tools to retrieve and inspect primary sources; record date, dataset revision, exact location, uncertainty, units and conflicts. A URL is not evidence that data was retrieved."}

def propose(catalog, directory):
    required={"schema_version","system","base_config","sources","levels","potential_candidates",
              "coupling_candidates","assessment","gaps","transfer"}
    require(set(catalog)==required and catalog["schema_version"]==1,"Unsupported evidence catalog schema")
    require(catalog["sources"] and catalog["assessment"],"Sources and an evidence-based assessment required")
    require(1<=len(catalog["potential_candidates"])<=3,"At most three sourced potential candidates")
    root=Path(directory)
    root.mkdir(parents=True,exist_ok=False)
    save(root/"evidence.json",catalog)
    configs=[]
    blocked=[]
    # Every parameter is explicitly supplied by the curated catalog, never guessed.
    for potential in catalog["potential_candidates"]:
        for choice in catalog["coupling_candidates"]:
            c=copy.deepcopy(catalog["base_config"])
            c["potential"]=copy.deepcopy(potential)
            c["couplings"]=copy.deepcopy(choice["couplings"])
            c["id"]=f"{potential['id']}__{choice['id']}"
            c["assumptions"]=c["assumptions"]+[choice["assumption"]]
            try:
                info=validate(c)
            except (ValidationError,KeyError,TypeError) as exc:
                blocked.append({"id":c["id"],"reason":str(exc)})
                continue
            configs.append(c)
    lines=[f"# {catalog['system']} 耦合与势参数建议单","","状态：待讨论，未运行研究计算。","",
           "## 物理判断","",*["- "+s for s in catalog["assessment"]],"",
           "## 核结构资料","","| 核 | 态 | 资料中的激发能 (MeV) | 候选计算值与来源 | 缺口 |",
           "|---|---|---:|---|---|"]
    for level in catalog["levels"]:
        energy=str(level["energy_MeV"])
        if level.get("energy_uncertainty_MeV") is not None:
            energy+=f" ± {level['energy_uncertainty_MeV']}"
        lines.append(f"| {level['nucleus']} | {level['spin_parity']} | {energy} ({level['energy_source']}) | {level['model_use']} | {level['gap']} |")
    lines += ["","### 已检索的跃迁强度","",
              "| 核与初态 | B(Eλ) | 方向 | 来源 |","|---|---|---|---|"]
    for level in catalog["levels"]:
        if level["B_value"] is not None:
            lines.append(f"| {level['nucleus']} {level['spin_parity']} | {level['B_value']} ± {level['B_uncertainty']} {level['B_unit']} | {level['B_direction']} | {level['B_source']} |")
    lines += ["","W.u.数值仅作为证据记录；beta-wu 工具可在明确方向与半径约定后换算，未自动替换候选模型参数。",
              "","能级数据与模型参数分别保留；空缺不是零。所有参数尚待讨论确认。","",
              "## 势参数候选","","R = r0 × (Ap^(1/3) + At^(1/3))；各组均为实 Woods–Saxon 势。",
              "","| ID | V0 (MeV) | r0 (fm) | a (fm) | 来源 |","|---|---:|---:|---:|---|"]
    for p in catalog["potential_candidates"]:
        lines.append(f"| {p['id']} | {p['V0_MeV']} | {p['r0_fm']} | {p['a_fm']} | {p['source']} |")
    lines += ["","## 耦合方案","","| ID | 方案与理由 | 模型假设 |","|---|---|---|"]
    for choice in catalog["coupling_candidates"]:
        lines.append(f"| {choice['id']} | {choice['rationale']} | {choice['assumption']} |")
    lines += ["","## 建议的讨论顺序","",
       "先确认能区和参数约定，再选择耦合与势参数组合。先固定一组势比较无耦合与逐级增加的耦合；再固定耦合比较势参数。互激发与截断作为单独的敏感性问题讨论。这里的候选不代表已获运行授权。",
       "","候选全集只是可审阅菜单。选择较小子集时重新 prepare，并对生成的新 plan.json 记录确认。",
       "","## 缺失信息与冲突","",*["- "+s for s in catalog["gaps"]],
       "","## 转移耦合提示","",catalog["transfer"],
       "","## 已检索来源",""]
    for source in catalog["sources"]:
        lines += [f"- **{source['id']}** [{source['title']}]({source['url']})。位置：{source['location']}；访问：{source['accessed']}。{source['limitations']}"]
    if blocked:
        lines += ["","## 暂不可运行的组合","",*["- "+b["id"]+": "+b["reason"] for b in blocked]]
    if configs:
        prepare(make_plan(configs,"候选耦合与势参数组合，先讨论并选择，再确认运行。"),root/"candidates")
    save(root/"proposal-status.json",{"status":"proposed","runnable_candidates":len(configs),"blocked":blocked})
    (root/"proposal.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    return {"report":str(root/"proposal.md"),"runnable_candidates":len(configs),"blocked":blocked,"status":"awaiting_discussion"}
