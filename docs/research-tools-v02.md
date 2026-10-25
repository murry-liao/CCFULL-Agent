# 使用指南

所有命令从项目根目录通过 `.venv/bin/python scripts/ccfull.py` 调用。使用各子命令的 `--help` 查看参数；以下路径为用户当前任务的输入/输出路径，不表示随包附带研究资料。

## 核结构资料

- `sources --nucleus NUCLEUS --out queries.json` 生成查询入口，不代表已检索到数据。
- `fetch-structure --nucleus NUCLEUS --out NEW_DIRECTORY` 获取NuDat页面并保存原文、下载时间、哈希和解析记录。
- `import-structure SNAPSHOT --nucleus NUCLEUS --manifest MANIFEST --out structure.json` 从已有页面读取数据。manifest是数组，每条记录包含file、url、retrieved_at、sha256；文件名和哈希必须匹配。
- `assess-structure structure.json --radius RADIUS --out assessment.json` 筛选候选态。`NUCLEUS` 使用质量数加元素符号的格式；`RADIUS` 为耦合半径系数，单位fm。

只自动解释支持的 adopted HTML 表和普通数值/对称误差。保留评价截止日期、原文、行号与未解决项；通过主表gamma分支匹配末态，不以激发能减gamma能量直接判断基态跃迁。上限、非对称误差和歧义值不自动变成确定参数。

默认筛选7 MeV以下的2⁺/3⁻态；明确的基态退激强度达到5 W.u.时列为优先候选。`--max-energy`、`--min-wu`可调整。这些是筛查规则，不是必须纳入耦合的判据，R4/2也不单独决定转子模型。

`beta` 接受显式物理单位的B值。`beta-wu`采用单粒子Weisskopf单位：

```text
B_W(Eλ) = (1/4π) [3/(λ+3)]² (1.2 A^(1/3))^(2λ) e² fm^(2λ)
B_up = (2λ+1) B_down       # 仅限 λ → 0 的退激
abs(βC) = 4π sqrt(B_up/e²) / (3 Z R^λ)
R = rc A^(1/3)
```

方向必须通过`--direction up|down`指定，耦合半径系数通过`--radius`指定。Weisskopf半径系数与耦合半径约定分开，方向转换只进行一次。误差仅传播B的对称误差并固定半径，不给出βN或符号。来源：[NNDC B(E2) Help](https://www.nndc.bnl.gov/be2/help.html)。

## 证据目录与候选方案

`propose --catalog CATALOG --out NEW_DIRECTORY` 要求显式提供当前研究任务的JSON目录，不自动载入特定体系。顶层字段为：

| 字段 | 内容 |
|---|---|
| schema_version | 1 |
| system | 体系名称 |
| base_config | 完整基础配置，字段结构参照官方输入配置；物理数值必须来自当前研究任务 |
| sources | 来源数组，每项含id、title、url、location、accessed、limitations |
| levels | 能级记录数组，见下文 |
| potential_candidates | 最多3组有来源的实势候选，与config中的potential结构相同 |
| coupling_candidates | 每项含id、couplings、rationale、assumption；couplings与配置结构相同 |
| assessment | 非空的判断说明字符串数组 |
| gaps | 缺失信息与冲突的字符串数组 |
| transfer | 转移通道证据/缺口说明字符串，不是运行参数 |

每项levels包含nucleus、spin_parity、energy_MeV、energy_source、model_use、gap，以及B_value、B_uncertainty、B_unit、B_direction、B_source。energy_uncertainty_MeV为可选字段；未知的数值可使用null，不能填零或凭空推断。

配置采用显式单位：能量MeV、长度fm、截面mb、无量纲β，能量坐标为cm。势半径为 `r0*(Ap^(1/3)+At^(1/3))`；各核耦合半径独立指定。振动模式需给energy_MeV、lambda、beta_C、beta_N、phonons及source。转子模式需给E2_MeV、beta2、beta4、levels及source。其他完整字段见`examples/official.json`和`model.py`校验定义。

候选生成不执行计算。最终选择必须包括βN/βC、形变符号、模型类型、截断、互激发、势参数和数值设置的依据。

## 研究运行

1. 用`prepare CONFIGS --purpose PURPOSE --out REVIEW_DIRECTORY`生成精确计划和可审阅输入。
2. 与研究者讨论具体配置；确认后用`approve PLAN --actor ACTOR --statement ACTUAL_USER_WORDS --out RECEIPT`记录实际意见。不能把这里的参数占位词写成真实确认。
3. 用`run PLAN --approval RECEIPT --out RUN_DIRECTORY`运行。`--keep-going`可继续尝试同一已确认批次的其余配置；失败状态仍保留。

任何输入或计划改变都会使旧确认不匹配。不得复用已完成运行目录。确认记录是审计证据，不是身份认证系统。

## 实验截面导入

`import-experiment CSV --metadata METADATA --out DATASET` 将明确单位的融合截面规范化。metadata字段：

| 字段 | 约定 |
|---|---|
| schema_version、id、source | 版本1、数据集标识、原始来源与具体位置 |
| kind | measured或synthetic，不能把合成数据标成测量 |
| system | projectile和target，各包含A与Z |
| observable | fusion；其他观测量需要先解决物理换算 |
| energy_unit | MeV或keV |
| energy_frame | cm或显式的lab_mass_number近似 |
| sigma_unit | mb、b或microbarn |
| uncertainty_convention | symmetric_1sigma或none |
| error_model | independent、covariance或none |

`independent`的CSV列为`energy,sigma,uncertainty`；误差必须为正的对称1σ。`covariance`和`none`的CSV列为`energy,sigma`。协方差模式另需metadata中的covariance二维矩阵，单位为原截面单位的平方、顺序与CSV一致，要求对称正定且最多1000点。none模式不计算χ²、不进行网格拟合。能量必须递增，非对称误差、能量误差和重复测点不静默处理。

`compare-experiment RUNS --data DATASET --out DIRECTORY` 输出残差CSV、图和原始χ²。统计量为 `rᵀ C⁻¹ r`，用Cholesky分解求解；不会自动给出自由度、p值或机制结论。来源：[NIST加权最小二乘](https://www.itl.nist.gov/div898/handbook/pmd/section4/pmd432.htm)。

默认要求能量测点与计算网格匹配。可以显式选`--interpolation linear`，禁止外推；插值误差未纳入协方差，也不自动拟合归一化或能量偏移。

## 扫描与拟合

`prepare-scan CONFIG --axes AXES --purpose PURPOSE --out DIRECTORY`准备网格，不运行求解器。AXES为“参数路径→明确数值数组”的JSON对象，例如键名可为`potential.a_fm`。最多3轴、30个笛卡尔积配置。

可扫描实势、现有耦合模式的连续参数和径向设置。βC变化不会自动带动βN变化；需要明确各参数的关系。通道类型或数量改变应另行生成配置，数值径向设置不能当作物理参数拟合。

`fit-grid RUNS --data DATASET --scan SCAN_RECORD --out DIRECTORY`仅分析已批准且严格通过校验的运行，报告最优采样点、离散profile、边界最优和缺失点。`--delta-chi2`定义离散支持集合，不是连续置信区间或后验。缺失点可能改变最优结论，不自动发起新计算。

## 比较、诊断与审阅

`compare RUNS --dimension couplings|potential --out DIRECTORY`要求其他因素固定。

`convergence RUNS --out FILE`比较仅径向设置不同的成功运行，默认截面相对差门槛为1%，不代表通道空间或物理模型收敛。

`diagnose RUNS --out FILE`检查存档文件和概率边界，保留原失败状态。原版`spin.dat`为平均角动量，`s-wave.dat`包含多个分波积分概率。

`review-c0 RUNS --statement ACTUAL_USER_WORDS --actor ACTOR --out DIRECTORY`在研究者明确接受后记录有限的C0概率上界例外。不得把占位词当成真实意见。该记录仅供带警告的截面对比，详见[已知限制](known-limitations.md)。

`benchmark --out DIRECTORY`用于官方固定算例的环境检查，不用于绕过研究方案的讨论步骤。
