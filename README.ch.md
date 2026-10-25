# CCFULL Agent

[中文](README.ch.md) | [English](README.md)

**A discussion-first research assistant for heavy-ion fusion with CCFULL.**

CCFULL Agent 通过对话整理核结构证据、提出耦合与势参数候选，在研究者确认具体方案后调用 ** Fortran77 CCFULL**，完成结果校验、比较和归档。

**版本：v0.2.0 · 交互方式：对话 agent + 本地 CLI**

研究推理由当前对话 agent 提供，Python 工具负责确定性的校验、执行和分析。无需额外 LLM API；项目不提供独立聊天网站或后台研究服务。

## 功能

- **核结构证据**：获取或导入 NuDat 页面，解析能级、跃迁强度及误差，记录来源和哈希。
- **耦合建议辅助**：筛选候选激发态、指出证据缺口；按明确方向和半径约定计算电磁形变幅值。
- **计算计划与执行**：生成输入、交互回答和批次计划；记录研究者确认后调用固定原版 CCFULL。
- **结果分析**：比较耦合或势参数，导入融合截面，计算独立/相关误差下的 χ²，输出图表与残差。
- **扫描与诊断**：准备有限参数网格，分析已批准的运行，检查径向收敛并诊断失败。

## 安装

需要 Python 3.10+、NumPy、Matplotlib；实际运行求解器还需要可从 PATH 调用的 `gfortran`。

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
gfortran --version
```

从项目根目录运行后续命令。当前采用源码目录或可编辑安装，需要保留 `vendor/` 和技能文件。`requirements-lock.txt` 提供固定依赖版本供复现环境时参考；不同Python/系统组合可能需要调整。

## 快速开始

### 查看命令与检查输入

```bash
.venv/bin/python scripts/ccfull.py --help
.venv/bin/python scripts/ccfull.py validate examples/official.json
```

### 校验官方算例

```bash
.venv/bin/python scripts/ccfull.py benchmark --out reports/official-check
```

该命令运行官方配置及两组径向数值检查，输出验证报告、原始结果和图表。它用于检查当前环境中的求解器，不构成新体系的物理验证。重复执行时使用新的输出目录。

### 开始研究任务

在支持项目技能的对话环境中打开本目录，让 agent 读取 `AGENTS.md` 和 `.agents/skills/ccfull-research/SKILL.md`。例如：

> 请使用 ccfull-research，根据我指定的弹靶体系和能区查找核结构资料，提出耦合及势参数候选，先讨论后运行。

也可以直接使用 CLI：

```bash
.venv/bin/python scripts/ccfull.py propose --catalog path/to/catalog.json --out reports/proposal
```

`path/to/catalog.json` 是当前任务自行准备的资料文件，不是随包文件。方案生成不会启动求解器；研究运行需要对具体计划的真实确认。

核结构、实验数据和扫描输入格式见 [使用指南](docs/research-tools-v02.md)。

## 命令

| 任务 | 命令 |
|---|---|
| 编译与校验 | `build`、`benchmark`、`validate` |
| 证据与候选 | `sources`、`fetch-structure`、`import-structure`、`assess-structure`、`beta`、`beta-wu`、`propose` |
| 研究计划与执行 | `prepare`、`prepare-scan`、`approve`、`run` |
| 比较与拟合 | `compare`、`import-experiment`、`compare-experiment`、`fit-grid` |
| 数值诊断与审阅 | `diagnose`、`convergence`、`review-c0` |

使用 `<命令> --help` 查看参数。

## 目录结构

```text
CCFULL-Agent/
├── README.md / README.en.md
├── AGENTS.md                  # agent 项目约定
├── pyproject.toml             # Python 项目配置
├── requirements-lock.txt      # 固定依赖版本
├── .agents/skills/            # 研究技能及物理判断说明
├── ccfull_agent/              # 执行、核数据与分析模块
├── scripts/ccfull.py          # CLI 入口
├── vendor/ccfull/             # 固定原版源码及官方参考文件
├── examples/official.json     # 官方输入配置
└── docs/                      # 功能与使用文档
```

运行时生成的 `build/`、`runs/` 和 `reports/` 已列入 `.gitignore`。输出目录不得覆盖复用。

## 支持范围

支持0⁺基态、一种弹核振动模式或转子、最多两种靶核振动模式或转子，以及显式互激发选择。βN、形变符号、振动/转动模型和截断需要结合证据判断，不能由单一指标自动决定。

暂不支持转移执行、自定义非谐矩阵、完整 ENSDF/EXFOR 自动入库、连续优化或贝叶斯推断。网格拟合只给出已采样点中的最优值，离散 Δχ² 集合不是连续置信区间。

程序正常结束不等于物理模型已获验证。失败结果必须保留诊断；有限范围的人工接受需另行记录，不能抹去原始状态。详见 [已知限制](docs/known-limitations.md) 和 [架构说明](docs/architecture.md)。

## 引用与许可

使用 CCFULL 开展研究时，请引用：

> K. Hagino, N. Rowley, A. T. Kruppa, *A program for coupled-channel calculations with all order couplings for heavy-ion fusion reactions*, Computer Physics Communications **123** (1999), 143–152. [DOI: 10.1016/S0010-4655(99)00243-X](https://doi.org/10.1016/S0010-4655(99)00243-X).

核数据和反应模型参数需引用各自来源。参见 [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES.md)。本项目尚未指定统一发行许可证，不为上游文件重新授予许可。
