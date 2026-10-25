# CCFULL Agent

[中文](README.md) | [English](README.en.md)

**A discussion-first research assistant for heavy-ion fusion with  CCFULL.**

CCFULL Agent organizes nuclear-structure evidence, proposes coupling schemes and potential parameters, and invokes ** Fortran77 CCFULL** after the researcher confirms a specific plan. It provides result validation, comparison, and traceable records.

**Version: v0.2.0 · Interface: conversational agent + local CLI**

The conversational agent supplies research reasoning; Python tools handle deterministic validation, execution, and analysis. No additional LLM API is required. The project does not provide a standalone chat website or background research service.

## Features

- **Nuclear-structure evidence:** fetch or import NuDat pages, extract levels, transition strengths and uncertainties, and retain provenance and hashes.
- **Coupling recommendations:** screen candidate excited states, identify evidence gaps, and derive electromagnetic deformation magnitudes using explicit transition-direction and radius conventions.
- **Planning and execution:** generate input decks, interactive answers, and batch plans; invoke the pinned original CCFULL after recording the researcher's confirmation.
- **Analysis:** compare couplings or potentials, import fusion cross sections, evaluate χ² with independent or correlated uncertainties, and produce figures and residuals.
- **Scans and diagnostics:** prepare finite parameter grids, analyze approved runs, check radial convergence, and diagnose failures.

## Installation

Requirements: Python 3.10+, NumPy, and Matplotlib. Solver execution also requires `gfortran` on PATH.

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
gfortran --version
```

Run subsequent commands from the project root. Use the source tree or an editable installation, retaining `vendor/` and the skill files. `requirements-lock.txt` provides fixed dependency versions for environment reproduction; other Python/platform combinations may need adjustments.

## Quick start

### Inspect commands and validate input

```bash
.venv/bin/python scripts/ccfull.py --help
.venv/bin/python scripts/ccfull.py validate examples/official.json
```

### Check the official example

```bash
.venv/bin/python scripts/ccfull.py benchmark --out reports/official-check
```

This runs the official configuration and two radial numerical checks, producing a verification report, raw results, and figures. It checks the solver in the current environment; it does not physically validate a new reaction system. Use a new output directory when repeating the command.

### Start a research task

Open the project in a conversational environment with skill support and have the agent read `AGENTS.md` and `.agents/skills/ccfull-research/SKILL.md`. For example:

> Use ccfull-research to retrieve nuclear-structure evidence for my specified projectile, target, and energy range. Propose coupling and potential candidates, and discuss them with me before running calculations.

You can also use the CLI directly:

```bash
.venv/bin/python scripts/ccfull.py propose --catalog path/to/catalog.json --out reports/proposal
```

`path/to/catalog.json` is an evidence file prepared for the current task, not a bundled file. Preparing a proposal does not launch the solver. Research execution requires actual acceptance of the specific plan.

See the [user guide](docs/research-tools-v02.md), currently in Chinese, for nuclear-structure, experimental-data, and scan schemas. CLI help and the [skill tool reference](.agents/skills/ccfull-research/references/tools.md) are available in English.

## Commands

| Task | Commands |
|---|---|
| Build and verification | `build`, `benchmark`, `validate` |
| Evidence and candidates | `sources`, `fetch-structure`, `import-structure`, `assess-structure`, `beta`, `beta-wu`, `propose` |
| Research planning and execution | `prepare`, `prepare-scan`, `approve`, `run` |
| Comparison and fitting | `compare`, `import-experiment`, `compare-experiment`, `fit-grid` |
| Numerical diagnostics and review | `diagnose`, `convergence`, `review-c0` |

Use `<command> --help` for arguments.

## Project structure

```text
CCFULL-Agent/
├── README.md / README.en.md
├── AGENTS.md                  # Project instructions for the agent
├── pyproject.toml             # Python project configuration
├── requirements-lock.txt      # Fixed dependency versions
├── .agents/skills/            # Research skill and physical decision guidance
├── ccfull_agent/              # Execution, nuclear data, and analysis modules
├── scripts/ccfull.py          # CLI entry point
├── vendor/ccfull/             # Pinned original source and official reference files
├── examples/official.json     # Official input configuration
└── docs/                      # Functionality and usage documentation
```

Runtime directories `build/`, `runs/`, and `reports/` are listed in `.gitignore`. Existing output directories must not be overwritten or reused.

## Scope and limitations

The wrapper supports 0⁺ ground states; one projectile vibrational mode or a rotor; up to two target vibrational modes or a rotor; and explicit target mutual-excitation choices. Nuclear deformation βN, deformation signs, vibrational or rotational assumptions, and truncation require evidence-based judgment rather than a single automatic criterion.

Transfer execution, custom anharmonic matrices, comprehensive ENSDF/EXFOR ingestion, continuous optimization, and Bayesian inference are not implemented. Grid fitting identifies the best evaluated point; a discrete Δχ² support set is not a continuous confidence interval.

Normal program termination does not establish physical validity. Failures must retain their diagnostics. A narrowly scoped researcher acceptance is recorded separately and cannot erase the original status. See [known limitations](docs/known-limitations.md) and [architecture](docs/architecture.md), currently in Chinese.

## Citation and licensing

When using CCFULL in research, cite:

> K. Hagino, N. Rowley, A. T. Kruppa, *A program for coupled-channel calculations with all order couplings for heavy-ion fusion reactions*, Computer Physics Communications **123** (1999), 143–152. [DOI: 10.1016/S0010-4655(99)00243-X](https://doi.org/10.1016/S0010-4655(99)00243-X).

Cite nuclear data and reaction-model parameters at their respective sources. See [THIRD_PARTY_NOTICES](THIRD_PARTY_NOTICES.md). A project-wide distribution license has not been assigned, and this project does not grant a new license for upstream files.
