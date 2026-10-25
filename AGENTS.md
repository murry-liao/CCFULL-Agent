# CCFULL research agent

For calculations, nuclear-structure lookup, coupling recommendations or potential comparisons, read `.agents/skills/ccfull-research/SKILL.md`.

Discuss research configurations before execution. Preparing proposals and running the fixed official benchmark do not authorize a new research model. Record the researcher’s actual acceptance of the exact plan before invoking research execution; never fabricate acceptance.

Preserve pinned upstream source and reference files. Use `.venv/bin/python scripts/ccfull.py` for supported operations. Preserve logs and explicit failure status. After solver-interface changes, verify the fixed official benchmark; use appropriate checks for other changes.

No reaction-specific research catalog is bundled. Obtain evidence for the current task and pass its path explicitly to `propose --catalog`. The conversational agent supplies research reasoning; there is no external LLM API or web service in this package.

Do not overwrite completed runs. Distinguish software verification, official numerical validation and approved research results.
