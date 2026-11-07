# Research Instructions

## Objective

Pursue the central mathematical problem described in `README.md` and `PROOF.md`. Read the primary literature, fix precise definitions and hypotheses, and seek a genuinely new theorem, explicit construction, or rigorous obstruction. Partial lemmas are steps toward the main problem, not automatic stopping points.

## Method

- Work freely on mathematics: derive lemmas, make constructions, test counterexamples, and change approaches when needed.
- Verify novelty and accurately credit existing results before making a research claim.
- Separate strict proofs, conditional theorems, certified finite computations, numerical observations, and conjectures.
- Use Python, symbolic algebra, or C++ when they help answer a precise mathematical question. Do not add code merely to decorate the repository.
- Keep researching until a substantial result, a decisive obstruction, or the available research budget is exhausted. No recurring progress reports are required.

## Repository

- `PROOF.md` is the main mathematical working document. Organize it by mathematical topic, update it in place, and keep genuine proofs complete.
- `README.md` briefly records the target and independently established results. Update the result summary when new theorems are proved.
- Add `scripts/` for necessary computation and `papers/` for actual paper sources only when useful.
- Keep generated artifacts in ignored `build/`. Do not commit large search outputs, temporary files, machine logs, JSON inventories, or binaries.
- No ceremonial CI/CD, SHA-256 attestations, release packaging, automated reports, heavyweight validation pipelines, or artificial workflow infrastructure.
- Prefer short, direct technical prose. Avoid defensive boilerplate and unsupported statements of novelty.

## Pull requests

Do the research in a branch and open a PR when there is substantive mathematical content. PR descriptions should state what has been established, the exact limitations, and where the argument lives. Do not force-push shared history.
