# Research Instructions

## Mandatory English and proof-quality policy

- English is the sole language for repository-authored prose: all Markdown documentation, proof notes, README, AGENTS, manuscript exposition, source comments and docstrings, and PR descriptions. Translate existing Chinese prose throughout tracked research files, including nested directories and archives. Keep original titles/quotations and notation only where necessary, with English context; do not create parallel Chinese documentation.
- Audit the **entire existing research corpus**, not only changed files. Check every mathematical claim against its exact definitions, quantifiers, hypotheses, boundary cases, proof dependencies, citations, computational evidence, and the stated scope of conclusions. Repair gaps with complete arguments, correct mistaken claims, and add essential missing definitions, lemmas, references, or counterexamples. Do not inflate conclusions or silently turn an incomplete argument into a theorem.
- A translation must preserve equations, proof structure, assumptions, labels, bibliographic facts, and mathematical strength. If a proof is incomplete or ambiguous, resolve it mathematically where possible; otherwise identify the precise unresolved step in clear English and keep the claim at its justified status.
- Submit focused, reviewable PRs per repository. Give exact file locations, substantive corrections, and validation steps; do not mix unrelated research with translation. Preserve cross-references and all functional code, run existing relevant checks, and manually inspect for residual Chinese prose. A character scan alone is insufficient.
- Apply these standards to every future contribution. Write concise, direct English; avoid ceremonial documentation, defensive boilerplate, and unsupported novelty or correctness claims.

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
