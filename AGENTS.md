# Research Instructions

## Zero-defensive-prose rule — mandatory

Apply this rule to all repository-authored English prose: README, mathematical research notes, `PROOF.md`, LaTeX exposition, documentation, pull-request descriptions, and code comments. Enforce it when auditing historical material and when writing anything new.

**Prohibited rhetorical constructions (including grammatical variants):**
- Chinese examples to recognize and eliminate in translated source material: `不是……`, `并不是……`, `不构成……`, `还没有做……`, `尚未……`, `而不是……`, `并未……`, `不意味着……`, `不代表……`, `不能据此……`, `这并不证明……`, `我们不声称……`, `不保证……`, `仅仅是……`.
- English equivalents to eliminate from authored explanations: `not a ...`, `is not ...`, `does not ...`, `do not ...`, `has not ...`, `have not ...`, `not yet ...`, `rather than ...`, `instead of ...`, `does not constitute ...`, `does not establish ...`, `does not imply ...`, `should not be interpreted as ...`, `we do not claim ...`, `no claim is made ...`, `merely ...`, `only a preliminary ...`, `beyond the scope ...`, `not intended to ...`.
- Delete disclaimers, self-exculpation, premature qualification, repetitive status caveats, and defensive contrasts. Do not relocate them to a different paragraph or replace them with synonyms serving the same rhetorical purpose. Zero such constructions are acceptable in explanatory prose.

**Required replacement style:** Lead with the actual theorem, deliverable, tested behavior, evidence, applicability conditions, exact scope, and concrete open mathematical question, using affirmative statements. Examples:
- `We prove the classification for n >= N under hypotheses H1 and H2.`
- `Verified cases: n = 2,...,7. General-case proof target: ...`
- `CI checks: 16 passed. Proof: Proposition 4. Research priority: external review.`
- `Algorithm: FPT for disjoint ordinary cardinality quotas. Next target: robust mixed-status source-SCC accessibility.`

**Mathematical truth takes precedence over stylistic cleanup.** Preserve essential logical negation in theorems, inequalities, counterexamples, constraints, source quotations, and executable code where changing it would alter correctness. In ordinary explanation, express boundaries by positive quantified scope and precise next objectives. If a claim is unsupported, correct or remove it; never invent evidence to avoid a negative sentence.

Before opening or updating a PR, review every added or modified paragraph for prohibited rhetoric and its semantic equivalents. Record the direct result and evidence; keep mathematical limitations in compact, exact statements of domain and open problems.

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
