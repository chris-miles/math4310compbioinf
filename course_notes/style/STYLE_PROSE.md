# MATH 4310 — Prose Style

Authors (human or LLM) must follow this. When in doubt, mimic an already-accepted chapter.

## Audience
Utah undergrad with linear algebra (MATH 2270), data structures (CS 2420), and one semester of cell bio (BIOL 2030). Goal: a **first honest taste of the theory** — complexity, dimensionality, noise, model assumptions. Not a graduate text (Durbin–Eddy density, Bishop-level formalism). If a passage feels graduate, pull back.

## Voice
- First-person plural: "we define," "we will show." Never "I." Use "you" only in direct exercise prompts.
- Warm but precise. Closer to a careful lecture than a paper.
- Present tense for math and algorithms; past tense only for history.
- Cut filler: "essentially," "basically," "in some sense," "note that," "it turns out that."

## Rigor
Define formally, then immediately give one sentence of intuition — in that order. State theorems precisely. Default to a **proof sketch** carrying the load-bearing idea; give a full proof only when it is short *and* the proof itself teaches a general lesson worth keeping. Mark deferred proofs explicitly. **Always** justify correctness and complexity claims, even if the justification is one sentence.

## Chapter Arc
1. Motivation (the biological problem, why naive fails).
2. Definitions → algorithm → small worked example.
3. Analysis: correctness sketch + complexity. This is where the theory taste lives.
4. Biology framing: what data, where assumptions break.
5. Pitfalls (short callout).
6. Exercises (4–8: derivation, small implementation, ≥1 failure-mode reflection).

## Theorem-likes, Callouts, and Sidenotes

For numbered, cross-referenceable environments use Quarto's crossref div syntax:

```
::: {#def-edit-distance}
The edit distance $d(x, y)$ is the minimum number of...
:::
```

Reference with `@def-edit-distance`. Quarto auto-numbers per chapter ("Definition 4.2").

| Prefix | Use |
|---|---|
| `def-` | Definition |
| `thm-` | Theorem |
| `lem-` | Lemma |
| `cor-` | Corollary |
| `prp-` | Proposition |
| `exm-` | Worked example with concrete numbers |
| `exr-` | End-of-chapter exercise |

Use Quarto callouts for student-facing emphasis:

| Callout | Use |
|---|---|
| `callout-note` | Important definitions, conventions, or model assumptions that do not need theorem numbering. |
| `callout-warning` | Pitfalls, common misconceptions, or failure modes. |
| `callout-tip` | Practical checks, coding hints, or verification habits. |
| `callout-important` | High-stakes constraints, such as data-use or reproducibility requirements. |

Reserve custom div classes for cases where Quarto callouts are a poor fit:

| Class | Use |
|---|---|
| `.algorithm` | Pseudocode / numbered steps. Title manually: "Algorithm 4.1". |
| `.biology` | Biological framing when the content should not look like a warning or note. |

Use sidenotes or margin notes for historical comments, optional connections, and interesting context that is not load-bearing for the argument. A sidenote should be short enough to skip without losing the mathematical thread.

## Math Notation
- Sequences: $x = x_1 x_2 \ldots x_n$, $|x| = n$.
- Score function $s(a,b)$; pick one symbol per matrix per chapter and stay consistent.
- $\log$ without base unless base matters.
- $\mathcal{O}(\cdot)$ for asymptotics; plain English when constants matter.
- New symbols added to `notation.qmd` appendix.

## Tables
- Markdown pipe tables for hand-written content.
- Computed tables: produce with `df.to_markdown(index=False)` inside a code chunk.
- Caption and label inline: `: My caption {#tbl-foo}`. Reference with `@tbl-foo`.
- Right-align numeric columns; left-align text. Bold the header row only (default for pipe tables).
- Keep tables narrow enough to render on mobile — no more than ~5 columns of normal text.

## Cross-References
Use Quarto labels everywhere; never "see the figure above." Standard prefixes: `@sec-`, `@fig-`, `@tbl-`, `@eq-`, plus the theorem-like prefixes above.

## Equations
Display equations get their own line. Number (with `{#eq-foo}`) only if referenced later. Always introduce a symbol the first time it appears in a chapter.

## Citations
BibTeX in `references.bib`. In-text: `@needleman1970`. Parenthetical: `[@needleman1970]`. Multiple: `[@a; @b]`. Cite all non-trivial historical or empirical claims.

## Forbidden Patterns (AI tells)
- "Not just X, but Y" — say what you mean.
- Rule of three lists where two would do.
- Inflated symbolism ("fundamental shift," "deep insight").
- Vague attribution ("widely believed," "many researchers") — cite or drop.
- Em dash overuse. One per paragraph cap; most paragraphs need none.
- Buzzwords: leverage, delve, tapestry, landscape, underscore, showcase, robust (unless statistical), comprehensive, seamless.
- Hedge stacks ("might perhaps sometimes").
- Bullet lists where a paragraph carries the same content.
- Closing recap paragraphs that restate the section.
- Student-facing notes should not explain internal course logistics unless the student needs the information to act.
- Avoid obvious contrast sentences such as "This is not X yet" or "That limitation is useful" when a direct statement of the current model is clearer.

When editing LLM drafts, invoke the `humanizer` skill.

## Stance on Software
We are not training engineers. Code serves understanding. Mention real tools (BWA, BLAST, scanpy) briefly and contrast with our pedagogical version; do not tutorial them.
