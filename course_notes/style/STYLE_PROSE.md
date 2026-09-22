# MATH 4310 — Prose Style

Authors (human or LLM) must follow this. When in doubt, mimic an already-accepted lesson; Lesson 1 is the reference implementation of these rules.

Rules come in three tiers. **Strict** rules are mechanical and are never bent: markup, ordering, one lesson per lecture. **Default** rules are the normal case; break one only with a reason you could say out loud, and if a lesson breaks two of them, split it. **Free** choices are yours.

## Audience

Upper-division undergraduates with linear algebra (MATH 2270), introductory programming (CS 1410), and basic cell biology (BIOL 2020 or 2021). Basic probability ideas are probably familiar, but a formal probability, statistics, algorithms, genetics, or optimization course is not assumed. Treat students as capable learners with different backgrounds.

Use a pedagogical textbook or lecture-notes tone: kind and generous to the reader, precise, and economical. Do not reteach elementary probability at length or announce that a concept is easy. Briefly recall notation when it matters, then spend the explanation on the new applied step: what is conditioned on, why a factorization holds, why a recurrence needs a state. Advanced ideas are welcome when the lesson supplies a concrete reason, a definition, a small calculation, and the load-bearing reasoning.

The goal is a **first honest taste of the theory** — complexity, dimensionality, noise, model assumptions. Not a graduate text (Durbin–Eddy density, Bishop-level formalism). If a passage feels graduate, pull back.

What we want students to leave with is not the list of techniques, which will change, but the habit of asking, for any new method: what is the data as a mathematical object, what question is being asked of it, what model is being assumed, and where that model breaks. We do not say this out loud on every page. The structure of every lesson says it.

## These Are Lecture Notes

Each numbered lesson corresponds to one class day: one ~75-minute meeting, including board work and discussion. Use actual meeting dates in `site/schedule.json` and `planning/lesson_calendar.md`; a holiday week can have only one meeting. Count exams and project sessions before assigning lessons.

The model is a good set of scribed lecture notes (MIT OCW 6.047 is the reference point): numbered sections whose headings name the topic, and under each heading the content itself, in blocks: definitions, a problem statement, an algorithm, a worked example, a result with its proof sketch, a figure, exercises. Prose supplies enough reasoning for students to follow the notes after class. Use short paragraphs to connect the data, objective, assumptions, and calculation. Avoid transcribing a lecture or leaving bare formulas that require an unstated lecture to explain them. A passage with three consecutive prose paragraphs deserves review, not automatic deletion.

Textbook structure, lecture-notes word count. A reader scanning the sidebar should know exactly what was covered; a reader with fifteen minutes should be able to read the whole page.

## Lesson Shape

Every lesson presents its topic in this order (strict):

1. **Core question.** One sentence under the title. Already the convention.
2. **The data as a mathematical object.** The definition of what we compute on. Where the data comes from, if it needs saying, is two sentences in a `.biology` div here, not a motivation essay.
3. **The problem.** Input, output, what is optimized. One display line or one short paragraph.
4. **The model.** The assumptions that make the problem answerable, usually as a definition (a scoring function, a Markov chain) rather than prose.
5. **The algorithm or derivation.** Pseudocode in an `.algorithm` div and the Python, or the derivation.
6. **A worked example.** One or two `exm-` blocks with concrete numbers, doable by hand in five minutes.
7. **Correctness and cost.** A numbered proposition or theorem with a proof sketch; the complexity in one line.
8. **Limitations.** A plain section of short numbered paragraphs, each naming the assumption that fails and the lesson that addresses it. Not warning boxes.
9. **Exercises.** Two to four `exr-` blocks.

The opening must briefly establish what the observed data are, why the biological question matters, and what output the lesson will produce. Lesson 1 is the benchmark. Usually one or two short paragraphs suffice; a continuing lesson can say what changes from the previous one. A named case or dataset is useful when it clarifies the question, but do not invent a historical hook or force one into every lesson. Cite empirical or historical claims. Hypothetical examples should be recognizably hypothetical. The opening should lead directly into the mathematical object.

### Titles and headings

Lesson titles name a recognizable purpose or biological question before specialized machinery: "Finding the Most Likely Genome Annotation" and "How Certain Are Our Genome Labels?" are more inviting than technique names alone. Name Viterbi, forward-backward, and other standard methods in the relevant section so students learn and can find the terminology.

Section headings name the topic, as a textbook would: "Hamming distance", "The alignment recurrence", "Eulerian paths", "Limitations". Never narrative or rhetorical headings ("What the rule assumed", "A better rule, and the first proof", "Three things you can do with a sequence"). Never headings that name the pedagogical step ("The model", "The question") unless that is literally the topic.

### Default budget per lesson

| Element | Default |
|---|---|
| Content sections (between Core question and Limitations) | 4 to 6, roughly 10 to 15 minutes of class each |
| Prose outside numbered environments, code, and exercises | 500 to 900 words |
| Whole page, everything included, excluding code | 1,400 to 1,900 words |
| Numbered definitions | 3 to 5 |
| Numbered theorem-likes | 1 or 2, each with a proof sketch |
| Worked examples | 2 or 3 |
| Code chunks | 2 or 3, short |
| Figures | 0 or 1 |
| Callouts of any kind | at most 2 |
| Exercises | 3 to 5 |
| "Going deeper" blocks | 0 to 2 |

These are pacing guides, not quotas. Never pad a lesson by repeating definitions, adding generic framing, or expanding every proof. Conversely, a correct page of formulas is not a full class. Check whether students have time to interpret the input, work a small example, follow the derivation, and test a limitation.

Length is calibrated to a 75-minute class with board work and discussion: Lesson 1 at five content sections and about 1,700 words is the reference. A page that could be covered in 35 minutes is missing a section; a page over 2,000 words needs a pacing review. Density matters more than any one number: if a paragraph can become a definition, an example, a table, or a display equation, it should.

### Lesson sources and front matter

```yaml
---
published: true
title: "Recovering and Checking an Alignment"
subtitle: "Lesson 4 · 2027-01-25 · Week 3"
nocite: |
  @durbin1998, @compeau2015
---
```

The subtitle carries the schedule. Active sources are `course_notes/lessons/lessonNN-short-slug.md`; `.qmd` files are legacy copies and are not built. See `site/AUTHORING.md`. A `published: true` flag makes a page available in a local build; committing/pushing remains a separate publication decision.

## Voice

- First-person plural: "we define," "we will show." Never "I." Use "you" only in direct exercise prompts.
- Warm but precise. Closer to a careful lecture than a paper.
- Present tense for math and algorithms; past tense only for history.
- Cut filler: "essentially," "basically," "in some sense," "note that," "it turns out that."
- Do not name the course's organizing idea as a theme. The lesson order embodies it; the prose does not announce it.

## Rigor

Define formally, then give one sentence of intuition — in that order, and the intuition goes in the prose after the definition block, not inside it. That keeps the numbered block quotable and the intuition skimmable.

State theorems precisely. Default to a **proof sketch** carrying the load-bearing idea (optimal substructure, conditional independence, orthogonality). Give a full proof only when it is short *and* the proof itself teaches something general; put it in a "Going deeper" block. Mark deferred proofs explicitly. **Always** justify correctness and complexity claims, even if the justification is one sentence.

Examples go numbers first, general statement second. If a worked example takes more than five minutes by hand, it is a problem-set problem.

Use just-in-time reminders: expectations and likelihood ratios in Lesson 2, big-O in Lesson 3, conditional and sequence probabilities in Lesson 9. Assume familiarity with basic ideas and explain the specific modeling use. Link to the combined computing and math reference for a reminder. Do not assume prior knowledge of constrained optimization, EM, or information theory; introduce the needed idea or put the advanced derivation in an optional block.

## Environments (strict)

Numbered, cross-referenceable environments use Quarto's crossref div syntax. The first line inside the div is a `##` heading that becomes the environment's name.

```
::: {#def-alphabet}
## Alphabet and sequence
An alphabet $\Sigma$ is a finite set of symbols. A sequence over $\Sigma$ ...
:::
```

Reference with `@def-alphabet`. Quarto numbers per lesson ("Definition 3.1"). Slugs are global, so name them after the concept, not the lesson: `def-alphabet`, not `def-lesson01-alphabet`. Lesson numbers change; concepts do not.

| Prefix | Use |
|---|---|
| `def-` | Definition |
| `thm-` | Theorem |
| `lem-` | Lemma |
| `cor-` | Corollary |
| `prp-` | Proposition |
| `exm-` | The lesson's worked example |
| `exr-` | End-of-lesson exercise, one task each; multi-part problems belong in problem sets |

Proofs and proof sketches go in `::: {.proof}` directly after the statement. Quarto labels it "Proof." Open a sketch with "Sketch." as the first word.

Unnumbered structure:

| Markup | Use |
|---|---|
| `::: {.algorithm}` | Pseudocode. First line bold: `**Algorithm: fill the alignment table.**` No manual numbering; refer to algorithms by name. |
| `::: {.biology}` | The "where the data comes from" note and other biological framing. At most one or two per lesson. |
| `::: {.callout-warning title="..."}` | One misconception worth interrupting the page for. Limitations go in the plain "Limitations" section, not here. |
| `::: {.callout-note title="..."}` | A convention or assumption that does not deserve definition numbering. |
| `::: {.callout-tip title="Check: ..."}` | A verification habit: the small input that would catch a wrong implementation. |
| `::: {.callout-note collapse="true" title="Going deeper: ..."}` | Full proofs, optional derivations, connections to the literature. |
| `::: {.callout-important}` | Reserved for data-use and reproducibility requirements. |

Use margin notes for history and optional connections that are not load-bearing. A margin note must be skippable without losing the thread.

Do not use `callout-note` with a "Definition:" title in place of a `def-` block. That was the pattern in early drafts and it does not number or cross-reference.

## Web Style and Branding

- The website uses University of Utah core colors for major UI: Utah Red `#BE0000`, black, white, and neutral greys.
- Accent colors are subordinate. Use them only for small elements such as callout borders, charts, or infographics. Do not use accent colors for headings, links, page backgrounds, or other major brand-bearing elements.
- Avoid blue or purple as stand-alone interface colors; they can read as other regional university colors.
- Headings use Montserrat Bold when available. Body copy uses Source Sans when available, with system sans-serif fallbacks.
- Use the website format deliberately: callouts, margin notes, collapsible code, and accessible figure captions can carry structure that a PDF cannot.

## Math Notation

- Sequences: $x = x_1 x_2 \ldots x_n$, $|x| = n$.
- Score function $s(a,b)$; pick one symbol per matrix per lesson and stay consistent.
- $\log$ without base unless base matters.
- $\mathcal{O}(\cdot)$ for asymptotics; plain English when constants matter.
- New symbols added to `notation.qmd` appendix.

## Tables

- Markdown pipe tables for hand-written content.
- Computed tables: produce with `df.to_markdown(index=False)` inside a code chunk.
- Caption and label inline: `: My caption {#tbl-foo}`. Reference with `@tbl-foo`.
- Right-align numeric columns; left-align text.
- Keep tables narrow enough to render on mobile — no more than ~5 columns of normal text.

## Cross-References

Use Quarto labels everywhere; never "see the figure above." Standard prefixes: `@sec-`, `@fig-`, `@tbl-`, `@eq-`, plus the theorem-like prefixes above.

## Equations

Display equations get their own line. Number (with `{#eq-foo}`) only if referenced later. Always introduce a symbol the first time it appears in a lesson.

## Citations

Keep bibliography entries in `references.bib`. Cite non-trivial historical or empirical claims where they occur. Credit books and resources actually consulted using front-matter `nocite` and a compact `#refs` div at the bottom:

```markdown
::: {#refs}
**References**
:::
```

The website renders the bibliography in small, readable type. Do not add a "Further reading" section or prose claiming exact chapter-by-chapter provenance. References give credit; they are not assigned reading. Only cite sources actually checked. Keep exact internal source notes in ignored working material when useful.

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
- Template headings that name the pedagogical step rather than the content ("## The Model", "## The Question"), and narrative headings that tease rather than name ("## What the rule assumed").
- Stacked callout boxes. Two boxes in a row is a list wearing a costume.
- Motivation essays. Prefer one or two short paragraphs that establish data, purpose, and output.
- Dramatic beats ("Except that it is not.", "Nothing went wrong in the arithmetic.").
- Student-facing notes should not explain internal course logistics unless the student needs the information to act.
- Obvious contrast sentences such as "This is not X yet" or "That limitation is useful" when a direct statement of the current model is clearer.

When editing LLM drafts, invoke the `humanizer` skill.

## Stance on Software

We are not training engineers. Code serves understanding. Mention real tools (BWA, BLAST, scanpy) briefly and contrast with our pedagogical version; do not tutorial them.

## Assignments and validation

Allow a student's preferred language unless a particular task needs a restriction. Python and scientific packages describe the lecture examples, not a universal submission requirement. Chris will rarely inspect code line by line; emphasize the reasoning, independently justified checks, actual outputs, and interpretation. Code should remain reproducible and explainable.

For an introductory problem set, aim for 3 or 4 problems and about one hour: small hand calculations or proofs, a compact coding or debugging task, and a concrete interpretation question. Let students use AI where permitted by the course policy. Ask for meaningful validation and a brief disclosure, not long process narratives. No take-home question is literally LLM-proof; unsupported claims of AI resistance do not belong in student material. Do not introduce new graded quizzes or assessment weights without a course-level decision.
