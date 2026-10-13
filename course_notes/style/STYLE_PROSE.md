# MATH 4310: Prose style

These rules govern drafting and revision by human and LLM authors. A student should be able to read a section after class and explain why its method is needed, how it works, and what its answer means. The prose must supply the reasoning needed to do that.

**Required** means a condition of acceptance. **Default** means a starting point that may change to suit the topic. Mathematical accuracy, explicit reasoning, and the markup conventions are required. Section order and length are defaults. When space is tight, reduce the scope before cutting an explanation the reader needs.

## Audience and tone

Write for upper-division undergraduates who know basic computational linear algebra and basic programming. Assume they can work with matrices and matrix multiplication; basic eigenvalue ideas may be familiar, but explain any new use. Do not reteach for loops or other elementary programming constructs. Their biology background is basic cell biology, with little or no bioinformatics. Explain the relevant biological context and data-producing technology rather than assuming familiarity with sequencing, genome annotation, expression assays, or genetic inference. A formal probability, statistics, algorithms, genetics, or optimization course is not assumed; introduce the probability and modeling ideas needed for the calculation where they arise.

The target is concise, self-contained lecture notes. Use enough connected prose to set up the problem and explain the reasoning, without expanding into book-length background or treating every topic as a survey. Concision does not permit unexplained formulas, missing assumptions, or skipped mathematical arguments. Preserve mathematical and applied depth through focused derivations, worked examples, and interpretation.

Use the voice of a patient instructor working through a precise question. Treat students as capable learners with different backgrounds. State what we know, explain the next step, and give them something concrete to check. Warmth comes from anticipating a reasonable difficulty and resolving it. Avoid calling a step "obvious," "trivial," or "easy."

Advanced ideas are welcome when the lesson supplies the needed background and reasoning. Explain a new use of conditioning or a new optimization constraint where it arises. Do not send students to an appendix to learn a prerequisite for the next line, or interrupt an application with a general review of elementary probability.

The target is a first careful treatment of the theory: complexity, noise, dimensionality, and model assumptions. Keep the main argument at undergraduate level. Do not assume constrained optimization, EM, or information theory; introduce the needed idea or make the advanced extension optional. The computing reference is for setup and notation conventions, not compressed mathematical prerequisite courses.

Students should learn to identify the data, the question, the model, and the limits of an answer. Develop those habits through the calculations and examples. Do not repeatedly announce them as the course's philosophy.

This may be a student's only course in computational bioinformatics methods. As a general design consideration, favor examples that give practice with ideas and technical skills transferable to other biological problems as tools change. When it helps the lesson, connect that practice to judging AI-assisted analyses: understanding inputs, checking assumptions and outputs, and interpreting evidence. Neither transfer nor AI needs an explicit discussion in every lesson; avoid forced connections and claims that a method is future-proof.

## The reference and what to learn from it

The writing reference is Peter Holderrieth and Ezra Erives, [*An Introduction to Flow Matching and Diffusion Models*](https://arxiv.org/abs/2506.02070) (MIT 6.S184 notes); study §1.3, §2.1, and §§3.1 to 3.3. For biological motivation and algorithmic intuition, Compeau and Pevzner's *Bioinformatics Algorithms* (chapter 5 in particular) is the model. What makes both clear:

- **Motivate before stating.** Compeau and Pevzner show that inserting gaps reveals conserved columns before naming alignment; Holderrieth and Erives raise the reader's doubt ("At first sight, a conditional vector field seems useless") before the theorem that resolves it.
- **Say where the argument stands.** Each section says what has been built and what is still missing, including the exact computational obstacle ("we cannot compute it efficiently").
- **Explain why, then digest.** A paragraph of intuition accompanies the proof, and formulas are read as operations ("First, draw a random time.... Second, ...").
- **Draw why.** Their figures show objects behaving (paths, flows, a game played step by step), not only the data structures.
- **Show what matters.** A few "Key Idea" boxes and section summaries separate the central moves from bookkeeping.

Their promotional openings, emphatic adjectives, and long technical proofs are not models. Our Python examples must also meet `STYLE_CODE.md`.

## Consistency across lessons

**Required: the lessons must read as parts of one book.** Before drafting or revising, read the relevant planning document, the neighboring accepted lessons, and the prose, code, figure, and schematic guides. `course_notes/lessons/lesson01-first-problem.md` is the reference for established course conventions. Use its concrete examples and direct explanations as a starting point, while applying the clarity requirements here to every lesson.

- Keep a common voice and level of detail. A change of topic should not produce a change from patient explanation to terse research prose or a chatty tutorial.
- Carry terminology and notation forward. Use the same name and symbol for the same object across lessons; explain any necessary change. Add shared notation to `course_notes/appendices/notation.md`.
- Reuse established examples across connected lessons. Preserve their strings, scores, state names, and index conventions unless the lesson explicitly changes them. The example, code, and figure must describe the same instance.
- Preserve the shared presentation: Core question, recognizable topic headings, numbered environments, proof placement, cross-references, executable code conventions, and end-of-lesson exercises. Use the common figure style and schematic helpers, including their established color meanings.
- Read the transition from the preceding lesson and into the following one. Recall what students already have, and identify what this lesson adds. Avoid silently redefining a model or presenting an established result as new. Keep references to other lessons mostly in these bridges; a forward reference does not replace a local explanation.

Section order may vary to explain the subject. Local formatting conventions should not. If an accepted lesson conflicts with a current style rule, follow the rule and address the inconsistency deliberately; do not copy the discrepancy into another lesson.

## Motivation and bridges within a lesson

**Required: every substantial section must make its purpose clear before asking the reader to absorb new machinery.** The purpose may come from biology, an unresolved computation, or a new question about an existing result. A short continuation can inherit its purpose from the preceding paragraph. Do not manufacture a fresh biological story for every subsection.

**Required: set the stage before introducing the mathematical method.** Explain which experiment or technology supplies the available data and what a record literally looks like: for example, measured bases become letter strings, reads are fragments without genomic coordinates, and an expression matrix records measurements for named genes across samples or cells. State the biological question and the inference we want, then translate it into a concrete mathematical problem with explicit inputs, desired output, and assumptions. Solve that problem and interpret what the result supports biologically. Explain only the technological and biological details needed for this chain of reasoning; a continuing lesson may recall an established setup instead of rebuilding it.

Use this progression in the lesson's main explanation, not solely in a biology callout. "Sequence analysis is important in modern biology" gives the reader no problem to solve. Mathematical definitions and algorithms should answer a question the reader already understands.

Use Lesson 1's sequencing-to-read comparison as the clarity benchmark throughout the book. At a topic's first appearance, briefly explain what the technology measures, show what one resulting record looks like, and state what we want to learn biologically. Terms such as read, reference, alignment, expression count, and reporter need that orientation before they can carry an explanation. A continuing lesson can recall the established data and pose its new question. Keep the Core question understandable before the lesson's technical definitions, and use a small sketch or comparison when it can show the problem more clearly than another paragraph.

**Required: keep the biological meaning connected to the work throughout when informative.** Explain how modeling choices represent the data, why a computational obstacle matters for the inference, what an algorithm's behavior means for the biological question, and what results and limitations allow us to conclude. An opening motivation followed by mathematics with no return to that question is insufficient. Do not force a biological aside into every algebraic step; place these connections where they help the reader understand a choice or an answer.

Use a named dataset or historical example when it clarifies that inference. Cite empirical and historical claims. Mark toy situations with language such as "suppose" or "in this example." Biological motivation belongs in the main explanation wherever it is needed; it is not confined to an opening `.biology` box.

A bridge should identify the specific relationship between two sections. Common relationships include:

- A biological question calls for an output we have not yet defined.
- We have defined an objective but cannot yet compute its optimum efficiently.
- We can compute a score but still need to recover the object that attains it.
- A modeling assumption excludes a case we now want to handle.
- We have one estimate but need to measure uncertainty about it.

State the relationship in the vocabulary of the topic. For example:

> The alignment table gives the optimal score. To identify which bases are paired, we also need an alignment that attains that score. We recover one by recording which predecessor supplies each table entry.

This explains why traceback is needed and what it will do. "Next, we explore traceback" only announces the itinerary.

Bridges can occur at the opening of a section, between two displays, after a worked example, or across lessons. Explain a limitation when it motivates the next step; do not save it for a final section. Equally, do not invent a failure in a sound method to introduce another topic. A different biological question is sufficient reason to change methods.

## Developing a section

Use the following sequence as a default, not as a set of headings:

1. Establish the local question and recall only what is needed to answer it.
2. Describe the proposed idea in ordinary language, possibly on a small instance.
3. Give the precise definition, equation, or algorithm.
4. Explain the reasoning that makes it work and show its use on an example.
5. Interpret the result, including the assumption that matters for its biological use.

A definition may follow a motivating example. A proof may appear during the derivation. Put these pieces where they answer the reader's next question. Do not defer all examples until after all the theory or repeat the entire sequence for each minor definition.

Give each paragraph a specific job: pose a problem, explain a choice, derive a step, or interpret a result. Develop that job in connected sentences. Paragraph length follows the reasoning; three consecutive prose paragraphs are acceptable when each is necessary. Blocks, tables, and equations organize information, but they do not replace explanations of why one idea follows from another.

At the end of a section, state the consequence the reader will use next, if it is not already clear. A useful short recap can gather a multi-step construction or distinguish similar quantities. Avoid recaps that merely repeat the heading and announce that the method is useful.

### Definitions and notation

Explain why a new object is needed before or alongside naming it. Keep a numbered definition precise enough to quote, with assumptions and domains included. Place motivation and interpretation in adjacent prose; neither has to fit into one sentence. The student should know both what the object is and what work it will do.

Introduce symbols as they become necessary. State what an index counts, what an array dimension represents, and whether a quantity is observed, chosen, estimated, or unknown. For probability expressions, identify what is random and what is conditioned on. Retain the same terms and symbols across the explanation, example, figure, and code.

Recall a definition or formula when the next argument depends on it. A cross-reference alone may require too much reconstruction. Repeat the needed part, then explain its new use. Do not repeat the full earlier derivation.

### Equations and derivations

**Required: state important mathematics explicitly when symbols or equations clarify the claim.** Objectives, assumptions, model definitions, recurrences, comparisons, and results should not remain only verbal ideas when a concrete mathematical statement makes them precise. Define symbols locally and explain the statement's meaning; keep it consistent across prose, worked examples, code, figures, and adjacent lessons. Use notation for clarity, not decoration or a quota.

**Required: a central display must have a stated purpose and an explanation of what it says.** Introduce what we are calculating or establishing. After the display, explain the new operation, relationship, or consequence. Reading every symbol aloud is unnecessary; interpreting a new sum, conditioning operation, or optimization is often essential.

For an algorithm, define what its stored state means before giving an update. Explain each case, the boundary values, and why the dependencies are available when needed. A recurrence without a state interpretation leaves students guessing what the table contains.

For a derivation, state the target before starting algebra. Give the reason for each step that uses a modeling assumption, theorem, change of variables, or optimization argument. Name the assumption at the point of use: "Given the hidden state, the model makes this emission independent of earlier observations." "By standard arguments" does not supply that reasoning.

Short annotations beside a display can identify cases or operations. Split a long derivation where the mathematical task changes. At its end, say what the result lets us calculate or conclude. Avoid an uninterrupted chain of equalities whose purpose appears only afterward.

### Examples and interpretation

Use one small running instance through the definition, calculation, code, and figure when possible. Reusing an instance lets the reader concentrate on the new idea. Change the instance when the change itself teaches something, and explain why.

A worked example must expose a decision or calculation, not just report a result. Show one table entry, one update, or one competing pair of candidates in detail. Keep instances small enough to trace by hand, but choose them so the new idea changes the answer; an example on which the old method already succeeds cannot show why the new one is needed. Explain what the answer means in the example, including what it supports biologically.

Most hand examples should take about five minutes. A longer example is appropriate only when its successive steps teach distinct parts of the method; otherwise move it to a problem set. Include a boundary case or small counterexample when it reveals an assumption that a successful example would hide.

### Correctness, assumptions, and cost

State claims precisely and justify them. Use numbered propositions or theorems for the central results. Default to a proof sketch containing the step on which correctness depends. For dynamic programming, saying "optimal substructure" is insufficient: explain why replacing a prefix by a better prefix would improve the whole solution. A short full proof belongs in the main text when that argument is the lesson. Put longer optional proofs in "Going deeper" blocks and identify any result used without proof.

Keep the model, the algorithm, and the biological conclusion distinct. An exact optimum under a scoring rule need not be the true evolutionary history. State such a qualification once, where it changes an interpretation, after saying what the result does support. A disclaimer after every result buries the claims it qualifies.

Justify complexity by counting the operations that dominate the computation. For a table algorithm, give the number of entries and the work per entry; include storage when it matters. State a fixed-parameter or unit-cost assumption if the bound depends on it. Do not use "efficient" as a substitute for this explanation.

## Sentences and voice

- Prefer concrete subjects and verbs: "The recurrence compares three possible final columns" rather than "A comparison of alternatives is performed."
- Begin with a familiar object when possible, then explain the new fact about it. Avoid making the reader resolve several new terms before reaching the verb.
- Use "we" for a shared mathematical action: "We condition on the final state." Use direct statements for facts. Avoid repeated "we will now" announcements. Do not use "I." Reserve "you" for direct exercise prompts.
- Keep technical names stable. A state, a label, and a class are not interchangeable unless the text explicitly identifies them. Repeat the noun when "this" or "it" could refer to several objects.
- Use "because," "therefore," and "under this assumption" when they state an actual logical relationship. A transition word cannot repair a missing reason.
- Prefer familiar wording: "use," "compute," "choose," "assume," and "estimate." Keep technical vocabulary when it adds precision; define it at first use.
- Be direct about uncertainty. Say what the model assumes or what evidence supports an interpretation. Avoid both overclaiming and stacks of hedges.
- Use present tense for mathematics and algorithms. Reserve boldface for a small number of terms or distinctions worth finding again.

Concrete questions can help when they express a real problem: "Which previous state could produce this observation?" Follow with the reasoning. Avoid theatrical questions, manufactured suspense, and reassurance that a topic is simple. The default voice is calm, interested, and explanatory.

Cut empty framing such as "It is important to note that" or "This powerful framework allows us to." Replace it with the actual fact. Do not delete useful reminders, comparisons, or explanations merely because they use "recall," "in other words," or "for example."

## Presenting algorithms and code

Code completes an explanation. Before the implementation, the reader should know its inputs, output, and mathematical rule. Use pseudocode in an `.algorithm` div before Python, as required by `STYLE_CODE.md`. Keep both in the order of the derivation and use names that make the correspondence visible.

Explain any difference between mathematical and programming conventions. For example, if `table[i][j]` stores the distance between prefixes of lengths $i$ and $j$, then the final letters of those prefixes are `x[i - 1]` and `y[j - 1]` in Python. That one sentence resolves a common obstacle more effectively than comments restating each assignment.

Code should read as a summary of the algorithm. Show a call on the running example with its actual output and interpret it; when the hand calculation already gives the answer, note the agreement in a sentence. Use code to explore as well as check: simulate from the model, vary a parameter, or plot a behavior.

Keep the essential algorithm visible. Put imports and drawing setup in setup chunks; use the shared figure helpers. Comments explain a modeling choice, an index convention, or the link to an equation. Do not bury the recurrence in a library call or turn the lesson into a tour of Python syntax. Introduce unfamiliar syntax only when students need it to understand or adapt the code.

## Examples of the expected prose

These examples illustrate explanatory choices, not text to repeat in every lesson.

### A biological reason for a definition

Weak: "We now introduce edit distance, a powerful measure of sequence similarity."

Better:

> Suppose a read contains one extra base relative to a reference segment. Comparing the strings position by position pairs the later bases with the wrong neighbors. We need a comparison that allows an insertion to shift the correspondence. Edit distance gives each insertion, deletion, or substitution a cost of one and finds the smallest total cost.

The revised opening identifies a concrete difficulty and explains the modeling choice. A precise definition can now follow.

### A recurrence with its reasoning

Weak: "By optimal substructure, the edit-distance recurrence follows."

Better:

> Let $D(i,j)$ be the edit distance between the first $i$ letters of $x$ and the first $j$ letters of $y$. An alignment of these prefixes ends in one of three ways: two letters paired, a letter of $x$ paired with a gap, or a gap paired with a letter of $y$. Removing that final column leaves a smaller alignment. Its cost must be minimal for its prefixes; otherwise replacing it with a cheaper alignment would reduce the total cost. For $i,j\geq 1$, this gives
>
> $$
> D(i,j)=\min\begin{cases}
> D(i-1,j-1)+\mathbf{1}[x_i\ne y_j], & \text{pair the final letters},\\
> D(i-1,j)+1, & \text{delete the final letter of }x,\\
> D(i,j-1)+1, & \text{insert the final letter of }y.
> \end{cases}
> $$
>
> The indicator $\mathbf{1}[x_i\ne y_j]$ is 1 for a mismatch and 0 for a match. Each candidate extends a valid smaller alignment, and the three cases cover every possible final column. Taking their minimum therefore gives the optimum. The boundaries are $D(i,0)=i$ and $D(0,j)=j$: converting a prefix to or from the empty string takes one edit per letter.

The explanation supplies the state meaning, the reason the cases suffice, and the boundary interpretation. A table and a hand calculation can now make the recurrence concrete.

### An output with an independent check

Weak: "The code returns 2, demonstrating the effectiveness of the algorithm."

Better:

> Deleting G and one T turns `ACGTT` into `ACT`, so two edits suffice. The lengths differ by two, so fewer than two edits cannot suffice. The distance is therefore 2. The function also returns 2 on this input, agreeing with the hand calculation.

This explains both the expected answer and what the code check establishes.

### A bridge to uncertainty

Weak: "Having explored Viterbi decoding, we now delve into the forward-backward algorithm."

Better:

> Viterbi decoding finds a most probable complete sequence of hidden states. For a particular base, we may instead want the probability that it belongs to each state, accounting for all possible state sequences. Forward-backward computes these probabilities under the same HMM.

The new question motivates a different computation without suggesting that Viterbi failed at its own task.

## Lesson shape and pacing

Write self-contained lecture notes with recognizable sections, concise prose, and clearly marked definitions, examples, and results. A reader scanning the headings should see what the lesson covers. A reader studying the argument should find the explanations needed to reconstruct it after class.

Each numbered lesson corresponds to one 80-minute meeting, including board work and discussion. Use the actual dates in `site/schedule.json` and `course_notes/planning/lesson_calendar.md`; count holidays, exams, and the flex meeting before assigning lessons.

Begin with a one-sentence **Core question** under the title. The opening should establish the biological or computational purpose, the observed data, and the output the lesson will develop. Usually one or two paragraphs suffice. A continuing lesson can recall the existing setup and explain the new question.

Across its content sections, a lesson must identify the mathematical objects, state the problem and assumptions, develop the method, work through an example, and justify its claims and relevant costs. Interleave these components as the explanation requires. Existing topic outlines can remain, but their headings do not excuse missing motivation or reasoning within a section.

Discuss limitations where they affect a claim or motivate another method. A final **Limitations** section is useful when several unresolved assumptions need to be gathered; it is optional if those assumptions are already clear in the discussion. Use plain paragraphs or a short numbered list. Name the assumption, explain its consequence, and point to a later treatment only when one exists.

End with exercises that ask students to use the central ideas. Three to five short tasks are a reasonable default. Each numbered exercise should have one coherent task; longer multi-part investigations belong in problem sets.

### Titles and headings

Retain recognizable topic names: "Global Alignment," "Markov Chains," "Viterbi Decoding," and "PCA." Add a purpose when it helps: "Uncertainty in HMMs: Forward-Backward." Section headings should identify their content, such as "The alignment recurrence" or "Recovering an alignment." Avoid generic headings such as "The model" and teaser headings such as "A better way." The prose beneath the heading supplies the motivation.

### Default budget

Prepare roughly one hour of taught material within the scheduled meeting. This is a teaching-duration target, not an hour of reading or a word-count target. Lesson 1 is a useful depth benchmark and may be slightly short. Some existing lessons need more development to support that hour. Concise, self-contained notes must supply enough derivation, worked examples, biological interpretation, and practice for a full lesson; they must not become skeletal summaries. Estimate time for explanation, board work, tracing an algorithm, and student calculation as well as prose. Add material that advances the lesson's question, without padding. Spend most of the lesson on its one or two central ideas; a side topic that would get only a paragraph is usually better cut or made an exercise.

Use these ranges to keep the scale of lessons consistent. They are pacing defaults, not quotas or acceptance tests.

| Element | Default per lesson |
|---|---|
| Numbered definitions | 2 to 4; define minor terms in a sentence |
| Numbered propositions or theorems | 1 or 2, each with justification |
| Worked examples | 2 or 3 |
| Code chunks | 2 or 3, short |
| Figures | 1 to 3; at least one shows the central idea at work, not only computed values |
| Interruptive callouts | At most 2, including Key idea boxes; numbered environments are not callouts |
| Exercises | 3 to 5 |
| "Going deeper" blocks | 0 to 2, genuinely optional |

Never add a definition, example, or figure to meet a count. Use more when the explanation requires them, and check that the lesson still fits its time budget. There is no separate word allowance for prose outside mathematical environments: essential reasoning belongs wherever the reader needs it.

If a lesson is too long, remove a secondary topic or move a technical extension to optional material. Preserve the sentences that explain why a definition appears, how a derivation works, or what a result means. If a lesson is short, check whether it leaves enough to interpret and practice before adding another topic. Lesson 1 provides examples of local conventions and reused instances, not an exemption from this guide's clarity checks.

### Sources and front matter

```yaml
---
published: true
title: "Local Alignment"
subtitle: "Lesson 4 · 2027-01-25 · Week 3"
nocite: |
  @durbin1998, @compeau2015
---
```

The subtitle carries the schedule. Active sources are `course_notes/lessons/lessonNN-short-slug.md`; `.qmd` files are legacy copies and are not built. See `site/AUTHORING.md`. Set `published: true` only after reviewing the page. Draft sources are still public if committed: keep private references and instructor-only working material in ignored locations. Committing and pushing are separate publication decisions.

## Environments (required)

Numbered, cross-referenceable environments retain Quarto's crossref div syntax in the current Pandoc/Python build. The first line inside a named environment is a `##` heading that becomes its name.

```
::: {#def-alphabet}
## Alphabet and sequence
An alphabet $\Sigma$ is a finite set of symbols. A sequence over $\Sigma$ ...
:::
```

Reference with `@def-alphabet`. The site builder numbers per lesson ("Definition 3.1"). Slugs are global, so name them after the concept, not the lesson: `def-alphabet`, not `def-lesson01-alphabet`. Lesson numbers change; concepts do not.

| Prefix | Use |
|---|---|
| `def-` | Definition |
| `thm-` | Theorem |
| `lem-` | Lemma |
| `cor-` | Corollary |
| `prp-` | Proposition |
| `exm-` | The lesson's worked example |
| `exr-` | End-of-lesson exercise, one task each; multi-part problems belong in problem sets |

Proofs and proof sketches go in `::: {.proof}` directly after the statement. The site builder labels it "Proof." Open a sketch with "Sketch." as the first word.

Unnumbered structure:

| Markup | Use |
|---|---|
| `::: {.algorithm}` | Pseudocode. First line bold: `**Algorithm: fill the alignment table.**` No manual numbering; refer to algorithms by name. |
| `::: {.callout-note title="Key idea"}` | One or two sentences stating a central idea of the lesson, placed where it has been earned. |
| `::: {.biology}` | A brief biological background note. At most one or two per lesson. Motivation needed to follow the argument belongs in the main prose. |
| `::: {.callout-warning title="..."}` | One misconception worth interrupting the page for. Discuss ordinary limitations in the main prose where they matter. |
| `::: {.callout-note title="..."}` | A convention or assumption that does not deserve definition numbering. |
| `::: {.callout-tip title="Check: ..."}` | A verification habit: the small input that would catch a wrong implementation. |
| `::: {.callout-note collapse="true" title="Going deeper: ..."}` | Longer proofs, optional derivations, connections to the literature. |
| `::: {.callout-important}` | Reserved for data-use and reproducibility requirements. |

Use margin notes for history and optional connections that are not load-bearing. A margin note must be skippable without losing the thread.

Do not use `callout-note` with a "Definition:" title in place of a `def-` block; it does not number or cross-reference.

## Web style and branding

- The website uses University of Utah core colors for major UI: Utah Red `#BE0000`, black, white, and neutral greys.
- Accent colors are subordinate. Use them only for small elements such as callout borders, charts, or infographics. Do not use accent colors for headings, links, page backgrounds, or other major brand-bearing elements.
- Avoid blue or purple as stand-alone interface colors; they can read as other regional university colors.
- Headings use Montserrat Bold when available. Body copy uses Source Sans when available, with system sans-serif fallbacks.
- Use the website format deliberately: callouts, margin notes, collapsible code, and accessible figure captions can carry structure that a PDF cannot.

## Math notation

- Sequences: $x = x_1 x_2 \ldots x_n$, $|x| = n$.
- Score function $s(a,b)$; pick one symbol per matrix per lesson and stay consistent.
- $\log$ without base unless base matters.
- $\mathcal{O}(\cdot)$ for asymptotics; plain English when constants matter.
- Add shared symbols to `course_notes/appendices/notation.md` and define them locally at first use.

## Tables

- Markdown pipe tables for hand-written content.
- Computed tables: produce with `df.to_markdown(index=False)` inside a code chunk.
- Caption and label inline: `: My caption {#tbl-foo}`. Reference with `@tbl-foo`.
- Right-align numeric columns; left-align text.
- Keep tables narrow enough to render on mobile: no more than about five columns of normal text.

## Cross-references

Use Quarto labels everywhere; never "see the figure above." Standard prefixes: `@sec-`, `@fig-`, `@tbl-`, `@eq-`, plus the theorem-like prefixes above.

## Equations

Display equations get their own line. Number (with `{#eq-foo}`) only if referenced later. Always introduce a symbol the first time it appears in a lesson.

## Citations

Keep bibliography entries in `course_notes/references.bib`. Cite non-trivial historical or empirical claims where they occur. Credit books and resources actually consulted using front-matter `nocite` and a compact `#refs` div at the bottom:

```markdown
::: {#refs}
**References**
:::
```

The website renders the bibliography in small, readable type. Do not add a "Further reading" section or prose claiming exact chapter-by-chapter provenance. References give credit; they are not assigned reading. Only cite sources actually checked. Keep exact internal source notes in ignored working material when useful.

## Patterns to remove during editing

Remove promotional claims, vague attribution, and filler that does not explain the subject. Avoid buzzwords such as "leverage," "delve," "tapestry," and "seamless." Use "robust" only with a precise technical meaning. Replace inflated judgments such as "a deep insight" with the insight itself.

Do not force lists into groups of three, fragment connected reasoning into bullets, or stack callout boxes. Prefer periods and commas to repeated dashes. Avoid dramatic fragments and stock contrasts such as "not just X, but Y." Keep comparisons that distinguish quantities students might actually confuse.

Remove generic motivation essays and closing praise. Remove the filler that rules applied mechanically produce: a disclaimer after every result, statements of what a check does not prove, answers certified before the method computes them, and mid-paragraph pointers to other lessons. Keep concrete biological motivation, bridges that say where the argument stands, and short syntheses that help a reader follow the argument. Internal course logistics belong in the notes only when students need the information to act.

When editing LLM drafts, invoke the `humanizer` skill. Removing verbal habits is a final editing pass; it cannot supply missing mathematical reasoning.

## Stance on software

Code serves mathematical understanding. Mention real tools such as BWA, BLAST, and scanpy briefly when they put the pedagogical method in context; explain the relevant difference without turning the lesson into a software tutorial.

## Assignments and validation

Allow a student's preferred language unless a particular task needs a restriction. Python and scientific packages describe the lecture examples, not a universal submission requirement. Chris will rarely inspect code line by line; emphasize the reasoning, independently justified checks, actual outputs, and interpretation. Code should remain reproducible and explainable.

For an introductory problem set, aim for 3 or 4 problems and about one hour: small hand calculations or proofs, a compact coding or debugging task, and a concrete interpretation question. Let students use AI where permitted by the course policy. Ask for meaningful validation and a brief disclosure, not long process narratives. No take-home question is literally LLM-proof; unsupported claims of AI resistance do not belong in student material. Do not introduce new graded quizzes or assessment weights without a course-level decision.

## Required review before accepting a lesson

Review the reasoning before polishing individual sentences. Revise by rewriting: cut what does not help a reader understand before adding anything. A lesson is not ready if any applicable check below fails.

1. **Purpose:** Each substantial section has a concrete reason to exist, supplied locally or by its immediate context. The opening of the lesson connects the data to a biological question and a specific output. A reader can state the lesson's central idea in plain words and why it is needed.
2. **Dependencies:** New objects are motivated and defined before use. The reader can follow the main argument with the stated prerequisites; necessary explanations are not hidden in optional material.
3. **Reasoning and mathematical precision:** Important objectives, assumptions, models, recurrences, comparisons, and results have explicit mathematical statements where useful, with locally defined and consistent symbols. Every central equation has an interpretation. Derivations explain the consequential steps, and correctness and cost claims have justifications.
4. **Examples and figures:** At least one small instance shows the method making a decision that matters, and at least one figure shows the central idea at work. Prose, mathematics, code, and figures agree.
5. **Code:** The implementation reads as a summary of the mathematics and shows actual output. Executable chunks run under the conventions in `STYLE_CODE.md`.
6. **Scope:** Assumptions and limitations appear once, where they affect conclusions. A model's answer is not presented as an unqualified biological fact.
7. **Continuity:** Terminology, notation, voice, example choices, and formatting agree with connected lessons and the shared guides. Bridges explain the relationship between topics.
8. **Pacing and prose:** The material fits the meeting with time to work examples. Remove expendable topics and filler before removing explanations. Read the prose aloud for awkward sentences and ambiguous references.

Finally, check the required markup, citations, and cross-references. For changes to lesson content, build and inspect the rendered page using `site/AUTHORING.md`; code and figure changes also require the checks in their respective guides.
