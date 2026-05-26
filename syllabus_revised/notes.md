# MATH 4310 Durable Skills Infusion: Design Notes

Working notes for the Spring 2027 MATH 4310 redesign under the 2026 MBE-CTE Durable Skills Faculty Fellows program. Internal notes; the syllabus itself is written in the instructor's voice and does not cite the fellowship materials.

## The five durable skills MATH 4310 will explicitly teach

Drawn from the 8-skill University of Utah Durable Skills Framework. The fellowship asks for 3 to 5; we land on 5 because each is genuinely doing work and each has at least one piece of gradable evidence.

1. **Critical Thinking & Problem Solving.** Bedrock of the course. Visible in correctness proofs, complexity analysis, and trade-off evaluation. The Failure-Mode Reflection on every problem set is its primary new instrument; the standard math content is the rest.

2. **Information & Technology Literacy.** Natural fit for bioinformatics: tool selection, data provenance, model assumptions, complex-system thinking, AI-assisted-coding verification. Operationalized through the Tools/Data/AI use policy, the reproducibility requirements on the final project (dependency versions, data provenance, validation check), and tool/data documentation as a separate rubric criterion on the project.

3. **Communication.** Writing crisp technical prose for readers from another discipline. The Failure-Mode Reflection on every problem set is the primary instrument: students must write it so a peer who has not seen their code can follow the argument, and audience-awareness is one of the two anchors of its single 0/1/2 rubric (the other anchor being specificity, which carries the Critical Thinking signal). The final project's README is the longer instance, graded against the same anchor pattern.

4. **Collaboration & Teamwork.** Math/CS courses default to individual problem sets. Problem Set 5 (Week 6) is restructured as a triad set with assigned roles, contribution log, and individual addenda; the final project is also a triad effort. Group sizes are deliberately small (triads): pairs are too fragile, quads dilute responsibility.

5. **Career & Self-Development.** Operationalized as a reproducible, version-controlled final-project repository plus a one-page revision memo. The framework's "receiving feedback" clause lives in the two graded checkpoints; the "adapting / grit" clause lives in the memo's required obstacle-and-adaptation prompt. The framing in the syllabus emphasizes reproducible computational scholarship rather than employer-portfolio language, which keeps the design clear of the orientation deck's "not vocational training" boundary.

## Skills considered and not included (and why)

- **Creativity & Innovation.** A defensible fit for the course's "build from scratch" signature, but naming a fifth analytic-flavored skill alongside the four others would over-load the framework block. The Algorithm-redesign work that would have lived under Creativity is folded into the final project's rubric and into the Failure-Mode Reflection prompt.
- **Emotional & Cultural Competence.** The cross-disciplinary translation work (math, CS, biology speaking past each other) is real, but routing it through Communication ("explain to a peer from another discipline") is more honest and assessable.
- **Leadership & Professionalism.** Could be claimed via group-project role assignments, but it dilutes focus. Rolled into the Collaboration rubric instead.

## Assessment architecture

| Component | Weight | Was | DS evidence |
|---|---:|---:|---|
| 10 Problem Sets (9 individual + 1 triad set in Week 6), each with Failure-Mode Reflection | 35% | 40% | CT (primary), ITL, Comm |
| Project workshops/checkpoints (proposal + repo skeleton; preliminary results) | 5% | 0% | Career&SD, Collaboration, Comm, ITL, CT |
| Final Project (group of 3): repo + README + writeup + contribution log + revision memo | 30% | 0% | Career&SD, Collaboration, Comm, ITL, CT |
| Short Midterm Exam 1 | 15% | 25% | CT, ITL |
| Short Midterm Exam 2 | 15% | 35% | CT, ITL |

Net effect: the old cumulative final exam is replaced by a larger final project, and the old single midterm/final-exam structure is replaced by two shorter in-class midterms. The weekly technical-work category remains 40% total when the two project workshops are included. This makes the final project the cumulative assessment while preserving individual checks on the algorithmic and probabilistic core.

### Rubrics

**Failure-Mode Reflection rubric (single 0/1/2, two skills).** Used on every problem set's reflection. It is the only writing item on a weekly PS, so one rubric carries both Critical Thinking and Communication evidence:

- 0 = absent, generic, or unclear to a peer reader
- 1 = present but partial: assumption named but failure scenario or fix is generic, or the argument is hard to follow
- 2 = specific, concrete, and clear to a peer reader: assumption identified, scenario realistic, failure traced to a specific place in the solution, fix plausible, written so someone who has not seen the code can follow it

**README and triad-PS addendum rubric.** Same 0/1/2 anchors applied to longer items where the audience-awareness signal dominates: *accurate / audience-aware / no handwaving*.

**Triad Problem Set (Week 6):**
- Shared technical solution: graded once for the team
- Individual addendum: graded against the Communication 0/1/2 rubric
- Contribution log: used to adjust individual grades only when documented non-participation appears

**Final Project (6 criteria, mapped to 30% of course grade; the two checkpoints contribute 5% through the Problem Sets and Project Workshops category):**

The deliverable is a single Git repository with three artifacts: implementation (with reproducibility built in), a project README (audience-aware explanation + method/results/limitations + a Team and Contributions paragraph at the top), and a one-page revision memo. Six rubric criteria evaluate it:

| Criterion | DS evidence | Points | Anchor |
|---|---|---:|---|
| Algorithm implementation applied to a bioinformatics problem | CT (and CT + ITL) | 0 to 5 | correct logic, edge cases, readable from-scratch code, real data, appropriate assumptions, interpretable result |
| Reproducibility | ITL | 0 to 3 | dependencies, run instructions, data source/provenance, validation check |
| Tool / data / AI documentation | ITL | 0 to 2 | tools and sources named, verification of any AI or external help |
| Project README (audience-aware explanation, method, results, limitations) | Communication + CT | 0 to 5 | accurate, audience-aware, no handwaving, method/results/limitations all present |
| Team and Contributions section + checkpoint engagement | Collaboration | 0 to 3 | roles, coordination, code review, evidence of dealing with conflicts or rebalancing |
| Revision / adaptation memo | Career & SD | 0 to 3 | feedback received, revision made, obstacle and adaptation, next step |

One grading pass per project at end of semester. Checkpoint 1 (Week 11, proposal + repo skeleton) and Checkpoint 2 (Week 14, preliminary results) contribute through the Problem Sets and Project Workshops category, so feedback gets read before the final submission. The contribution log was originally a separate artifact; it is now a paragraph at the top of the README, which captures the same accountability signal without adding a deliverable.

### Workload defenses

- **Failure-Mode Reflection grading drift.** Realistic estimate is closer to 2 minutes per student per PS than 1 minute. Mitigate with 2 to 3 anchor exemplars per assignment posted in the rubric (a strong response, a borderline response, a generic response) so the calibration is shared and student-visible.
- **Paired PS in Week 6.** Triads, one-week scope, charter signed at start, contribution log required, individual addendum is the unit of grading for participation. Identical-submission risk neutralized by the addendum.
- **Final project grading at end of semester.** One pass per project against the 24-point rubric, no prose feedback unless requested. Use GitHub Classroom (free for educators) so submissions auto-clone into a single folder; if institutional preference dictates otherwise, accept zipped repos via Canvas with the same checklist.
- **No detached reflection essays.** Reflection is embedded into homework (Failure-Mode Reflection) and into the final project (revision memo), each tied to concrete technical work. The "100 reflection essays at finals week" failure mode does not arise.

## Compliance tracking

The PPM 6-100 / HB261 checklist is maintained separately in `requirements_checklist.md` so the student-facing syllabus does not include internal compliance rationale.

## Canvas Toolkit contributions (small, reusable artifacts)

For the cohort's shared toolkit. Each is one page, drops into any computational STEM course.

1. **0/1/2 Technical Communication Rubric** with two worked examples (a strong response, a thin response).
2. **Failure-Mode Reflection prompt template** with three anchor exemplars and the 0/1/2 rubric.
3. **Triad Problem Set template**: charter, role list (driver / verifier / scribe), contribution log, individual addendum prompt.
4. **AI / Tool / Data Use Disclosure Policy** for coding-heavy STEM assignments: when permitted, what to disclose, what to verify.
5. **Final Project rubric (8 criteria)** plus the GitHub Classroom setup guide and the README template.

## Marking durable-skills language in the syllabus

Durable-skills framing language in `syllabus.md` is wrapped in pandoc's `[text]{.underline}` syntax so it renders as real Word underlining in `syllabus.docx`. Search the markdown for `{.underline}` to see every DS-tagged passage. In the docx, the underlines are visually obvious and easy to highlight further with the Word highlighter tool.

## Open questions for the instructor

1. **Group size.** Triads recommended. Confirm or override.
2. **Final project platform.** GitHub Classroom is the cleanest path. If institutional preference is Canvas-only or U of U CHPC GitLab, the same rubric and checklists apply.
3. **Reflection essays as a separate component.** Currently absent. If you want a low-stakes 0/1 completion grade once per module (4 short reflections, ~150 words each), that is a clean addition; just say the word and I'll fold it in.
4. **U Career Success guest visit.** Optional one-period visit in Week 14 from a UCS career coach to help students frame their final-project repos for biotech or bioinformatics roles. Not graded. Worth keeping?
