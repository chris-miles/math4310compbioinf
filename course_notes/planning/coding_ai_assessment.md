# Coding and AI Assessment Model

## Core Position

Coding should be required in MATH 4310, but the coding should be small, scaffolded, and tied directly to the mathematics.

The goal is not software engineering. The goal is for students to understand an algorithm well enough to:

- Translate a recurrence or model into code.
- Create small examples with known answers.
- Detect silent failures.
- Explain why an implementation matches the mathematical object.
- Use AI and software tools without trusting them blindly.

This is the right modern skillset for bioinformatics. Real bioinformatics work often uses complex tools written by other people. The professional question is not "Did I type every line myself?" The professional question is "Do I know enough to verify and interpret what this tool produced?"

## Student-Facing Principle

Suggested wording:

> You may use AI coding assistants, but AI output is unverified output. You are responsible for testing, explaining, and defending every submitted algorithm. A correct-looking answer from an AI tool is not evidence of correctness unless you can validate it against a known case, derivation, or independent check.

## What Students Should Code

Students should implement core ideas from scratch when implementation forces them to understand the method.

| Topic | Coding expectation |
|---|---|
| k-mers and scoring | Implement k-mer counting and toy log-odds calculations. |
| Global alignment | Implement DP table and traceback. |
| Local alignment | Implement Smith-Waterman; derive affine gaps or implement if time allows. |
| Seed-and-extend | Implement a tiny seed index or search demonstration. |
| de Bruijn assembly | Implement graph construction and simple path/contig logic on toy reads. |
| Markov chains | Implement sequence likelihood and likelihood-ratio classification. |
| HMM/Viterbi | Implement Viterbi in log space, preferably as the triad assignment. |
| Forward-backward | Implement a small forward algorithm or posterior table, or trace by hand plus partial code. |
| Phylogenetics | Implement UPGMA or Fitch scoring on small examples. |
| PCA/SVD | Use `numpy.linalg.svd` or equivalent; implement centering, projection, and interpretation checks. |
| FDR | Implement Benjamini-Hochberg from scratch. |
| t-SNE/UMAP | Do not implement from scratch; interpret and critique outputs. |
| ODEs | Use a solver or simple Euler method; focus on model behavior and assumptions. |

## What To Scaffold

Provide:

- Data loading.
- File parsing.
- Plotting.
- Starter tests.
- Function signatures.
- Toy datasets.
- Recurrence/model statement.
- High-level pseudocode when useful.

Do not provide:

- The exact assessed inner loop.
- Complete traceback logic for alignment assignments.
- Complete Viterbi implementation.
- Complete BH implementation.
- Written explanation of why the algorithm is correct.

The line is: scaffold the environment, not the thinking.

## Assignment Template

Each coding assignment should include these sections.

### 1. Mathematical Object

Students identify the recurrence, graph, model, matrix factorization, test procedure, or ODE system.

Example prompt:

> State the recurrence your code implements. Define the table entries and boundary conditions.

### 2. Implementation

Students write a compact implementation of the core method.

Example prompt:

> Implement `needleman_wunsch(x, y, score, gap)` so that it returns both the optimal score and one optimal alignment.

### 3. Verification

Students provide at least two checks.

Recommended checks:

- Empty or one-character input.
- A hand-computed toy case.
- A deterministic edge case.
- Symmetry check where appropriate.
- Comparison to a known answer from the assignment handout.
- A failure case showing where assumptions break.

Example prompt:

> Include two tests where the expected answer is known before running your code. For each test, explain why the expected answer is known.

### 4. Explanation

Students explain how code matches math.

Example prompt:

> In 5-8 sentences, explain how your loops correspond to the recurrence and how traceback recovers an alignment.

### 5. AI / Tool Disclosure

Students disclose assistance.

Example prompt:

> List any AI tools, libraries, websites, or classmates that materially helped with your implementation. For AI-assisted code, say what you used it for and how you verified the result.

### 6. Failure-Mode Reflection

Students identify an assumption and how the method degrades.

Example prompt:

> Identify one assumption your solution depends on. Describe a realistic biological or computational situation where that assumption fails, and say what you would change.

## Grading Template

Recommended rubric for a coding-heavy problem:

| Criterion | Points |
|---|---:|
| Mathematical recurrence/model stated correctly | 20 |
| Core implementation matches the method | 30 |
| Verification tests are meaningful and pass | 20 |
| Explanation connects code to math | 15 |
| AI/tool disclosure is complete | 5 |
| Failure-mode reflection is concrete | 10 |

For lighter coding problems, collapse this to:

- Math/model: 25%.
- Code: 30%.
- Verification: 20%.
- Explanation: 15%.
- Disclosure/reflection: 10%.

## AI Use Policy for Assignments

Allowed:

- Asking AI to explain a bug.
- Asking AI to draft helper code.
- Asking AI to suggest tests.
- Asking AI to translate pseudocode into a first draft, if disclosed and verified.
- Using AI to improve comments or README clarity, if the technical claims remain yours.

Not acceptable:

- Submitting AI-generated code you cannot explain.
- Submitting AI-generated proofs or derivations as your own reasoning.
- Using an AI tool to bypass the core algorithmic task.
- Claiming code is correct because an AI tool said so.
- Omitting AI/tool use that materially affected the submission.

Boundary examples:

- Acceptable: "I used ChatGPT to help debug my traceback loop. I verified it on the hand example from the assignment and added a tie-case test."
- Not acceptable: "ChatGPT wrote my Viterbi implementation and I submitted it because it ran without errors."

## Exam Connection

Because AI use is allowed outside class, exams should check individual understanding without tools.

Exam questions should ask students to:

- Fill a small DP table.
- Complete a recurrence.
- Find a bug in pseudocode.
- Interpret a tiny output table.
- Explain a boundary condition.
- Compute one step of Viterbi, forward, UPGMA, Fitch, PCA, or BH.

This keeps the course honest: students can use modern tools on homework, but still need enough fluency to verify output independently.

