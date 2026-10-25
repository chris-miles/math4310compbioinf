---
published: false
title: "Problem Set 7: Evolutionary Distances and Trees"
subtitle: "Lessons 13 and 14 · About 80 minutes"
nocite: |
  @durbin1998, @compeau2015, @felsenstein1978
---

Submit a short write-up with your calculations, actual program outputs, and explanations, plus the code needed to reproduce those outputs. The problems total 100 points. Use the due date posted in Canvas.

Python, R, Julia, MATLAB, or another language is fine. Lesson 13's `upgma` function and Lesson 14's `fitch_site` function may be copied, translated, or adapted. Library tree-building routines are not a substitute for the hand calculations below. You may use AI to help write or debug code. Briefly disclose that help and state how you checked the result. The [AI policy](../course-info.html#tools-data-and-ai-use) applies.

All sequences and distances in this set are synthetic. Distances are in substitutions per site unless a problem states other units.

## 1. UPGMA when the distances are not ultrametric (25 points, about 22 minutes)

Five sampled sequences give the following estimated distances, in units of 0.01 substitutions per site.

|   | A | B | C | D | E |
|---|---:|---:|---:|---:|---:|
| A | 0 | 4 | 14 | 20 | 36 |
| B | 4 | 0 | 14 | 20 | 36 |
| C | 14 | 14 | 0 | 32 | 36 |
| D | 20 | 20 | 32 | 0 | 25 |
| E | 36 | 36 | 36 | 25 | 0 |

a. Before running UPGMA, decide whether the matrix is ultrametric (@def-ultrametric). If it is not, give one triple of samples that violates the condition. Predict the first two merges.

b. Run UPGMA by hand with the size-weighted update from Lesson 13. For each merge, record the two clusters, the merge distance, the parent height, and both branch lengths below the parent. Check your merge history with Lesson 13's `upgma` function or your own implementation.

c. Compute all ten tip-to-tip path distances in your tree and tabulate them beside the input distances. Which distances does the tree reproduce exactly? Which pair has the largest absolute misfit? Relate the misfits to your answer in part (a).

d. In @exm-upgma-four-samples, the choice of averaging rule changed a merge distance but not the tree. Here it can do more. A classmate replaces the size-weighted update with a simple average at every merge: when clusters $X$ and $Y$ merge, the distance to another cluster $Z$ becomes $\bigl(d(X,Z)+d(Y,Z)\bigr)/2$. Find the first merge at which this changes a computed distance and the first merge at which it changes the topology. Write both trees in nested-parenthesis form. What does the disagreement tell a biologist about the evidence for where D belongs?

**Evidence to submit:** the violating triple, a merge table with heights and branch lengths, the comparison table of input and tree distances, and both trees.

## 2. The clock, the root, and four tips (20 points, about 15 minutes)

Part (c) is adapted from @durbin1998 and @compeau2015.

a. Choose nonnegative branch lengths for the rooted tree $((A,B),C)$ so that UPGMA, given the exact path distances, first joins a pair other than A and B. Certify your example: give the branch lengths, the three distances, and UPGMA's merges. Then write the terminal branch lengths as $a,b,c$ and the length of the edge from the A–B ancestor to the root as $t$. State an inequality on these lengths under which UPGMA joins the wrong pair first.

b. Two classmates discuss your example. One says: "With three tips there is only one unrooted tree, so UPGMA cannot get the tree wrong." The other says: "UPGMA gave the wrong answer to Lesson 13's core question." Say precisely what is correct and what is incorrect in each statement.

c. Consider an unrooted binary tree on four tips with split $ab\mid cd$: tips $a$ and $b$ attach to one internal node, $c$ and $d$ attach to another, and the edge joining the two internal nodes has length $\beta\geq0$. With additive path distances, prove that
$$
d_{ac}+d_{bd}=d_{ad}+d_{bc}=d_{ab}+d_{cd}+2\beta .
$$
Explain how this identity can identify the split from a distance matrix without assuming a molecular clock. You will use it in Problem 4.

**Evidence to submit:** your three-tip tree with its distances and merges, the inequality, the assessment of both statements, and the proof.

## 3. Diagnose a Fitch implementation (25 points, about 20 minutes)

This version of Lesson 14's `fitch_site` has one changed passage in the postorder step.

```python
def fitch_site(children, root, leaf_states):
    """Score one site and return one assignment; this version has a defect."""
    sets = {}

    def postorder(node):
        if node not in children:
            sets[node] = {leaf_states[node]}
            return 0
        left, right = children[node]
        score = postorder(left) + postorder(right)
        shared = sets[left] & sets[right]
        sets[node] = sets[left] | sets[right]
        if shared:
            return score
        return score + 1

    score = postorder(root)
    assigned = {root: min(sets[root])}

    def preorder(node):
        if node not in children:
            return
        for child in children[node]:
            if assigned[node] in sets[child]:
                assigned[child] = assigned[node]
            else:
                assigned[child] = min(sets[child])
            preorder(child)

    preorder(root)
    return score, sets, assigned
```

a. Before running it, predict whether this version reproduces Lesson 14's results: site scores $[1,1]$ on $((a,b),(c,d))$ for the two-site alignment, and totals 2, 3, and 3 for the three splits. Then run it and report the results.

b. Identify the change from Lesson 14's recurrence. Which statement fails for the sets this version stores: the attachment cost @eq-fitch-attachment or @thm-fitch-correctness?

c. Find a rooted binary tree with as few leaves as possible, and one site, for which this version returns the wrong score. Certify the true minimum by exhaustive search over the internal states, by hand or in code. Then argue that no tree with three or fewer leaves exposes the defect, and that the four-leaf shape $((a,b),(c,d))$ does not expose it either.

d. Show that this version's score is never larger than the true minimum. Hint: compare the stored sets with the correct Fitch sets node by node. Use this fact to explain why counting the changed edges in the returned assignment, as Lesson 14's code does, exposes every input on which the score is wrong.

**Evidence to submit:** the prediction and outputs, the identified change, the exposing tree and site with the exhaustive minimum, the minimality argument, and the inequality argument.

## 4. Long branches and three methods (30 points, about 25 minutes)

Adapted from @felsenstein1978.

The following 1,000-site alignment was simulated under the Jukes–Cantor model of Lesson 13 on the unrooted tree with split $AB\mid CD$. The terminal branches have lengths 0.5 to A, 0.05 to B, 0.5 to C, and 0.05 to D, and the internal edge has length $\beta=0.05$ substitutions per site. Lineages A and C therefore change ten times as fast as B and D.

The table records how many columns show each pattern of shared bases. A vertical bar separates groups of taxa with the same base, and different groups have different bases. For example, AB | CD means A and B share one base while C and D share another; BCD | A means only A differs; AB | C | D means A and B agree while C and D each have a different base.

| Pattern | Columns | Pattern | Columns | Pattern | Columns |
|---|---:|---|---:|---|---:|
| ABCD | 343 | AB \| CD | 31 | AD \| B \| C | 14 |
| BCD \| A | 204 | AC \| BD | 44 | BC \| A \| D | 10 |
| ACD \| B | 22 | AD \| BC | 6 | BD \| A \| C | 79 |
| ABD \| C | 192 | AB \| C \| D | 9 | CD \| A \| B | 20 |
| ABC \| D | 21 | AC \| B \| D | 1 | A \| B \| C \| D | 4 |

a. Before computing, predict which of the three unrooted splits parsimony will prefer, and explain why the branch lengths suggest it.

b. Show that the Fitch score of a single column differs among the three splits only for the patterns AB | CD, AC | BD, and AD | BC. You may score one representative column of each pattern on each split with `fitch_site`. Use the table to compute the total parsimony score of each split.

c. From the table, compute the six pairwise mismatch fractions, then the Jukes–Cantor distances. Apply the identity from Problem 2(c) to the raw fractions and to the corrected distances. Which split does each favor? Estimate $\beta$ from the corrected distances and compare it with 0.05. Run UPGMA on the corrected distances and report its first merge.

d. Which of the four analyses (parsimony, UPGMA on corrected distances, the four-tip identity on raw fractions, and the identity on corrected distances) recovers the split used in the simulation? For each failure, name the assumption responsible. If the alignment instead had 10,000 sites simulated from the same tree, which conclusions would you expect to change, and which would become more firmly held? Base your answer on the pattern proportions, not on a new simulation.

**Evidence to submit:** the prediction, the pattern-by-split scores and totals, the mismatch fractions, corrected distances, both four-tip comparisons, the UPGMA first merge, and the interpretation.

## Submission check

Include enough output to reproduce each conclusion without opening the code. A compact table or pasted text is preferable to screenshots. End with two or three sentences naming any code, tools, AI assistance, or people you used and describing one check you performed independently. State "No outside assistance" if applicable.

Credit follows the evidence: hand calculations, constructions whose stated property you verified, checks with independently known answers, and conclusions about trees limited to what the distances, characters, and model assumptions support.

::: {#refs}
**References**
:::
