---
published: false
title: "Problem Set 2: Global Alignment"
---

## Draft Role

This assignment develops a complete global-alignment calculation from the recurrence through a checked alignment. The required work uses Lesson 3's main text; optional extensions are not prerequisites.

## Planned Sections

1. Define the global-alignment recurrence (Needleman–Wunsch).
2. Fill a small DP table by hand.
3. Prove correctness with an optimal-substructure argument.
4. Trace an optimal path by hand, then adapt the lesson's score-table code to accept a supplied substitution-score function rather than fixed match/mismatch values. Supply a working traceback function that takes this table, the two sequences, the substitution-score function, and the gap score.
5. Write a short independent check that rescores the returned alignment and verifies that removing gaps recovers both inputs. State the check's rules in the assignment, including rejecting gap-gap columns. Include an empty-string boundary and a tiny tie case; accept any maximizing traceback.
6. Explain what an optimal alignment supports biologically, and disclose material AI/tool assistance.

Optional enrichment: implement traceback from scratch or compare how two scoring systems change the chosen alignment. Counting optimal paths and extended parameter-sensitivity analysis are not required.
