# Reproduction: optimal sumset–difference-set exponent

This repository independently checks the explicit construction in Lin and Li,
*Settling the Optimal Exponent Relating Sumsets and Difference Sets*
(arXiv:2607.27199).

The fixed experiment entry point is:

```bash
python -u reproduce.py
```

It verifies:

1. the base-12 digit gadget and its three-state carry automaton;
2. the recurrence for the modular difference count;
3. the symmetric additive basis in `Z / s^2 Z`;
4. the exact CRT sum/difference fiber counts;
5. the integer lifting inequalities by direct enumeration of a toy instance;
6. the theorem's quantitative bounds for a range of even `K`.

The baseline Kubernetes run also uses all 16 cluster GPUs for an independent
exhaustive search over half-density digit gadgets in neighboring moduli.

The carry-spectrum child goes further: it exhaustively finds every
minimum-difference, full-sum half-density gadget in even bases 12 through 26,
constructs the reachable carry-subset automaton for every minimizer, and ranks
its Perron growth factor normalized by the base. Both eight-GPU nodes repeat
the complete sweep, providing an independent cross-node replication while the
leader prints the complete result table.

The improved-construction child promotes the best gadget
`{0,1,3,4,5,8}`. It checks an exact rational Perron certificate reducing the
paper's digit depth from `22(K+2)` to `15(K+2)`, and uses GPU FFTs to count
the actual integer sumsets and difference sets for original/improved toy CRT
constructions.

The multi-block child generalizes the final integer lift to
`A_h = R + q{0,...,h-1}`. It checks the generalized lift constants and scans
`h=1,...,8` with exact GPU FFT support counts to measure which finite lift
maximizes the exponent ratio.
