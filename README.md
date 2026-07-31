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

The infinite-block child observes that both support counts are exactly affine
in `h`. It measures the two intercepts with one exact FFT per construction,
then evaluates arbitrarily large block widths and the closed-form
`h -> infinity` exponent limit for substantially deeper/larger CRT cases.

The parameter-frontier child drops the special restriction `s=2^K+1` and
searches every odd `s` coprime to 12 up to ten million, for digit depths
through 256. It evaluates the exact infinite-block formula on GPUs and
rechecks the leading candidates with arbitrary-precision integers.

The asymptotic-depth child extends this optimization through depth 2048 on a
dense logarithmic `s` grid, snaps the floating candidates to valid integers,
and exact-checks both the global leaders and fixed-depth checkpoints.

The deep-scaling child uses a two-million-point log grid through depth 8192.
All objective calculations stay in the log domain, and exact verification uses
bit-length logarithms so even thousand-digit construction parameters do not
overflow floating point.

The ultra-deep child extends the frontier to depth 32768 on an eight-million
point grid.  It derives every logarithmic digit-growth factor from the exact
integer recurrence, avoiding both underflow and cancellation, and constructs
the enormous admissible outer parameters with arbitrary-precision decimals.

The convergence-constant child pushes to depth 131072 with a sixteen-million
point grid and reports the scaled gap `d(2-C)`, which should approach a finite
constant if the optimized exponent converges to two at rate `1/d`.

The asymptotic-constant child reaches depth 524288 and compares the measured
scaled gap against `3 log(2) / -log(lambda/12)`, with `lambda` computed from
the exact characteristic polynomial of the improved carry automaton.

The million-depth child reaches depth 2097152.  Its GPU search uses the
dominant-root log-growth model calibrated against an exact recurrence value,
while the parent concurrently computes exact counts for every candidate depth.

The eight-million-depth child quadruples this range once more, using a
128-million-point grid and exact checkpoints through depth 8388608.
