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

The Kubernetes run also uses all 16 cluster GPUs for an independent exhaustive
search over half-density digit gadgets in moduli 10 through 25. The search is
split over two indexed pods with eight GPUs each.
