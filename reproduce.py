#!/usr/bin/env python3
"""Exact and computational reproduction of arXiv:2607.27199.

The proof itself is exact integer arithmetic. CUDA is used for a neighboring
exhaustive gadget search, which is an independent stress test of the paper's
base-12 choice rather than part of the proof.
"""

from __future__ import annotations

import itertools
import json
import math
import multiprocessing as mp
import os
import queue
import sys
import time
from decimal import Decimal, localcontext
from fractions import Fraction
from typing import Iterable


W = frozenset({0, 1, 2, 4, 5, 9})
IMPROVED_W = frozenset({0, 1, 3, 4, 5, 8})
DELTA = frozenset(
    {-9, -8, -7, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 7, 8, 9}
)
BALANCED = tuple(range(-6, 6))
STATES = (frozenset({0}), frozenset({-1, 0}), frozenset({0, 1}))
M_EXPECTED = ((5, 3, 3), (4, 5, 3), (4, 4, 4))


def modular_sumset(a: Iterable[int], b: Iterable[int], modulus: int) -> set[int]:
    aa, bb = tuple(a), tuple(b)
    return {(x + y) % modulus for x in aa for y in bb}


def modular_diffset(a: Iterable[int], b: Iterable[int], modulus: int) -> set[int]:
    aa, bb = tuple(a), tuple(b)
    return {(x - y) % modulus for x in aa for y in bb}


def integer_sumset(a: Iterable[int], b: Iterable[int]) -> set[int]:
    aa, bb = tuple(a), tuple(b)
    return {x + y for x in aa for y in bb}


def integer_diffset(a: Iterable[int], b: Iterable[int]) -> set[int]:
    aa, bb = tuple(a), tuple(b)
    return {x - y for x in aa for y in bb}


def digit_set(j: int) -> set[int]:
    values = {0}
    place = 1
    for _ in range(j):
        values = {x + w * place for x in values for w in W}
        place *= 12
    return values


def recurrence_values(max_j: int) -> list[int]:
    vals = [1]
    if max_j:
        vals.append(11)
    while len(vals) <= max_j:
        vals.append(13 * vals[-1] - 16 * vals[-2])
    return vals


def phi(epsilon: int, carries: frozenset[int]) -> frozenset[int]:
    return frozenset(
        next_carry
        for next_carry in (-1, 0, 1)
        if any(epsilon + carry - 12 * next_carry in DELTA for carry in carries)
    )


def verify_digit_gadget() -> dict:
    sums = modular_sumset(W, W, 12)
    diffs = modular_diffset(W, W, 12)
    assert sums == set(range(12))
    assert diffs == set(range(12)) - {6}
    assert {x - y for x in W for y in W} == set(DELTA)

    transition = []
    for current in STATES:
        row = []
        for nxt in STATES:
            row.append(sum(phi(epsilon, current) == nxt for epsilon in BALANCED))
        assert all(
            not phi(epsilon, current) or phi(epsilon, current) in STATES
            for epsilon in BALANCED
        )
        transition.append(tuple(row))
    assert tuple(transition) == M_EXPECTED

    recurrence = recurrence_values(16)
    direct = []
    for j in range(5):
        y = digit_set(j)
        assert len(y) == 6**j
        if j:
            assert modular_sumset(y, y, 12**j) == set(range(12**j))
        t_direct = len(modular_diffset(y, y, 12**j))
        direct.append(t_direct)
        assert t_direct == recurrence[j]

    lam = (13 + math.sqrt(105)) / 2
    mu = (13 - math.sqrt(105)) / 2
    assert lam / 12 < 31 / 32
    assert (31 / 32) ** 22 < 1 / 2
    for j, exact in enumerate(recurrence[1:], 1):
        closed = (
            (math.sqrt(105) + 9) / (2 * math.sqrt(105)) * lam**j
            + (math.sqrt(105) - 9) / (2 * math.sqrt(105)) * mu**j
        )
        assert abs(closed - exact) < max(1e-7, exact * 1e-12)
        assert exact < lam**j

    return {
        "W": sorted(W),
        "modular_sum_size": len(sums),
        "modular_difference_size": len(diffs),
        "transition_matrix": transition,
        "direct_t_0_to_4": direct,
        "recurrence_t_0_to_16": recurrence,
        "lambda": lam,
        "lambda_over_12": lam / 12,
    }


def cyclic_basis(k_value: int) -> tuple[int, int, set[int], set[int]]:
    s = 2**k_value + 1
    q_outer = s * s
    m = (s - 1) // 2
    h = {x % q_outer for x in range(-m, m + 1)}
    v = {(j * s) % q_outer for j in range(s)}
    i_set = h | v
    b_set = i_set - {0}
    assert h & v == {0}
    assert len(i_set) == 2 * s - 1
    assert {-x % q_outer for x in i_set} == i_set
    assert modular_sumset(i_set, i_set, q_outer) == set(range(q_outer))
    bb_sum = modular_sumset(b_set, b_set, q_outer)
    bb_diff = modular_diffset(b_set, b_set, q_outer)
    outside = set(range(q_outer)) - i_set
    assert outside <= bb_sum
    assert outside <= bb_diff
    return s, q_outer, i_set, b_set


def crt_toy_check() -> dict:
    """Directly enumerate K=2 with d=2 instead of the theorem's d=88.

    Lemmas 2.5--2.7 do not depend on the large-d inequalities, so the reduced
    instance validates the exact CRT fiber formulas and integer lift directly.
    """

    k_value, d = 2, 2
    s, q_outer, i_set, b_set = cyclic_basis(k_value)
    rho = len(i_set)
    n = 12**d
    q_total = q_outer * n
    assert math.gcd(q_outer, n) == 1
    y = digit_set(d)
    t = len(modular_diffset(y, y, n))

    r_set = {
        x
        for x in range(q_total)
        if x % q_outer == 0 or (x % q_outer in b_set and x % n in y)
    }
    expected_r = n + (rho - 1) * 6**d
    assert len(r_set) == expected_r
    rsum_mod = modular_sumset(r_set, r_set, q_total)
    rdiff_mod = modular_diffset(r_set, r_set, q_total)
    expected_diff = rho * n + (q_outer - rho) * t
    assert len(rsum_mod) == q_total
    assert len(rdiff_mod) == expected_diff

    a_set = r_set | {x + q_total for x in r_set}
    sums = integer_sumset(a_set, a_set)
    diffs = integer_diffset(a_set, a_set)
    assert len(a_set) == 2 * expected_r
    assert len(sums) >= 3 * len(rsum_mod)
    assert len(diffs) <= 4 * len(rdiff_mod)

    sigma = Fraction(len(sums), len(a_set))
    delta = Fraction(len(diffs), len(a_set))
    c_value = math.log(float(sigma)) / math.log(float(delta))
    return {
        "K": k_value,
        "reduced_d": d,
        "Q": q_outer,
        "n": n,
        "|R|": len(r_set),
        "|R+R mod q|": len(rsum_mod),
        "|R-R mod q|": len(rdiff_mod),
        "difference_formula": expected_diff,
        "|A|": len(a_set),
        "|A+A|": len(sums),
        "|A-A|": len(diffs),
        "sigma": float(sigma),
        "delta": float(delta),
        "C": c_value,
        "sum_lift_ratio": len(sums) / len(rsum_mod),
        "diff_lift_ratio": len(diffs) / len(rdiff_mod),
    }


def theorem_certificates() -> list[dict]:
    rows = []
    for k_value in range(2, 66, 2):
        s = 2**k_value + 1
        q_outer = s * s
        rho = 2 * s - 1
        d = 22 * (k_value + 2)
        n = 12**d
        assert math.gcd(q_outer, n) == 1
        t = recurrence_values(d)[-1]
        alpha = Fraction((rho - 1) * 6**d, n)
        beta = Fraction(t, n)
        assert alpha == Fraction(1, 2 ** (21 * k_value + 43))
        assert alpha < Fraction(1, 2)
        assert beta < Fraction(1, 2 ** (k_value + 2))

        cardinality_a = 2 * (n + (rho - 1) * 6**d)
        sigma_lower = Fraction(3 * q_outer * n, cardinality_a)
        delta_upper = Fraction(
            4 * (rho * n + (q_outer - rho) * t), cardinality_a
        )
        assert sigma_lower > q_outer
        assert delta_upper < 5 * s - 2
        assert 5 * s - 2 < 8 * s
        theorem_bound = Fraction(2 * k_value, k_value + 3)
        with localcontext() as decimal_context:
            decimal_context.prec = 100
            log_s = Decimal(s).ln()
            c_cert_decimal = 2 * log_s / (log_s + 3 * Decimal(2).ln())
            theorem_bound_decimal = Decimal(theorem_bound.numerator) / Decimal(
                theorem_bound.denominator
            )
            assert c_cert_decimal > theorem_bound_decimal
        rows.append(
            {
                "K": k_value,
                "d": d,
                "digits_|A|": len(str(cardinality_a)),
                "beta_log2": math.log2(t) - d * math.log2(12),
                "C_certificate": float(c_cert_decimal),
                "C_certificate_decimal": str(c_cert_decimal),
                "theorem_bound": float(theorem_bound),
            }
        )
    return rows


def improved_theorem_certificates() -> tuple[dict, list[dict]]:
    automaton = carry_automaton_metrics(sorted(IMPROVED_W), 12)
    matrix = automaton["transition_matrix"]
    expected_matrix = [[7, 2, 2], [6, 4, 2], [6, 3, 3]]
    assert matrix == expected_matrix

    # Exact Collatz--Wielandt upper certificate:
    # M v <= r v for v=(1,53/48,53/48), r=183/16.
    vector = (Fraction(1), Fraction(53, 48), Fraction(53, 48))
    rational_radius = Fraction(183, 16)
    matrix_vector = tuple(
        sum(Fraction(matrix[row][col]) * vector[col] for col in range(3))
        for row in range(3)
    )
    assert all(
        matrix_vector[index] <= rational_radius * vector[index]
        for index in range(3)
    )
    normalized = rational_radius / 12
    assert normalized == Fraction(61, 64)
    assert normalized**15 < Fraction(1, 2)

    counts = automaton["t_0_to_8"]
    for j, count in enumerate(counts):
        assert count <= rational_radius**j
    for j in range(6):
        if len(counts) <= j + 3:
            break
        assert (
            counts[j + 3]
            == 14 * counts[j + 2] - 31 * counts[j + 1] + 18 * counts[j]
        )

    rows = []
    for k_value in range(2, 98, 2):
        d = 15 * (k_value + 2)
        s = 2**k_value + 1
        q_outer = s * s
        rho = 2 * s - 1
        alpha = Fraction(2 ** (k_value + 1), 2**d)
        beta_upper = normalized**d
        assert alpha < Fraction(1, 2)
        assert beta_upper < Fraction(1, 2 ** (k_value + 2))
        assert q_outer * beta_upper < Fraction(s, 2)
        rows.append(
            {
                "K": k_value,
                "d_improved": d,
                "d_paper": 22 * (k_value + 2),
                "depth_reduction_fraction": Fraction(7, 22),
                "alpha_log2": k_value + 1 - d,
                "beta_upper_log2": float(
                    d * math.log2(float(normalized))
                ),
            }
        )
    certificate = {
        "digits": sorted(IMPROVED_W),
        "transition_matrix": matrix,
        "t_0_to_8": counts,
        "characteristic_recurrence": "t[j+3]=14t[j+2]-31t[j+1]+18t[j]",
        "collatz_vector": [str(value) for value in vector],
        "rational_radius_bound": str(rational_radius),
        "normalized_bound": str(normalized),
        "normalized_bound_power_15": str(normalized**15),
        "computed_spectral_radius": automaton["spectral_radius"],
        "computed_normalized_growth": automaton["normalized_growth"],
    }
    return certificate, rows


def rotate_left(values, shift: int, bits: int, full_mask: int):
    if shift == 0:
        return values
    return ((values << shift) | (values >> (bits - shift))) & full_mask


def carry_automaton_metrics(digits: list[int], base: int) -> dict:
    """Build the reachable carry-subset automaton for an even base."""

    assert base % 2 == 0
    delta = {x - y for x in digits for y in digits}
    balanced = tuple(range(-(base // 2), base // 2))

    def transition(state: frozenset[int], epsilon: int) -> frozenset[int]:
        return frozenset(
            next_carry
            for next_carry in (-1, 0, 1)
            if any(
                epsilon + carry - base * next_carry in delta for carry in state
            )
        )

    initial = frozenset({0})
    states = [initial]
    cursor = 0
    while cursor < len(states):
        state = states[cursor]
        cursor += 1
        for epsilon in balanced:
            next_state = transition(state, epsilon)
            if next_state and next_state not in states:
                states.append(next_state)

    matrix = [[0 for _ in states] for _ in states]
    for row, state in enumerate(states):
        for epsilon in balanced:
            next_state = transition(state, epsilon)
            if next_state:
                matrix[row][states.index(next_state)] += 1

    # Perron iteration for the nonnegative transition matrix. Reachable
    # carry automata here have positive diagonal entries, avoiding periodicity.
    perron_vector = [1.0 for _ in states]
    spectral_radius = 0.0
    for _ in range(500):
        next_vector = [
            sum(matrix[row][col] * perron_vector[col] for col in range(len(states)))
            for row in range(len(states))
        ]
        spectral_radius = max(next_vector)
        assert spectral_radius > 0
        perron_vector = [value / spectral_radius for value in next_vector]
    counts = [1]
    vector = [1 if state == initial else 0 for state in states]
    for _ in range(8):
        vector = [
            sum(vector[row] * matrix[row][col] for row in range(len(states)))
            for col in range(len(states))
        ]
        counts.append(sum(vector))

    # Independently enumerate the first three digit products.
    direct_counts = []
    values = {0}
    for j in range(4):
        modulus = base**j
        direct_counts.append(
            len(modular_diffset(values, values, modulus)) if j else 1
        )
        assert direct_counts[-1] == counts[j]
        values = {x + digit * modulus for x in values for digit in digits}

    return {
        "states": [sorted(state) for state in states],
        "transition_matrix": matrix,
        "spectral_radius": spectral_radius,
        "normalized_growth": spectral_radius / base,
        "t_0_to_8": counts,
        "direct_t_0_to_3": direct_counts,
    }


def gpu_gadget_worker(local_rank: int, pod_index: int, result_queue) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    # Both pods independently repeat the full sweep, so the leader's log
    # contains all bases while the second node provides cross-node replication.
    base = 12 + 2 * local_rank
    target_size = base // 2
    full_mask = (1 << base) - 1
    total_masks = 1 << base
    chunk = 1 << 20
    additive_count = 0
    best_diff = base + 1
    best_count = 0
    best_mask = None
    best_masks: list[int] = []
    started = time.time()

    for start in range(0, total_masks, chunk):
        stop = min(total_masks, start + chunk)
        masks = torch.arange(start, stop, dtype=torch.int64, device=device)
        # Translation-normalize by requiring zero in W.
        keep = (masks & 1) == 1
        popcount = torch.zeros_like(masks)
        for bit in range(base):
            popcount += (masks >> bit) & 1
        masks = masks[keep & (popcount == target_size)]
        if masks.numel() == 0:
            continue

        sum_masks = torch.zeros_like(masks)
        diff_masks = torch.zeros_like(masks)
        for bit in range(base):
            selected = (masks >> bit) & 1
            sum_masks |= rotate_left(masks, bit, base, full_mask) * selected
            diff_masks |= rotate_left(
                masks, (base - bit) % base, base, full_mask
            ) * selected
        valid = masks[sum_masks == full_mask]
        valid_diffs = diff_masks[sum_masks == full_mask]
        additive_count += int(valid.numel())
        if valid.numel() == 0:
            continue

        diff_sizes = torch.zeros_like(valid_diffs)
        for bit in range(base):
            diff_sizes += (valid_diffs >> bit) & 1
        chunk_min = int(diff_sizes.min().item())
        if chunk_min <= best_diff:
            winners = valid[diff_sizes == chunk_min]
            count = int(winners.numel())
            example = int(winners.min().item())
            if chunk_min < best_diff:
                best_diff = chunk_min
                best_count = count
                best_mask = example
                best_masks = [int(value) for value in winners.cpu().tolist()]
            else:
                best_count += count
                best_mask = min(best_mask, example)
                best_masks.extend(int(value) for value in winners.cpu().tolist())

    torch.cuda.synchronize(device)
    example_set = [bit for bit in range(base) if (best_mask >> bit) & 1]
    spectra = []
    for mask in best_masks:
        digits = [bit for bit in range(base) if (mask >> bit) & 1]
        metrics = carry_automaton_metrics(digits, base)
        spectra.append((metrics["normalized_growth"], digits, metrics))
    spectra.sort(key=lambda item: (item[0], item[1]))
    best_spectrum = spectra[0]
    worst_spectrum = spectra[-1]
    result = {
        "pod_index": pod_index,
        "gpu": local_rank,
        "base": base,
        "subset_size": target_size,
        "normalized_candidates": math.comb(base - 1, target_size - 1),
        "sum_full_candidates": additive_count,
        "minimum_difference_size": best_diff,
        "minimizer_count": best_count,
        "lexicographic_mask_example": example_set,
        "best_carry_digits": best_spectrum[1],
        "best_carry_automaton": best_spectrum[2],
        "worst_normalized_growth": worst_spectrum[0],
        "distinct_normalized_growth_values": len(
            {round(item[0], 12) for item in spectra}
        ),
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("GPU_GADGET_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def run_gpu_search(pod_index: int) -> list[dict]:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(target=gpu_gadget_worker, args=(rank, pod_index, result_queue))
        for rank in range(visible)
    ]
    for process in processes:
        process.start()

    results = []
    while len(results) < len(processes):
        try:
            results.append(result_queue.get(timeout=30))
        except queue.Empty:
            failed = [p.exitcode for p in processes if p.exitcode not in (None, 0)]
            if failed:
                raise RuntimeError(f"GPU search worker failed: {failed}")
            print(
                f"GPU_SEARCH_PROGRESS pod={pod_index} "
                f"finished={len(results)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0
    return sorted(results, key=lambda row: row["base"])


def custom_digit_set(depth: int, digits: Iterable[int], base: int = 12) -> set[int]:
    values = {0}
    place = 1
    digit_tuple = tuple(digits)
    for _ in range(depth):
        values = {x + digit * place for x in values for digit in digit_tuple}
        place *= base
    return values


def build_crt_set(k_value: int, depth: int, digits: Iterable[int]) -> tuple:
    s, q_outer, i_set, b_set = cyclic_basis(k_value)
    n = 12**depth
    q_total = q_outer * n
    y = custom_digit_set(depth, digits)
    inverse_q = pow(q_outer, -1, n)
    r_set = {q_outer * inner for inner in range(n)}
    for outer in b_set:
        for inner in y:
            lift = outer + q_outer * (((inner - outer) * inverse_q) % n)
            r_set.add(lift)
    assert len(r_set) == n + (len(i_set) - 1) * len(y)
    return s, q_outer, len(i_set), n, q_total, y, r_set


def fft_support_counts(a_set: set[int], device) -> tuple[int, int, int]:
    import torch

    maximum = max(a_set)
    transform_length = 1 << (2 * maximum + 1).bit_length()
    indicator = torch.zeros(
        transform_length, dtype=torch.float64, device=device
    )
    indices = torch.tensor(sorted(a_set), dtype=torch.int64, device=device)
    indicator[indices] = 1.0
    frequency = torch.fft.rfft(indicator)
    sum_convolution = torch.fft.irfft(
        frequency * frequency, n=transform_length
    )
    diff_correlation = torch.fft.irfft(
        frequency * torch.conj(frequency), n=transform_length
    )
    sum_count = int((sum_convolution > 0.5).sum().item())
    diff_count = int((diff_correlation > 0.5).sum().item())
    del indicator, indices, frequency, sum_convolution, diff_correlation
    torch.cuda.empty_cache()
    return transform_length, sum_count, diff_count


def fft_comparison_worker(local_rank: int, pod_index: int, result_queue) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    k_value = 2 if local_rank < 4 else 4
    depth = local_rank % 4 + 1
    case_rows = []
    started = time.time()
    for label, digits in (("paper", W), ("improved", IMPROVED_W)):
        s, q_outer, rho, n, q_total, y, r_set = build_crt_set(
            k_value, depth, digits
        )
        a_set = r_set | {value + q_total for value in r_set}
        transform_length, sum_count, diff_count = fft_support_counts(
            a_set, device
        )
        t_value = len(modular_diffset(y, y, n))
        modular_diff_formula = rho * n + (q_outer - rho) * t_value
        assert sum_count >= 3 * q_total
        assert diff_count <= 4 * modular_diff_formula
        sigma = Fraction(sum_count, len(a_set))
        delta = Fraction(diff_count, len(a_set))
        c_value = math.log(float(sigma)) / math.log(float(delta))
        row = {
            "replica": pod_index,
            "gpu": local_rank,
            "gadget": label,
            "digits": sorted(digits),
            "K": k_value,
            "d": depth,
            "Q": q_outer,
            "n": n,
            "|Y-Y mod n|": t_value,
            "|R|": len(r_set),
            "|A|": len(a_set),
            "|A+A|": sum_count,
            "|A-A|": diff_count,
            "sigma": float(sigma),
            "delta": float(delta),
            "C": c_value,
            "fft_length": transform_length,
        }
        print("FFT_INTEGER_CASE " + json.dumps(row, sort_keys=True), flush=True)
        case_rows.append(row)
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "K": k_value,
        "d": depth,
        "paper": case_rows[0],
        "improved": case_rows[1],
        "delta_C": case_rows[1]["C"] - case_rows[0]["C"],
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("FFT_COMPARISON_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def run_fft_comparisons(pod_index: int) -> list[dict]:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=fft_comparison_worker,
            args=(rank, pod_index, result_queue),
        )
        for rank in range(visible)
    ]
    for process in processes:
        process.start()
    results = []
    while len(results) < len(processes):
        try:
            results.append(result_queue.get(timeout=30))
        except queue.Empty:
            failed = [p.exitcode for p in processes if p.exitcode not in (None, 0)]
            if failed:
                raise RuntimeError(f"FFT comparison worker failed: {failed}")
            print(
                f"FFT_PROGRESS replica={pod_index} "
                f"finished={len(results)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0
    return sorted(results, key=lambda row: (row["K"], row["d"]))


def main() -> None:
    pod_index = int(os.environ.get("JOB_COMPLETION_INDEX", "0"))
    print(
        "CONFIG "
        + json.dumps(
            {
                "paper": "arXiv:2607.27199v1",
                "pod_index": pod_index,
                "python": sys.version.split()[0],
                "fixed_command": "python -u reproduce.py",
            },
            sort_keys=True,
        ),
        flush=True,
    )

    gadget = verify_digit_gadget()
    print("DIGIT_GADGET " + json.dumps(gadget, sort_keys=True), flush=True)

    basis_rows = []
    for k_value in (2, 4, 6):
        s, q_outer, i_set, _ = cyclic_basis(k_value)
        basis_rows.append(
            {"K": k_value, "s": s, "Q": q_outer, "rho": len(i_set)}
        )
    print("CYCLIC_BASIS " + json.dumps(basis_rows, sort_keys=True), flush=True)

    toy = crt_toy_check()
    print("DIRECT_CRT_TOY " + json.dumps(toy, sort_keys=True), flush=True)

    certificates = theorem_certificates()
    selected = [
        row for row in certificates if row["K"] in {2, 4, 8, 16, 32, 64}
    ]
    print("THEOREM_CERTIFICATES " + json.dumps(selected, sort_keys=True), flush=True)

    improved_certificate, improved_rows = improved_theorem_certificates()
    print(
        "IMPROVED_THEOREM_CERTIFICATE "
        + json.dumps(improved_certificate, sort_keys=True),
        flush=True,
    )
    print(
        "IMPROVED_DEPTH_ROWS "
        + json.dumps(
            [
                row
                for row in improved_rows
                if row["K"] in {2, 4, 8, 16, 32, 64, 96}
            ],
            sort_keys=True,
            default=str,
        ),
        flush=True,
    )

    fft_rows = run_fft_comparisons(pod_index)
    print("FFT_COMPARISON_SUMMARY " + json.dumps(fft_rows, sort_keys=True), flush=True)

    summary = {
        "status": "PASS",
        "pod_index": pod_index,
        "exact_claim_groups_passed": 6,
        "theorem_K_values_checked": len(certificates),
        "improved_gadget": sorted(IMPROVED_W),
        "paper_normalized_growth": gadget["lambda_over_12"],
        "improved_normalized_growth": improved_certificate[
            "computed_normalized_growth"
        ],
        "paper_depth_multiplier": 22,
        "improved_depth_multiplier": 15,
        "fft_cases_checked": 2 * len(fft_rows),
        "fft_improved_C_wins": sum(
            row["delta_C"] > 0 for row in fft_rows
        ),
        "paper_W_difference_size": 11,
        "toy_C": toy["C"],
        "largest_C_certificate": certificates[-1]["C_certificate"],
        "largest_theorem_bound": certificates[-1]["theorem_bound"],
    }
    print("FINAL_REPRODUCTION_SUMMARY " + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
