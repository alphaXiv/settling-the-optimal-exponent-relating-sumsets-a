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
    sum_count = int((sum_convolution > 0.5).sum().item())
    del sum_convolution
    diff_correlation = torch.fft.irfft(
        frequency * torch.conj(frequency), n=transform_length
    )
    diff_count = int((diff_correlation > 0.5).sum().item())
    del indicator, indices, frequency, diff_correlation
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


def block_lift_worker(local_rank: int, pod_index: int, result_queue) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    k_value = 2 if local_rank < 4 else 4
    depth = local_rank % 4 + 1
    s, q_outer, rho, n, q_total, y, r_set = build_crt_set(
        k_value, depth, IMPROVED_W
    )
    t_value = len(modular_diffset(y, y, n))
    modular_diff_formula = rho * n + (q_outer - rho) * t_value
    rows = []
    started = time.time()
    for block_count in range(1, 9):
        a_set = {
            value + block * q_total
            for block in range(block_count)
            for value in r_set
        }
        assert len(a_set) == block_count * len(r_set)
        transform_length, sum_count, diff_count = fft_support_counts(
            a_set, device
        )
        assert sum_count >= (2 * block_count - 1) * q_total
        assert diff_count <= 2 * block_count * modular_diff_formula
        sigma = Fraction(sum_count, len(a_set))
        delta = Fraction(diff_count, len(a_set))
        c_value = math.log(float(sigma)) / math.log(float(delta))
        row = {
            "replica": pod_index,
            "gpu": local_rank,
            "K": k_value,
            "d": depth,
            "h": block_count,
            "Q": q_outer,
            "n": n,
            "|R|": len(r_set),
            "|A_h|": len(a_set),
            "|A_h+A_h|": sum_count,
            "|A_h-A_h|": diff_count,
            "sigma": float(sigma),
            "delta": float(delta),
            "C": c_value,
            "sum_lift_per_residue_lower": 2 * block_count - 1,
            "diff_lift_per_residue_upper": 2 * block_count,
            "fft_length": transform_length,
        }
        print("BLOCK_LIFT_CASE " + json.dumps(row, sort_keys=True), flush=True)
        rows.append(row)
    best = max(rows, key=lambda row: row["C"])
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "K": k_value,
        "d": depth,
        "best_h": best["h"],
        "best_C": best["C"],
        "C_at_h2": rows[1]["C"],
        "gain_over_h2": best["C"] - rows[1]["C"],
        "monotone_through_h8": all(
            rows[index]["C"] >= rows[index - 1]["C"]
            for index in range(1, len(rows))
        ),
        "rows": rows,
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("BLOCK_LIFT_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def run_block_lift_audit(pod_index: int) -> list[dict]:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=block_lift_worker,
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
                raise RuntimeError(f"block-lift worker failed: {failed}")
            print(
                f"BLOCK_LIFT_PROGRESS replica={pod_index} "
                f"finished={len(results)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0
    return sorted(results, key=lambda row: (row["K"], row["d"]))


INFINITE_BLOCK_CASES = (
    (2, 5),
    (2, 6),
    (4, 5),
    (6, 1),
    (6, 2),
    (6, 3),
    (8, 1),
    (8, 2),
)


def infinite_block_worker(local_rank: int, pod_index: int, result_queue) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    k_value, depth = INFINITE_BLOCK_CASES[local_rank]
    s, q_outer, rho, n, q_total, y, r_set = build_crt_set(
        k_value, depth, IMPROVED_W
    )
    improved_counts = carry_automaton_metrics(sorted(IMPROVED_W), 12)[
        "t_0_to_8"
    ]
    t_value = improved_counts[depth]
    if depth <= 4:
        assert t_value == len(modular_diffset(y, y, n))
    modular_diff_count = rho * n + (q_outer - rho) * t_value

    transform_length, sum_count_h1, diff_count_h1 = fft_support_counts(
        r_set, device
    )
    two_quotient_sum_fibers = sum_count_h1 - q_total
    two_quotient_diff_fibers = diff_count_h1 - modular_diff_count
    assert 0 <= two_quotient_sum_fibers <= q_total
    assert 0 <= two_quotient_diff_fibers <= modular_diff_count

    def affine_row(block_count: int) -> dict:
        sum_count = (
            q_total * (2 * block_count - 1) + two_quotient_sum_fibers
        )
        diff_count = (
            modular_diff_count * (2 * block_count - 1)
            + two_quotient_diff_fibers
        )
        cardinality = block_count * len(r_set)
        sigma = Fraction(sum_count, cardinality)
        delta = Fraction(diff_count, cardinality)
        return {
            "h": block_count,
            "|A_h|": cardinality,
            "|A_h+A_h|": sum_count,
            "|A_h-A_h|": diff_count,
            "sigma": float(sigma),
            "delta": float(delta),
            "C": math.log(float(sigma)) / math.log(float(delta)),
        }

    h_values = (1, 2, 4, 8, 16, 32, 64, 128, 1024, 1_000_000)
    rows = [affine_row(block_count) for block_count in h_values]
    sigma_limit = Fraction(2 * q_total, len(r_set))
    delta_limit = Fraction(2 * modular_diff_count, len(r_set))
    c_limit = math.log(float(sigma_limit)) / math.log(float(delta_limit))
    assert all(
        rows[index]["C"] >= rows[index - 1]["C"]
        for index in range(1, len(rows))
    )
    assert c_limit >= rows[-1]["C"]
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "K": k_value,
        "d": depth,
        "Q": q_outer,
        "n": n,
        "|R|": len(r_set),
        "q": q_total,
        "modular_difference_count": modular_diff_count,
        "two_quotient_sum_fibers": two_quotient_sum_fibers,
        "two_quotient_diff_fibers": two_quotient_diff_fibers,
        "fft_length": transform_length,
        "rows": rows,
        "sigma_limit": float(sigma_limit),
        "delta_limit": float(delta_limit),
        "C_limit": c_limit,
        "gain_limit_over_h2": c_limit - rows[1]["C"],
        "gap_h1e6_to_limit": c_limit - rows[-1]["C"],
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("INFINITE_BLOCK_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def run_infinite_block_audit(pod_index: int) -> list[dict]:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=infinite_block_worker,
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
                raise RuntimeError(f"infinite-block worker failed: {failed}")
            print(
                f"INFINITE_BLOCK_PROGRESS replica={pod_index} "
                f"finished={len(results)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0
    return sorted(results, key=lambda row: (row["K"], row["d"]))


def parameter_frontier_worker(local_rank: int, pod_index: int, result_queue) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    maximum_s = 10_000_000
    maximum_depth = 256
    candidates = torch.arange(
        5, maximum_s + 1, dtype=torch.int64, device=device
    )
    candidates = candidates[(candidates % 2 == 1) & (candidates % 3 != 0)]
    s_values = candidates[local_rank::8].to(torch.float64)
    q_values = s_values * s_values
    rho_values = 2 * s_values - 1
    top_rows = []
    beta_values = [1.0, 11.0 / 12.0, 125.0 / (12.0**2)]
    while len(beta_values) <= maximum_depth:
        beta_values.append(
            (14.0 / 12.0) * beta_values[-1]
            - (31.0 / (12.0**2)) * beta_values[-2]
            + (18.0 / (12.0**3)) * beta_values[-3]
        )
    started = time.time()
    for depth in range(1, maximum_depth + 1):
        beta = beta_values[depth]

        alpha = (2 * s_values - 2) * (2.0 ** (-depth))
        denominator = 1.0 + alpha
        sigma = 2 * q_values / denominator
        delta = (
            2
            * (rho_values + (q_values - rho_values) * beta)
            / denominator
        )
        c_values = torch.log(sigma) / torch.log(delta)
        values, indices = torch.topk(c_values, k=8)
        for value, index in zip(values.cpu().tolist(), indices.cpu().tolist()):
            s_integer = int(s_values[index].item())
            top_rows.append(
                {
                    "replica": pod_index,
                    "gpu": local_rank,
                    "s": s_integer,
                    "d": depth,
                    "C_float": value,
                    "sigma_float": float(sigma[index].item()),
                    "delta_float": float(delta[index].item()),
                    "log10_R": (
                        depth * math.log10(12)
                        + math.log10(float(denominator[index].item()))
                    ),
                }
            )
    top_rows.sort(key=lambda row: row["C_float"], reverse=True)
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "s_candidates": int(s_values.numel()),
        "depths": maximum_depth,
        "top": top_rows[:32],
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("PARAMETER_SHARD_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def improved_t_values(maximum_depth: int) -> list[int]:
    values = [1, 11, 125]
    while len(values) <= maximum_depth:
        values.append(
            14 * values[-1] - 31 * values[-2] + 18 * values[-3]
        )
    return values[: maximum_depth + 1]


def exact_frontier_row(row: dict, t_values: list[int]) -> dict:
    s = row["s"]
    depth = row["d"]
    q_outer = s * s
    rho = 2 * s - 1
    n = 12**depth
    r_size = n + (rho - 1) * 6**depth
    modular_diff_count = rho * n + (q_outer - rho) * t_values[depth]
    sigma = Fraction(2 * q_outer * n, r_size)
    delta = Fraction(2 * modular_diff_count, r_size)
    c_exact = math.log(float(sigma)) / math.log(float(delta))
    assert abs(c_exact - row["C_float"]) < 1e-10
    return {
        **row,
        "sigma_exact_float": float(sigma),
        "delta_exact_float": float(delta),
        "C_exact": c_exact,
        "|R|_digits": len(str(r_size)),
        "gcd_s_12": math.gcd(s, 12),
    }


def run_parameter_frontier(pod_index: int) -> list[dict]:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=parameter_frontier_worker,
            args=(rank, pod_index, result_queue),
        )
        for rank in range(visible)
    ]
    for process in processes:
        process.start()
    shard_results = []
    while len(shard_results) < len(processes):
        try:
            shard_results.append(result_queue.get(timeout=30))
        except queue.Empty:
            failed = [p.exitcode for p in processes if p.exitcode not in (None, 0)]
            if failed:
                raise RuntimeError(f"parameter-frontier worker failed: {failed}")
            print(
                f"PARAMETER_PROGRESS replica={pod_index} "
                f"finished={len(shard_results)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0

    combined = [
        row for shard in shard_results for row in shard["top"]
    ]
    combined.sort(key=lambda row: row["C_float"], reverse=True)
    t_values = improved_t_values(256)
    exact_top = [exact_frontier_row(row, t_values) for row in combined[:64]]
    print(
        "PARAMETER_FRONTIER_EXACT_TOP "
        + json.dumps(exact_top, sort_keys=True),
        flush=True,
    )
    return exact_top


ASYMPTOTIC_CHECKPOINTS = frozenset({512, 1024, 1536, 2048})


def asymptotic_frontier_worker(
    local_rank: int, pod_index: int, result_queue
) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    grid_size = 200_000
    log_grid = torch.linspace(
        math.log(5), 100.0, grid_size, dtype=torch.float64, device=device
    )
    log_s_values = log_grid[local_rank::8]
    s_values = torch.exp(log_s_values)
    q_values = s_values * s_values
    rho_values = 2 * s_values - 1
    beta_values = [1.0, 11.0 / 12.0, 125.0 / (12.0**2)]
    while len(beta_values) <= 2048:
        beta_values.append(
            (14.0 / 12.0) * beta_values[-1]
            - (31.0 / (12.0**2)) * beta_values[-2]
            + (18.0 / (12.0**3)) * beta_values[-3]
        )

    top_rows = []
    checkpoint_rows = []
    started = time.time()
    for depth in range(257, 2049):
        beta = beta_values[depth]
        log_alpha = (
            torch.log(2 * s_values - 2) - depth * math.log(2)
        )
        log_denominator = torch.nn.functional.softplus(log_alpha)
        denominator = torch.exp(log_denominator)
        sigma = 2 * q_values / denominator
        delta = (
            2
            * (rho_values + (q_values - rho_values) * beta)
            / denominator
        )
        c_values = torch.log(sigma) / torch.log(delta)
        values, indices = torch.topk(c_values, k=4)
        depth_rows = []
        for value, index in zip(values.cpu().tolist(), indices.cpu().tolist()):
            row = {
                "replica": pod_index,
                "gpu": local_rank,
                "log_s": float(log_s_values[index].item()),
                "s_float": float(s_values[index].item()),
                "d": depth,
                "C_float": value,
                "sigma_float": float(sigma[index].item()),
                "delta_float": float(delta[index].item()),
            }
            top_rows.append(row)
            depth_rows.append(row)
        if depth in ASYMPTOTIC_CHECKPOINTS:
            checkpoint_rows.append(depth_rows[0])

    top_rows.sort(key=lambda row: row["C_float"], reverse=True)
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "log_s_samples": int(log_s_values.numel()),
        "depth_count": 2048 - 256,
        "top": top_rows[:32],
        "checkpoints": checkpoint_rows,
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("ASYMPTOTIC_SHARD_RESULT " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def nearest_admissible_s(s_float: float) -> int:
    center = int(round(s_float))
    candidates = [
        value
        for offset in range(-12, 13)
        if (value := center + offset) >= 5
        and value % 2 == 1
        and value % 3 != 0
    ]
    assert candidates
    return min(candidates, key=lambda value: abs(value - s_float))


def exact_asymptotic_row(row: dict, t_values: list[int]) -> dict:
    s = nearest_admissible_s(row["s_float"])
    depth = row["d"]
    q_outer = s * s
    rho = 2 * s - 1
    n = 12**depth
    r_size = n + (rho - 1) * 6**depth
    modular_diff_count = rho * n + (q_outer - rho) * t_values[depth]
    sigma = Fraction(2 * q_outer * n, r_size)
    delta = Fraction(2 * modular_diff_count, r_size)
    c_exact = math.log(float(sigma)) / math.log(float(delta))
    assert abs(c_exact - row["C_float"]) < 1e-7
    return {
        **row,
        "s": s,
        "C_exact": c_exact,
        "sigma_exact_float": float(sigma),
        "delta_exact_float": float(delta),
        "gap_to_2": 2 - c_exact,
        "|R|_digits": len(str(r_size)),
        "gcd_s_12": math.gcd(s, 12),
    }


def run_asymptotic_frontier(pod_index: int) -> dict:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=asymptotic_frontier_worker,
            args=(rank, pod_index, result_queue),
        )
        for rank in range(visible)
    ]
    for process in processes:
        process.start()
    shards = []
    while len(shards) < len(processes):
        try:
            shards.append(result_queue.get(timeout=30))
        except queue.Empty:
            failed = [p.exitcode for p in processes if p.exitcode not in (None, 0)]
            if failed:
                raise RuntimeError(f"asymptotic worker failed: {failed}")
            print(
                f"ASYMPTOTIC_PROGRESS replica={pod_index} "
                f"finished={len(shards)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0

    combined_top = [row for shard in shards for row in shard["top"]]
    combined_top.sort(key=lambda row: row["C_float"], reverse=True)
    combined_checkpoints = [
        row for shard in shards for row in shard["checkpoints"]
    ]
    checkpoint_winners = []
    for depth in sorted(ASYMPTOTIC_CHECKPOINTS):
        checkpoint_winners.append(
            max(
                (row for row in combined_checkpoints if row["d"] == depth),
                key=lambda row: row["C_float"],
            )
        )
    t_values = improved_t_values(2048)
    result = {
        "top": [
            exact_asymptotic_row(row, t_values)
            for row in combined_top[:64]
        ],
        "checkpoints": [
            exact_asymptotic_row(row, t_values)
            for row in checkpoint_winners
        ],
    }
    print(
        "ASYMPTOTIC_FRONTIER_EXACT " + json.dumps(result, sort_keys=True),
        flush=True,
    )
    return result


DEEP_CHECKPOINTS = frozenset({65536, 98304, 131072})


def improved_log_beta_values(maximum_depth: int) -> list[float]:
    """Compute log(t_d/12**d) without a floating recurrence."""
    values = [0.0] * (maximum_depth + 1)
    initial = [1, 11, 125]
    for depth, count in enumerate(initial):
        values[depth] = log_bigint(count) - depth * math.log(12)
    previous = initial
    for depth in range(3, maximum_depth + 1):
        count = 14 * previous[-1] - 31 * previous[-2] + 18 * previous[-3]
        values[depth] = log_bigint(count) - depth * math.log(12)
        previous = [previous[-2], previous[-1], count]
    return values


def improved_t_at_depths(depths: set[int]) -> dict[int, int]:
    """Evaluate the exact recurrence while retaining only requested depths."""
    wanted = set(depths)
    result = {}
    previous = [1, 11, 125]
    for depth, count in enumerate(previous):
        if depth in wanted:
            result[depth] = count
    for depth in range(3, max(wanted) + 1):
        count = 14 * previous[-1] - 31 * previous[-2] + 18 * previous[-3]
        if depth in wanted:
            result[depth] = count
        previous = [previous[-2], previous[-1], count]
    assert result.keys() == wanted
    return result


def deep_scaling_worker(
    local_rank: int,
    pod_index: int,
    log_beta_values: list[float],
    result_queue,
) -> None:
    import torch

    torch.cuda.set_device(local_rank)
    device = torch.device("cuda", local_rank)
    grid_size = 16_000_000
    full_grid = torch.linspace(
        1400.0, 6600.0, grid_size, dtype=torch.float64, device=device
    )
    log_s_values = full_grid[local_rank::8]

    best_rows = []
    checkpoints = []
    log_two = math.log(2)
    started = time.time()
    for depth in range(32769, 131073):
        log_beta = log_beta_values[depth]
        # At log(s)>=90, replacing log(2s-2) by log(2s) and
        # log(s^2-2s+1) by 2log(s) changes less than 1e-38.
        log_alpha = log_two + log_s_values - depth * log_two
        log_denominator = torch.nn.functional.softplus(log_alpha)
        log_sigma = log_two + 2 * log_s_values - log_denominator
        log_inner = torch.logaddexp(
            log_two + log_s_values,
            2 * log_s_values + log_beta,
        )
        log_delta = log_two + log_inner - log_denominator
        c_values = log_sigma / log_delta
        value, index = torch.max(c_values, dim=0)
        row = {
            "replica": pod_index,
            "gpu": local_rank,
            "log_s": float(log_s_values[index].item()),
            "d": depth,
            "C_float": float(value.item()),
            "log_sigma_float": float(log_sigma[index].item()),
            "log_delta_float": float(log_delta[index].item()),
        }
        best_rows.append(row)
        if depth in DEEP_CHECKPOINTS:
            checkpoints.append(row)

    best_rows.sort(key=lambda row: row["C_float"], reverse=True)
    result = {
        "replica": pod_index,
        "gpu": local_rank,
        "log_s_samples": int(log_s_values.numel()),
        "depth_count": 131072 - 32768,
        "top": best_rows[:16],
        "checkpoints": checkpoints,
        "seconds": time.time() - started,
        "device": torch.cuda.get_device_name(local_rank),
    }
    print("DEEP_SCALING_SHARD " + json.dumps(result, sort_keys=True), flush=True)
    result_queue.put(result)


def log_bigint(value: int) -> float:
    assert value > 0
    shift = max(0, value.bit_length() - 53)
    return math.log(value >> shift) + shift * math.log(2)


def exact_deep_row(row: dict, t_values: list[int]) -> dict:
    with localcontext() as context:
        context.prec = 100
        s = int(Decimal.from_float(row["log_s"]).exp())
    while math.gcd(s, 12) != 1:
        s += 1
    depth = row["d"]
    q_outer = s * s
    rho = 2 * s - 1
    n = 12**depth
    r_size = n + (rho - 1) * 6**depth
    modular_diff_count = rho * n + (q_outer - rho) * t_values[depth]
    log_sigma = (
        math.log(2)
        + log_bigint(q_outer)
        + log_bigint(n)
        - log_bigint(r_size)
    )
    log_delta = (
        math.log(2)
        + log_bigint(modular_diff_count)
        - log_bigint(r_size)
    )
    c_exact = log_sigma / log_delta
    assert abs(c_exact - row["C_float"]) < 1e-7
    return {
        **row,
        "s": s,
        "C_exact": c_exact,
        "gap_to_2": 2 - c_exact,
        "scaled_gap_d_times_2_minus_C": depth * (2 - c_exact),
        "log_sigma_exact": log_sigma,
        "log_delta_exact": log_delta,
        "|R|_digits": int(log_bigint(r_size) / math.log(10)) + 1,
        "gcd_s_12": math.gcd(s, 12),
    }


def run_deep_scaling(pod_index: int) -> dict:
    import torch

    visible = torch.cuda.device_count()
    assert visible == 8, f"manifest promised 8 visible GPUs, found {visible}"
    print(
        f"DEEP_SCALING_PRECOMPUTE replica={pod_index} max_depth=131072",
        flush=True,
    )
    log_beta_values = improved_log_beta_values(131072)
    ctx = mp.get_context("spawn")
    result_queue = ctx.Queue()
    processes = [
        ctx.Process(
            target=deep_scaling_worker,
            args=(rank, pod_index, log_beta_values, result_queue),
        )
        for rank in range(visible)
    ]
    for process in processes:
        process.start()
    shards = []
    while len(shards) < len(processes):
        try:
            shards.append(result_queue.get(timeout=30))
        except queue.Empty:
            failed = [p.exitcode for p in processes if p.exitcode not in (None, 0)]
            if failed:
                raise RuntimeError(f"deep-scaling worker failed: {failed}")
            print(
                f"DEEP_SCALING_PROGRESS replica={pod_index} "
                f"finished={len(shards)}/{len(processes)}",
                flush=True,
            )
    for process in processes:
        process.join()
        assert process.exitcode == 0

    combined_top = [row for shard in shards for row in shard["top"]]
    combined_top.sort(key=lambda row: row["C_float"], reverse=True)
    combined_checkpoints = [
        row for shard in shards for row in shard["checkpoints"]
    ]
    checkpoint_winners = [
        max(
            (row for row in combined_checkpoints if row["d"] == depth),
            key=lambda row: row["C_float"],
        )
        for depth in sorted(DEEP_CHECKPOINTS)
    ]
    exact_depths = {
        row["d"] for row in combined_top[:64] + checkpoint_winners
    }
    t_values = improved_t_at_depths(exact_depths)
    result = {
        "top": [
            exact_deep_row(row, t_values) for row in combined_top[:64]
        ],
        "checkpoints": [
            exact_deep_row(row, t_values) for row in checkpoint_winners
        ],
    }
    print("DEEP_SCALING_EXACT " + json.dumps(result, sort_keys=True), flush=True)
    return result


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

    deep_scaling = run_deep_scaling(pod_index)
    print(
        "DEEP_SCALING_SUMMARY "
        + json.dumps(deep_scaling, sort_keys=True),
        flush=True,
    )

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
        "deep_d_max": 131072,
        "deep_log_s_max": 6600.0,
        "deep_grid_size": 16_000_000,
        "deep_best": deep_scaling["top"][0],
        "deep_checkpoints": deep_scaling["checkpoints"],
        "paper_W_difference_size": 11,
        "toy_C": toy["C"],
        "largest_C_certificate": certificates[-1]["C_certificate"],
        "largest_theorem_bound": certificates[-1]["theorem_bound"],
    }
    print("FINAL_REPRODUCTION_SUMMARY " + json.dumps(summary, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
