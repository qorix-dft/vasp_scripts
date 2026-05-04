#!/usr/bin/env python3
"""
Lightweight comparison for two Phonopy band.yaml files.

This intentionally reads only metadata, atom positions, and frequencies. That is
enough to catch common failures before doing expensive eigenvector-overlap work.
"""

from __future__ import annotations

import argparse
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np


NUMBER_RE = re.compile(r"[-+]?\d+(?:\.\d*)?(?:[Ee][-+]?\d+)?")


def numbers(text: str) -> list[float]:
    return [float(x) for x in NUMBER_RE.findall(text)]


def parse_band_yaml(path: Path) -> dict:
    data = {
        "path": path,
        "natom": None,
        "lattice": [],
        "symbols": [],
        "frac_positions": [],
        "q_positions": [],
        "freqs_by_q": [],
    }

    in_lattice = False
    lattice_left = 0
    in_points = False
    current_symbol = None
    current_freqs = None

    with path.open(encoding="utf-8") as handle:
        for line in handle:
            stripped = line.strip()

            if stripped.startswith("natom:"):
                data["natom"] = int(stripped.split(":", 1)[1])
            elif stripped.startswith("lattice:"):
                in_lattice = True
                lattice_left = 3
            elif in_lattice and lattice_left:
                values = numbers(stripped.split("#", 1)[0])
                if len(values) >= 3:
                    data["lattice"].append(values[:3])
                    lattice_left -= 1
                if lattice_left == 0:
                    in_lattice = False
            elif stripped == "points:":
                in_points = True
            elif stripped == "phonon:":
                in_points = False
            elif in_points and stripped.startswith("- symbol:"):
                current_symbol = stripped.split()[2]
            elif in_points and stripped.startswith("coordinates:"):
                data["symbols"].append(current_symbol)
                data["frac_positions"].append(numbers(stripped)[:3])
            elif stripped.startswith("- q-position:"):
                data["q_positions"].append(numbers(stripped)[:3])
                current_freqs = []
                data["freqs_by_q"].append(current_freqs)
            elif current_freqs is not None and stripped.startswith("frequency:"):
                current_freqs.append(float(stripped.split(":", 1)[1].split()[0]))

    data["lattice"] = np.array(data["lattice"], dtype=float)
    data["frac_positions"] = np.array(data["frac_positions"], dtype=float)
    data["cart_positions"] = data["frac_positions"] @ data["lattice"]
    return data


def minimum_image(lattice: np.ndarray, dr_cart: np.ndarray) -> np.ndarray:
    dr_frac = dr_cart @ np.linalg.inv(lattice)
    dr_frac -= np.round(dr_frac)
    return dr_frac @ lattice


def summarize_frequencies(label: str, freqs_by_q: list[list[float]]) -> None:
    for q_index, freqs in enumerate(freqs_by_q):
        sorted_freqs = np.sort(np.array(freqs, dtype=float))
        n_modes = sorted_freqs.size
        negative = int(np.sum(sorted_freqs < -0.01))
        strongly_negative = int(np.sum(sorted_freqs < -0.1))
        near_zero = int(np.sum(np.abs(sorted_freqs) <= 0.1))
        print(f"{label} q{q_index}: modes={n_modes}")
        print(
            "  "
            f"min={sorted_freqs[0]:.6f} THz, "
            f"median={np.median(sorted_freqs):.6f} THz, "
            f"max={sorted_freqs[-1]:.6f} THz"
        )
        print(
            "  "
            f"negative<-0.01={negative}, "
            f"negative<-0.1={strongly_negative}, "
            f"|freq|<=0.1={near_zero}"
        )
        print("  first10=" + ", ".join(f"{x:.4f}" for x in sorted_freqs[:10]))


def compare_position_embedding(small: dict, large: dict, tolerance: float) -> None:
    used_large = np.zeros(len(large["symbols"]), dtype=bool)
    distances = []
    unmatched = []
    pairs = []

    for small_index, (symbol, small_pos) in enumerate(zip(small["symbols"], small["cart_positions"])):
        candidates = [
            large_index
            for large_index, large_symbol in enumerate(large["symbols"])
            if large_symbol == symbol and not used_large[large_index]
        ]
        if not candidates:
            unmatched.append((small_index, symbol, math.inf))
            continue

        dr = minimum_image(large["lattice"], large["cart_positions"][candidates] - small_pos[None, :])
        candidate_distances = np.sqrt(np.einsum("ij,ij->i", dr, dr))
        best_local = int(np.argmin(candidate_distances))
        best_distance = float(candidate_distances[best_local])
        best_large = candidates[best_local]

        if best_distance <= tolerance:
            used_large[best_large] = True
            distances.append(best_distance)
            pairs.append((small_index + 1, best_large + 1, best_distance))
        else:
            unmatched.append((small_index + 1, symbol, best_distance))

    print("Position embedding:")
    print(
        "  "
        f"matched={len(distances)}/{len(small['symbols'])}, "
        f"unmatched={len(unmatched)}, "
        f"unused_large_atoms={int(np.sum(~used_large))}"
    )
    if distances:
        d = np.array(distances)
        print(
            "  "
            f"distance_A min={np.min(d):.8f}, "
            f"mean={np.mean(d):.8f}, "
            f"p95={np.quantile(d, 0.95):.8f}, "
            f"max={np.max(d):.8f}"
        )
        print("  first10 small->large=" + ", ".join(f"{i}->{j} ({dist:.4f} A)" for i, j, dist in pairs[:10]))
    if unmatched:
        print("  first_unmatched=" + repr(unmatched[:5]))


def nearest_frequency_report(small_freqs: list[float], large_freqs: list[float]) -> None:
    small_positive = np.array([x for x in small_freqs if x > 0.1], dtype=float)
    large_positive = np.sort(np.array([x for x in large_freqs if x > 0.1], dtype=float))
    if small_positive.size == 0 or large_positive.size == 0:
        return

    diffs = []
    for freq in small_positive:
        index = int(np.searchsorted(large_positive, freq))
        options = []
        if index < large_positive.size:
            options.append(abs(large_positive[index] - freq))
        if index > 0:
            options.append(abs(large_positive[index - 1] - freq))
        diffs.append(min(options))

    diffs = np.sort(np.array(diffs))
    print("Frequency nearest-neighbour check, positive modes only:")
    print(f"  small_positive={small_positive.size}, large_positive={large_positive.size}")
    print(
        "  "
        f"mean_diff={np.mean(diffs):.6f} THz, "
        f"p50={np.quantile(diffs, 0.50):.6f}, "
        f"p90={np.quantile(diffs, 0.90):.6f}, "
        f"p99={np.quantile(diffs, 0.99):.6f}, "
        f"max={np.max(diffs):.6f}"
    )


def print_metadata(label: str, data: dict) -> None:
    lengths = [np.linalg.norm(row) for row in data["lattice"]]
    print(f"{label}: {data['path']}")
    print(f"  natom={data['natom']}, expected_modes_per_q={3 * data['natom']}, qpoints={len(data['q_positions'])}")
    print("  species=" + repr(dict(Counter(data["symbols"]))))
    print("  lattice_lengths_A=" + ", ".join(f"{x:.6f}" for x in lengths))
    print("  q_positions=" + repr(data["q_positions"]))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("small_band_yaml", type=Path)
    parser.add_argument("large_band_yaml", type=Path)
    parser.add_argument("--position-tolerance", type=float, default=0.15)
    args = parser.parse_args()

    small = parse_band_yaml(args.small_band_yaml)
    large = parse_band_yaml(args.large_band_yaml)

    print_metadata("SMALL", small)
    print_metadata("LARGE", large)
    print()

    summarize_frequencies("SMALL", small["freqs_by_q"])
    summarize_frequencies("LARGE", large["freqs_by_q"])
    print()

    compare_position_embedding(small, large, args.position_tolerance)
    print()

    nearest_frequency_report(small["freqs_by_q"][0], large["freqs_by_q"][0])


if __name__ == "__main__":
    main()
