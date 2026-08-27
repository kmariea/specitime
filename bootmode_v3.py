#!/usr/bin/env python3
"""
BOOTMODE Python v3

Port of the Perl BOOTMODE Version 1.0 payload extracted from the original
Windows executable. It calculates the mean of bootstrapped half-range modes
(M-HRMB), optionally calculates a bootstrap standard error and percentile
confidence interval, and prints reference statistics for the original data.

v3 performance-oriented changes:
- Optimized recursive half-range mode calculation. The original port counted
  values in each candidate half-range with a nested scan, making mode
  estimation roughly O(n^2) per recursive level. This version uses sorted-list
  index windows plus bisect operations, preserving the original tie-breaking
  behavior while reducing the main scan to near O(n) per level.
- Avoids re-sorting inside recursive mode calls.
- Uses random.choices() for bootstrap resampling.
- Throttles terminal progress updates so large runs do not spend excessive
  time flushing stdout.
- When bootstrap standard error is requested, prompts separately for outer
  and inner bootstrap iteration counts. Press Enter at the inner prompt to use
  the outer count, matching the original behavior.

Input format: plain text file with one numeric value per line. Lines that do
not contain a parseable number are ignored.
"""

from __future__ import annotations

import bisect
import csv
import math
import os
import random
import re
import time
from dataclasses import dataclass
from pathlib import Path
from statistics import mean, median
from typing import Iterable, List, Optional, Sequence, Tuple


NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


@dataclass
class BootmodeResult:
    label: str
    n: int
    source_mean: float
    source_sem: Optional[float]
    source_median: float
    source_sd: Optional[float]
    source_half_range_mode: float
    boot_mean_of_means: Optional[float]
    mode: Optional[float]
    bootstrap_std_error: Optional[float]
    ci_low: Optional[float]
    ci_high: Optional[float]
    modes: List[float]
    outer_iterations: int
    inner_iterations: Optional[int]


class ProgressPrinter:
    """Single-line progress renderer with time-based throttling."""

    def __init__(self, label: str, min_interval_seconds: float = 0.25) -> None:
        self.label = label
        self.min_interval_seconds = min_interval_seconds
        self.last_update = 0.0
        self.last_width = 0

    def update(self, completed: int, total: int, detail: str = "", force: bool = False) -> None:
        if total <= 0:
            return
        now = time.monotonic()
        if not force and completed < total and (now - self.last_update) < self.min_interval_seconds:
            return
        pct = (completed / total) * 100.0
        message = f"{self.label}: {pct:6.2f}% complete ({completed}/{total})"
        if detail:
            message += f" {detail}"
        width = max(self.last_width, len(message), 80)
        print("\r" + message.ljust(width), end="", flush=True)
        self.last_update = now
        self.last_width = width

    def finish(self, completed: int, total: int) -> None:
        self.update(completed, total, force=True)
        print("", flush=True)


def sample_standard_deviation(values: Sequence[float]) -> Optional[float]:
    """Sample standard deviation, matching Perl Statistics::Descriptive behavior."""
    n = len(values)
    if n < 2:
        return None
    m = mean(values)
    return math.sqrt(sum((x - m) ** 2 for x in values) / (n - 1))


def _leftmost_indexes(sorted_values: Sequence[float]) -> List[int]:
    """For each position, return the index of the first equal value."""
    n = len(sorted_values)
    leftmost = [0] * n
    i = 0
    while i < n:
        j = i + 1
        while j < n and sorted_values[j] == sorted_values[i]:
            j += 1
        for k in range(i, j):
            leftmost[k] = i
        i = j
    return leftmost


def findmode(values: Iterable[float]) -> float:
    """Optimized recursive half-range mode estimator.

    This preserves the Perl payload's key interval conventions:
    - candidate bin counts use [times[j], times[j] + width], inclusive;
    - the selected recursive subset usually uses [times[j], times[j] + width),
      with an exclusive upper bound;
    - tied candidate bins are resolved by the narrowest inclusive occupied
      range, then by the union of tied ranges if still tied.
    """
    times = sorted(float(x) for x in values)
    return _findmode_sorted(times)


def _findmode_sorted(times: List[float]) -> float:
    # Iterative form avoids Python recursion overhead while applying exactly
    # the same successive narrowing logic.
    while True:
        num = len(times)
        if num == 0:
            raise ValueError("No datapoints")
        if num == 1:
            return times[0]
        if num == 2:
            return (times[0] + times[1]) / 2.0

        low = times[0]
        high = times[-1]
        width = (high - low) / 2.0

        leftmost = _leftmost_indexes(times)
        bin_sizes: List[int] = []
        right = 0

        # Match original loop: for j in range(num - 1). The upper edge is
        # inclusive for bin-size counting.
        for j in range(num - 1):
            upper = times[j] + width
            while right < num and times[right] <= upper:
                right += 1
            bin_sizes.append(right - leftmost[j])

        max_bin_size = max(bin_sizes)
        max_bin_indexes = [i for i, b in enumerate(bin_sizes) if b == max_bin_size]

        if len(max_bin_indexes) == 1:
            jtemp = max_bin_indexes[0]
            upper = times[jtemp] + width
            lo = bisect.bisect_left(times, times[jtemp])
            hi = bisect.bisect_left(times, upper)
            times = times[lo:hi]
        else:
            ranges: List[float] = []
            mins: List[float] = []
            maxes: List[float] = []

            for idx in max_bin_indexes:
                upper = times[idx] + width
                lo = bisect.bisect_left(times, times[idx])
                hi = bisect.bisect_right(times, upper)
                # hi should be greater than lo because times[idx] is included.
                ranges.append(times[hi - 1] - times[lo])
                mins.append(times[lo])
                maxes.append(times[hi - 1])

            width_new = min(ranges)
            min_range_indexes = [i for i, r in enumerate(ranges) if r == width_new]

            if len(min_range_indexes) == 1:
                j = min_range_indexes[0]
                jtemp = max_bin_indexes[j]
                upper = times[jtemp] + width
                lo = bisect.bisect_left(times, times[jtemp])
                hi = bisect.bisect_left(times, upper)
                times = times[lo:hi]
            else:
                vmin = min(mins)
                vmax = max(maxes)
                lo = bisect.bisect_left(times, vmin)
                hi = bisect.bisect_right(times, vmax)
                times = times[lo:hi]

        if len(times) == num:
            one = times[1] - times[0]
            two = times[-1] - times[-2]
            if one < two:
                times = times[1:]
            elif one > two:
                times = times[:-1]
            else:
                times = times[1:-1]


def bootstrap_sample(values: Sequence[float]) -> List[float]:
    return random.choices(values, k=len(values))


def percentile_index(rounds: int, percentile: float) -> int:
    """Replicate the original ceil(rounds * percentile / 100) index convention."""
    idx = math.ceil((rounds * percentile) / 100.0)
    if idx >= rounds:
        idx = rounds - 1
    if idx < 0:
        idx = 0
    return idx


def run_bootmode(
    data: Sequence[float],
    label: str,
    want_se: bool,
    outer_iterations: int,
    inner_iterations: Optional[int] = None,
    progress: bool = True,
) -> BootmodeResult:
    if not data:
        raise ValueError(f"{label}: no datapoints")

    original = sorted(float(x) for x in data)
    n = len(original)
    source_mean = mean(original)
    source_median = median(original)
    source_sd = sample_standard_deviation(original)
    source_sem = (source_sd / math.sqrt(n)) if source_sd is not None else None
    source_half_range_mode = _findmode_sorted(original[:])

    modes: List[float] = []
    means: List[float] = []

    if want_se:
        if inner_iterations is None:
            inner_iterations = outer_iterations
        total_steps = outer_iterations * inner_iterations
    else:
        total_steps = outer_iterations

    completed_steps = 0
    progress_printer = ProgressPrinter(label) if progress else None

    for outer in range(1, outer_iterations + 1):
        boot_data = bootstrap_sample(original)
        means.append(mean(boot_data))

        if want_se:
            sub_modes: List[float] = []
            assert inner_iterations is not None
            for inner in range(1, inner_iterations + 1):
                sub_data = bootstrap_sample(boot_data)
                sub_modes.append(findmode(sub_data))
                completed_steps += 1
                if progress_printer and completed_steps < total_steps:
                    progress_printer.update(
                        completed_steps,
                        total_steps,
                        detail=f"outer {outer}/{outer_iterations}, inner {inner}/{inner_iterations}",
                    )
            mode_i = mean(sub_modes)
        else:
            mode_i = findmode(boot_data)
            completed_steps += 1
            if progress_printer and completed_steps < total_steps:
                progress_printer.update(
                    completed_steps,
                    total_steps,
                    detail=f"iteration {outer}/{outer_iterations}",
                )

        modes.append(mode_i)

    if progress_printer:
        progress_printer.finish(completed_steps, total_steps)

    boot_mean_of_means: Optional[float] = None
    mode: Optional[float] = None
    bootstrap_std_error: Optional[float] = None
    ci_low: Optional[float] = None
    ci_high: Optional[float] = None

    if outer_iterations > 0 and modes:
        boot_mean_of_means = mean(means)
        mode = mean(modes)

        if want_se and len(modes) > 1:
            bootstrap_std_error = math.sqrt(sum((m - mode) ** 2 for m in modes) / (len(modes) - 1))
            sorted_modes = sorted(modes)
            ci_low = sorted_modes[percentile_index(len(sorted_modes), 2.5)]
            ci_high = sorted_modes[percentile_index(len(sorted_modes), 97.5)]

    return BootmodeResult(
        label=label,
        n=n,
        source_mean=source_mean,
        source_sem=source_sem,
        source_median=source_median,
        source_sd=source_sd,
        source_half_range_mode=source_half_range_mode,
        boot_mean_of_means=boot_mean_of_means,
        mode=mode,
        bootstrap_std_error=bootstrap_std_error,
        ci_low=ci_low,
        ci_high=ci_high,
        modes=modes,
        outer_iterations=outer_iterations,
        inner_iterations=inner_iterations if want_se else None,
    )


def read_numeric_data(path: Path) -> List[float]:
    values: List[float] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = NUMBER_RE.search(line)
            if match:
                values.append(float(match.group(0)))
    return values


def make_value_width_bins(values: Sequence[float], bin_size: float) -> List[Tuple[str, Optional[float], Optional[float], List[float]]]:
    """Partition values into zero-anchored contiguous bins of numeric width bin_size.

    For non-negative values, the first bin is [0, bin_size), the next is
    [bin_size, 2 * bin_size), and so on. Negative values, if present, are placed
    into matching bins below zero, such as [-bin_size, 0).
    """
    if bin_size <= 0:
        raise ValueError("Bin size must be greater than zero")
    if not values:
        return []

    grouped: dict[int, List[float]] = {}
    for value in values:
        bin_index = math.floor(value / bin_size)
        grouped.setdefault(bin_index, []).append(float(value))

    bins: List[Tuple[str, Optional[float], Optional[float], List[float]]] = []
    for index in sorted(grouped):
        start = index * bin_size
        end = start + bin_size
        label = f"Bin {index + 1}: [{start:g}, {end:g})"
        bins.append((label, start, end, sorted(grouped[index])))

    return bins


def format_optional(value: Optional[float]) -> str:
    return "N/A" if value is None else f"{value:g}"


def render_result(result: BootmodeResult, want_se: bool) -> str:
    lines = [
        f"\n===== {result.label} =====",
        f"Datapoints = {result.n}",
        f"Outer bootstrap iterations = {result.outer_iterations}",
    ]
    if want_se:
        lines.append(f"Inner bootstrap iterations per outer sample = {result.inner_iterations}")

    lines.extend(
        [
            "\nRESULTS:",
            f" Mode = {format_optional(result.mode)}",
        ]
    )
    if want_se:
        lines.extend(
            [
                f"Bootstrap Std. Error of Modes = {format_optional(result.bootstrap_std_error)}",
                f"95% Confidence Interval (range) = {format_optional(result.ci_low)} - {format_optional(result.ci_high)}",
            ]
        )
    lines.extend(
        [
            "\nSTATISTICS OF ORIGINAL DATA (for reference):",
            f" Mean = {result.source_mean:g}",
            f" Standard Error of the Mean = {format_optional(result.source_sem)}",
            f" Median = {result.source_median:g}",
            f" Half-Range Mode = {result.source_half_range_mode:g}",
        ]
    )
    return "\n".join(lines) + "\n"


def prompt_csv_output_path() -> Path:
    while True:
        raw = input("CSV results filename: ").strip().strip('"')
        if not raw:
            print("Please enter a CSV filename.")
            continue

        path = Path(raw)
        if path.suffix == "":
            path = path.with_suffix(".csv")

        parent = path.parent if str(path.parent) else Path(".")
        if not parent.exists():
            print(f"Directory does not exist: {parent}")
            continue
        if path.exists() and path.is_dir():
            print("That path is a directory. Please enter a file path.")
            continue

        try:
            with path.open("a", encoding="utf-8", newline=""):
                pass
            if not os.access(path, os.W_OK):
                raise PermissionError(f"Not writable: {path}")
        except OSError as exc:
            print(f"Cannot write to {path}: {exc}")
            continue

        return path


def csv_value(value: Optional[float]) -> str:
    return "" if value is None else f"{value:.12g}"


def result_to_csv_row(
    result: BootmodeResult,
    source_file: str,
    bin_start: Optional[float],
    bin_end: Optional[float],
) -> dict[str, str | int]:
    return {
        "source_file": source_file,
        "bin_label": result.label,
        "bin_start_inclusive": csv_value(bin_start),
        "bin_end_exclusive": csv_value(bin_end),
        "n": result.n,
        "outer_iterations": result.outer_iterations,
        "inner_iterations": "" if result.inner_iterations is None else result.inner_iterations,
        "source_mean": csv_value(result.source_mean),
        "source_sem": csv_value(result.source_sem),
        "source_median": csv_value(result.source_median),
        "source_sd": csv_value(result.source_sd),
        "source_half_range_mode": csv_value(result.source_half_range_mode),
        "bootstrap_mean_of_means": csv_value(result.boot_mean_of_means),
        "mode": csv_value(result.mode),
        "bootstrap_std_error": csv_value(result.bootstrap_std_error),
        "ci_low_2_5_percent": csv_value(result.ci_low),
        "ci_high_97_5_percent": csv_value(result.ci_high),
    }


def write_csv_results(path: Path, rows: Sequence[dict[str, str | int]]) -> None:
    fieldnames = [
        "source_file",
        "bin_label",
        "bin_start_inclusive",
        "bin_end_exclusive",
        "n",
        "outer_iterations",
        "inner_iterations",
        "source_mean",
        "source_sem",
        "source_median",
        "source_sd",
        "source_half_range_mode",
        "bootstrap_mean_of_means",
        "mode",
        "bootstrap_std_error",
        "ci_low_2_5_percent",
        "ci_high_97_5_percent",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


EXPLANATION = """
EXPLANATION:
The "Mode" is the mean of the bootstrapped half-range modes (M-HRMB). The
"Standard Error of the Mode" is a bootstrap standard error, calculated by
taking the standard deviation of the bootstrapped mode estimates. During the
error calculation, higher-level bootstrap iterations resample the original
data, and lower-level bootstrap iterations treat each of those higher-level
bootstrapped data sets as a starting point for lower-level resampling to
calculate the individual M-HRMB estimates. The "95% confidence interval" range
accounts for both symmetric and asymmetric distribution, as it specifies the
range between the 2.5% and 97.5% bootstrap values generated.

For application of this program, and for citation, refer to:
Hedges, S. B. and P. Shah. 2003. Comparison of mode estimation methods and
application in molecular clock analysis. BMC Bioinformatics 4:xx.
"""


def prompt_yes_no(prompt: str) -> bool:
    while True:
        answer = input(prompt).strip().lower()
        if answer in {"y", "yes"}:
            return True
        if answer in {"n", "no"}:
            return False
        print("Please answer y or n.")


def prompt_positive_int(prompt: str) -> int:
    while True:
        raw = input(prompt).strip()
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a positive integer.")
            continue
        if value > 0:
            return value
        print("Please enter a positive integer.")


def prompt_positive_int_default(prompt: str, default: int) -> int:
    while True:
        raw = input(prompt).strip()
        if raw == "":
            return default
        try:
            value = int(raw)
        except ValueError:
            print("Please enter a positive integer or press Enter for the default.")
            continue
        if value > 0:
            return value
        print("Please enter a positive integer or press Enter for the default.")


def prompt_optional_bin_size() -> Optional[float]:
    while True:
        raw = input(
            "Optional bin size by numeric value width; leave blank to analyze the entire dataset: "
        ).strip()
        if raw == "":
            return None
        try:
            value = float(raw)
        except ValueError:
            print("Please enter a numeric bin size or leave blank.")
            continue
        if value > 0:
            return value
        print("Bin size must be greater than zero.")


def main() -> None:
    print(
        "\nBOOTMODE Python Version 3\n"
        "This program calculates the mode, standard error of the mode, 95% confidence\n"
        "interval of the mode, and other reference statistics. Input = text file\n"
        "with one number per row.\n"
    )

    datafile = input("Data file name : ").strip().strip('"')
    data_path = Path(datafile)
    if not data_path.exists():
        raise SystemExit(f"Cannot open datafile: {datafile}")

    csv_output_path = prompt_csv_output_path()

    data = read_numeric_data(data_path)
    if not data:
        raise SystemExit("No numeric datapoints were found in the data file.")

    bin_size = prompt_optional_bin_size()

    want_se = prompt_yes_no("Calculate Bootstrap Standard Error? (y/n) ")
    if want_se:
        outer_iterations = prompt_positive_int("Outer bootstrap iterations? (100 or more is recommended) ")
        inner_iterations = prompt_positive_int_default(
            f"Inner bootstrap iterations per outer sample? [default: {outer_iterations}] ",
            outer_iterations,
        )
    else:
        outer_iterations = prompt_positive_int("Bootstrap iterations? (1000 or more is recommended) ")
        inner_iterations = None

    if bin_size is None:
        datasets = [("Entire dataset", None, None, sorted(data))]
    else:
        datasets = make_value_width_bins(data, bin_size)
        if not datasets:
            raise SystemExit("No non-empty bins were created.")

    all_output: List[str] = [f"------------------Data File: {datafile}----------------\n"]
    all_modes_for_file: List[Tuple[str, List[float]]] = []
    csv_rows: List[dict[str, str | int]] = []

    for label, bin_start, bin_end, values in datasets:
        print(f"\nProcessing {label} ({len(values)} datapoints)")
        result = run_bootmode(
            values,
            label=label,
            want_se=want_se,
            outer_iterations=outer_iterations,
            inner_iterations=inner_iterations,
        )
        rendered = render_result(result, want_se)
        print(rendered)
        all_output.append(rendered)
        all_modes_for_file.append((label, result.modes))
        csv_rows.append(result_to_csv_row(result, datafile, bin_start, bin_end))

    write_csv_results(csv_output_path, csv_rows)

    print(EXPLANATION)
    all_output.append(EXPLANATION)

    with Path("Bootmode.out").open("a", encoding="utf-8") as out:
        out.write("\n".join(all_output))
        out.write("\n")

    if want_se:
        with Path("modes.txt").open("w", encoding="utf-8") as out:
            for label, modes in all_modes_for_file:
                if len(all_modes_for_file) > 1:
                    out.write(f"# {label}\n")
                for mode_value in sorted(modes):
                    out.write(f"{mode_value:g}\n")

    print(f"CSV summary written to {csv_output_path}")
    print("Detailed text results appended to Bootmode.out")
    if want_se:
        print("Bootstrap mode values written to modes.txt")


if __name__ == "__main__":
    main()
