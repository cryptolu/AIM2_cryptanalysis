#!/usr/bin/env python3
import re
import glob
import os
from collections import defaultdict

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
import matplotlib.pyplot as plt

# Methods to extract from logs
METHODS = {
    "Our total time": "Our total time",
    "determinant_via_evaluation_geometric": "determinant_via_evaluation_geometric",
    "determinant_via_evaluation_geometric_sparse": "determinant_via_evaluation_geometric_sparse",
    "determinant_via_evaluation_FFT": "determinant_via_evaluation_FFT",
    "determinant_via_linsolve": "determinant_via_linsolve",
    "determinant_generic_knowing_degree": "determinant_generic_knowing_degree",
}

LOG_RE = re.compile(r"logs/main_W(\d+)_DEG(\d+)\.log$")
TIME_RE = re.compile(r"^\[T\]\s+([^:]+):\s+(\d+)\s+ms")


def parse_log(path):
    times = {}
    with open(path, "r", errors="ignore") as f:
        for line in f:
            m = TIME_RE.match(line.strip())
            if not m:
                continue
            label, ms = m.group(1), int(m.group(2))
            times[label] = ms / 1000.0
    return times


def main():
    logs = glob.glob("logs/main_W*_DEG*.log")
    if not logs:
        print("No logs found under ./logs")
        return

    data = defaultdict(lambda: defaultdict(dict))  # data[W][method][DEG] = seconds

    for path in logs:
        m = LOG_RE.search(path)
        if not m:
            continue
        W = int(m.group(1))
        DEG = int(m.group(2))
        times = parse_log(path)
        if not times:
            continue
        for method, label in METHODS.items():
            if label in times:
                data[W][method][DEG] = times[label]

    os.makedirs("plots", exist_ok=True)

    if not data:
        print("No timing lines found in logs. Ensure logs contain lines like '[T] ...: <ms> ms'.")
        return

    styles = [
        ("-", "o"),
        ("--", "s"),
        (":", "^"),
        ("-.", "D"),
        ("-", "x"),
        ("--", "+"),
    ]
    styles_map = {}
    for i, method in enumerate(sorted(METHODS)):
        styles_map[method] = styles[i]

    for W, methods in sorted(data.items()):
        plt.figure(figsize=(8, 5))
        for idx, (method, series) in enumerate(sorted(methods.items())):
            if not series:
                continue
            print("   ", W, ":", idx, method)
            xs = sorted(series.keys())
            ys = [series[x] for x in xs]
            ls, mk = styles_map[method]
            if method == "Our total time":
                method = "Our sparse method"
            else:
                method = "PML:" + method
            # ls, mk = styles[idx % len(styles)]
            plt.plot(xs, ys, linestyle=ls, marker=mk, label=method, markerfacecolor="none")

        plt.title(f"Benchmark timings for W={W}")
        plt.xlabel("Degree")
        plt.ylabel("Time (s)")
        plt.xscale("log", base=2)
        plt.yscale("log", base=2)
        plt.grid(True, which="both", ls="--", alpha=0.3)
        plt.legend()
        out = f"plots/W{W}.png"
        plt.tight_layout()
        plt.savefig(out, dpi=300)
        plt.close()
        print(f"Wrote {out}")


if __name__ == "__main__":
    main()
