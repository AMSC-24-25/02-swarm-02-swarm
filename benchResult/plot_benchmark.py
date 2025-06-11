import argparse
import os
import re
import pandas as pd
import matplotlib.pyplot as plt

def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate scaling and speedup plots from Google Benchmark results."
    )
    parser.add_argument(
        "--csv", "-c",
        default="results.csv",
        help="Input CSV file produced by Google Benchmark (default: results.csv)"
    )
    parser.add_argument(
        "--speedup-creatures", "-s",
        nargs="+",
        type=int,
        default=[],
        help="List of 'creatures' values for which to compute the strong speedup plot (e.g., 4 512 1024)."
    )
    parser.add_argument(
        "--output-dir", "-o",
        default=".",
        help="Destination folder for the generated plots (default: current directory)."
    )
    return parser.parse_args()

def main():
    args = parse_args()
    df = pd.read_csv(args.csv)

    # Filter only rows for the DeOpenMP_Rosenbrock benchmark
    df = df[df['name'].str.startswith('DeOpenMP_Rosenbrock')]

    # Extract creatures and threads from the name
    def extract_params(name):
        m = re.search(r'creatures:(\d+)/threads:(\d+)', name)
        return int(m.group(1)), int(m.group(2))
    df[['creatures', 'threads']] = df['name'].apply(lambda n: pd.Series(extract_params(n)))

    # Convert real_time from nanoseconds to seconds
    df['time_s'] = df['real_time'] * 1e-9

    os.makedirs(args.output_dir, exist_ok=True)

    # Plot 1: Time vs number of creatures for each number of threads
    plt.figure()
    for t in sorted(df['threads'].unique()):
        sub = df[df['threads'] == t]
        plt.plot(sub['creatures'], sub['time_s'], marker='o', label=f'{t} threads')
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.xlabel('Number of creatures')
    plt.ylabel('Time [s]')
    plt.title('Scaling: time vs number of creatures')
    plt.legend()
    plt.grid(True, which="both", ls="--", linewidth=0.5)
    plt.tight_layout()
    path1 = os.path.join(args.output_dir, 'time_vs_creatures.png')
    plt.savefig(path1, dpi=300)
    print(f"Plot saved: {path1}")

    # Plot 2: Speedup vs number of threads for selected creatures values
    if args.speedup_creatures:
        base = df[df['threads'] == 1].set_index('creatures')['time_s']
        plt.figure()
        for c in args.speedup_creatures:
            sub = df[df['creatures'] == c]
            if c not in base.index:
                print(f"WARNING: creatures value {c} not found for thread=1, skipping.")
                continue
            speedup = base[c] / sub['time_s'].values
            plt.plot(sub['threads'], speedup, marker='o', label=f'{c} creatures')
        plt.xscale('log', base=2)
        plt.xlabel('Number of threads')
        plt.ylabel('Speedup (t1 / tn)')
        plt.title('Strong speedup vs number of threads')
        plt.legend()
        plt.grid(True, which="both", ls="--", linewidth=0.5)
        plt.tight_layout()
        path2 = os.path.join(args.output_dir, 'speedup_vs_threads.png')
        plt.savefig(path2, dpi=300)
        print(f"Plot saved: {path2}")
    else:
        print("No speedup values specified, speedup plot not generated.")

if __name__ == "__main__":
    main()
