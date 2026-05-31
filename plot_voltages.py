#!/usr/bin/env python3
"""
Plot three-phase voltages from simulation CSV output.
"""

import pandas as pd
import matplotlib.pyplot as plt
import sys


def list_columns(csv_file):
    """Print available CSV columns with their index."""
    df = pd.read_csv(csv_file, nrows=1)
    print("Available columns:")
    for idx, col in enumerate(df.columns):
        print(f"  [{idx}] {col}")


def plot_voltages(csv_file, time_column=None, value_columns=None):
    """
    Read CSV file and plot selected columns over time.
    
    Args:
        csv_file: Path to CSV file.
        time_column: Name of the column to use as x-axis.
        value_columns: List of column names to plot on y-axis.
    """
    df = pd.read_csv(csv_file)

    time_col = time_column if time_column else df.columns[0]
    if time_col not in df.columns:
        raise ValueError(f"Time column '{time_col}' not found in CSV.")

    if value_columns:
        missing = [col for col in value_columns if col not in df.columns]
        if missing:
            raise ValueError(f"Value columns not found in CSV: {missing}")
        y_cols = value_columns
    else:
        # Default for a 3-column CSV: first is time, second/third are signals (d and q).
        if len(df.columns) >= 3:
            y_cols = [df.columns[1], df.columns[2]]
        else:
            y_cols = [col for col in df.columns if col != time_col][:2]

    if not y_cols:
        raise ValueError("No value columns selected for plotting.")

    plt.figure(figsize=(12, 6))

    for col in y_cols:
        plt.plot(df[time_col], df[col], label=str(col), linewidth=1.5)
    
    plt.xlabel('Time (s)', fontsize=12)
    plt.ylabel('Voltage (V)', fontsize=12)
    plt.title('d-q Voltages Over Time', fontsize=14, fontweight='bold')
    plt.legend(loc='best')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    
    # Save the plot
    output_file = csv_file.replace('.csv', '_plot.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Plot saved to: {output_file}")
    print(f"Using time column: {time_col}")
    print(f"Plotted columns: {y_cols}")
    
    # Also try to show the plot (will only work if display is available)
    try:
        plt.show()
    except:
        pass

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_voltages.py <csv_file> [--list-columns] [--time-column <name>] [--columns <name1,name2,...>]")
        print("Example: python plot_voltages.py simulation_output.csv --time-column time --columns v_d,v_q")
        sys.exit(1)

    csv_file = sys.argv[1]
    args = sys.argv[2:]

    if '--list-columns' in args:
        list_columns(csv_file)
        sys.exit(0)

    time_column = None
    value_columns = None

    if '--time-column' in args:
        idx = args.index('--time-column')
        if idx + 1 >= len(args):
            print("Error: --time-column requires a column name")
            sys.exit(1)
        time_column = args[idx + 1]

    if '--columns' in args:
        idx = args.index('--columns')
        if idx + 1 >= len(args):
            print("Error: --columns requires a comma-separated list")
            sys.exit(1)
        value_columns = [c.strip() for c in args[idx + 1].split(',') if c.strip()]

    plot_voltages(csv_file, time_column=time_column, value_columns=value_columns)
