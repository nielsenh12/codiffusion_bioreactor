#!/usr/bin/env python3
"""
Linear interpolation script for filling gaps in Summary.csv

This script reads Summary.csv, identifies gaps (consecutive missing values)
in each numeric column, and fills them using linear interpolation between
the endpoints of each gap.
"""

import pandas as pd
import numpy as np
from pathlib import Path


def interpolate_gaps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Linearly interpolate gaps in numeric columns.

    For each numeric column, finds gaps (sequences of NaN values) and
    fills them with linearly interpolated values between the last valid
    value before the gap and the first valid value after the gap.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe with potential gaps

    Returns
    -------
    pd.DataFrame
        Dataframe with gaps filled via linear interpolation
    """
    df_interpolated = df.copy()

    for col in df.columns:
        # Skip non-numeric columns
        if not pd.api.types.is_numeric_dtype(df[col]):
            continue

        # Get the column data
        series = df_interpolated[col]

        # Skip if no missing values
        if not series.isna().any():
            continue

        # Skip if all values are missing (can't interpolate)
        if series.isna().all():
            continue

        # Use pandas interpolate with linear method
        # limit_direction='both' ensures we interpolate in gaps between valid values
        # but won't extrapolate beyond the first/last valid values
        df_interpolated[col] = series.interpolate(
            method='linear',
            limit_direction='forward',
            limit_area='inside'  # Only fill NaN values inside valid values (not at edges)
        )

    return df_interpolated


def main():
    # Define paths
    script_dir = Path(__file__).parent
    input_file = script_dir / "Summary.csv"
    output_file = script_dir / "Summary_interpolated.csv"

    print(f"Reading {input_file}...")

    # Read CSV - the file has a multi-level header
    df = pd.read_csv(input_file, header=[0, 1])

    print(f"Original shape: {df.shape}")

    # Count missing values before interpolation
    missing_before = df.isna().sum().sum()
    print(f"Total missing values before interpolation: {missing_before}")

    # Show missing values per column (only columns with gaps)
    print("\nColumns with gaps (before interpolation):")
    missing_per_col = df.isna().sum()
    for col, count in missing_per_col.items():
        if count > 0:
            print(f"  {col[1][:50]:50s}: {count:5d} missing")

    # Perform interpolation
    print("\nInterpolating gaps...")
    df_interpolated = interpolate_gaps(df)

    # Count missing values after interpolation
    missing_after = df_interpolated.isna().sum().sum()
    print(f"\nTotal missing values after interpolation: {missing_after}")
    print(f"Values filled: {missing_before - missing_after}")

    # Show remaining gaps (at edges where interpolation isn't possible)
    print("\nRemaining gaps (at edges, cannot interpolate):")
    remaining_per_col = df_interpolated.isna().sum()
    for col, count in remaining_per_col.items():
        if count > 0:
            print(f"  {col[1][:50]:50s}: {count:5d} remaining")

    # Save interpolated data
    print(f"\nSaving to {output_file}...")
    df_interpolated.to_csv(output_file, index=False)
    print("Done!")

    return df_interpolated


if __name__ == "__main__":
    main()
