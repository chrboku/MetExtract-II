from __future__ import print_function, division, absolute_import
import polars as pl
from .utils import add_sheet_to_excel


import pprint

pp = pprint.PrettyPrinter(indent=1)


# defines a group (input data) for which statistic columns (_Stat_) are calculated
def addGroup(to, groupName, minCount, cols):
    to[groupName] = {}

    to[groupName]["minCount"] = minCount
    to[groupName]["colsNames"] = cols
    to[groupName]["cols"] = []


# matches the input columns to the respective statistic column
def matchRows(groups, rowNames):
    for groupName, groupProps in groups.items():
        for rowName in groupProps["colsNames"]:
            for i in range(len(rowNames)):
                if rowNames[i] == rowName:
                    groupProps["cols"].append(i)


# re-arrange data and add statistic columns
def addStatsColumnToResults(metaFile, groups, outputOrder, sheet_name="Sheet1", new_sheet_name="Sheet2", commentStartingCharacter="#"):
    """
    Add statistics columns to results matrix.

    Args:
        metaFile: Input Excel file path
        groups: Dictionary of group definitions
        toFile: Output Excel file path
        outputOrder: List of group names in desired output order
        sheet_name: Name of the sheet to read from and write to
        commentStartingCharacter: Character that marks comment rows
    """
    # Read Excel file using polars
    df = pl.read_excel(metaFile, sheet_name=sheet_name)
    headers = df.columns

    # Match columns to groups
    matchRows(groups, headers)

    # Find MZ and RT column indices for sorting
    mzInd = headers.index("MZ") if "MZ" in headers else -1
    rtInd = headers.index("RT") if "RT" in headers else -1

    # Add group statistics columns
    for groupName in outputOrder:
        # Count non-empty values in group columns
        col_names = [headers[i] for i in groups[groupName]["cols"]]

        # Create a new column that counts non-null, non-empty values
        count_expr = pl.lit(0, pl.Int64)
        for col_name in col_names:
            # Check if column value is not null and not empty string
            count_expr = count_expr + pl.when((pl.col(col_name).is_not_null()) & (pl.col(col_name).cast(pl.Utf8).str.strip_chars() != "")).then(1).otherwise(0)

        df = df.with_columns(count_expr.alias(groupName))

    # Sort by RT and MZ if both columns exist
    if rtInd != -1 and mzInd != -1:
        df = df.sort(["RT", "MZ"])

    # Write the dataframe to Excel
    add_sheet_to_excel(metaFile, df, new_sheet_name)


# omit those columns that have less results than required by the user input
def performGroupOmit(infile, groupStats, sheet_name="Sheet1", new_sheet_name="FilteredResults"):
    """
    Filter rows based on group statistics and create separate sheets for omitted and false positive features.

    Args:
        infile: Input Excel file path
        groupStats: List of tuples (group_name, min_count, omit_flag, false_positive_flag)
        outfile: Output Excel file path
        sheet_name: Name of the sheet to read from
    """
    # Read Excel file using polars
    df = pl.read_excel(infile, sheet_name=sheet_name)
    headers = df.columns

    # Initialize filter expressions
    use_filter = pl.lit(False)
    all_gomit = True
    false_positive_filter = pl.lit(False)

    for gname, gmin, gomit, gremoveAsFalsePositive in groupStats:
        if gname in headers:
            if gomit:
                # Include rows where this group has >= min_count
                use_filter = use_filter | (pl.col(gname).cast(pl.Int64) >= gmin)
                all_gomit = False
            if gremoveAsFalsePositive:
                # Mark as false positive if any value > 0 in this group
                false_positive_filter = false_positive_filter | (pl.col(gname).cast(pl.Int64) > 0)

    # Create filtered dataframes
    if all_gomit:
        # No omit filters applied, keep all data
        data_df = df
        not_used_df = pl.DataFrame()
    else:
        data_df = df.filter(use_filter)
        not_used_df = df.filter(~use_filter)

    false_positives_df = df.filter(false_positive_filter)

    add_sheet_to_excel(infile, data_df, new_sheet_name)
    if len(not_used_df) > 0:
        add_sheet_to_excel(infile, not_used_df, new_sheet_name + "_Omitted")

    if len(false_positives_df) > 0:
        add_sheet_to_excel(infile, false_positives_df, new_sheet_name + "_FalsePositives")
