from __future__ import print_function, division, absolute_import
import polars as pl
import logging
from .utils import add_sheet_to_excel, getDBSuffix, getDBFormat
from .PolarsDB import PolarsDB
from .resultsPostProcessing import searchDatabases
from .formulaTools import formulaTools


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


def annotateWithDatabases(
    file,
    sheet_name,
    new_sheet_name,
    dbFiles,
    useAdducts,
    ppm,
    correctppmPosMode,
    correctppmNegMode,
    rtError,
    useRt,
    checkXnInHits,
    processedElement,
    pwMaxSet=None,
    pwValSet=None,
):
    """
    Annotate metabolites by searching in databases using PolarsDB.

    This function:
    1. Loads the results table from the specified sheet using PolarsDB
    2. Imports database entries from provided files
    3. Adds all necessary database annotation columns
    4. Searches for database hits for each metabolite
    5. Updates only rows with hits using the [] operator
    6. Saves the annotated results to a new sheet

    Args:
        file: Path to the results file (PolarsDB format)
        sheet_name: Name of the sheet to read from (e.g., "5_afterReintegration")
        new_sheet_name: Name of the sheet to write to (e.g., "6_afterDBSearch")
        dbFiles: List of database file paths
        useAdducts: List of adduct definitions [[name, mzoffset, polarity, charge, mCount], ...]
        ppm: Mass accuracy in ppm
        correctppmPosMode: PPM correction for positive mode
        correctppmNegMode: PPM correction for negative mode
        rtError: Maximum RT error for matching
        useRt: Whether to use RT in matching
        checkXnInHits: How to check Xn ("Exact", "Minimum", "Don't use", or "PlusMinus_X")
        processedElement: Element to check in formulas (e.g., "C")
        pwMaxSet: Progress callback for max value
        pwValSet: Progress callback for current value

    Returns:
        List of annotation column names added
    """
    logging.info(f"Starting database search annotation using PolarsDB")
    logging.info(f"Reading file {file}")

    # Load the PolarsDB
    plDB = PolarsDB(file, format="xlsx", load_all_tables=True)

    # Load the results dataframe from the specified sheet
    try:
        results_df = plDB.get_table(sheet_name)
    except Exception as e:
        logging.error(f"Failed to load sheet '{sheet_name}': {e}")
        plDB.close()
        raise

    # Initialize database search
    db = searchDatabases.DBSearch()
    dbNames = []

    # Import database files and collect database names
    logging.info(f"Importing {len(dbFiles)} database file(s)")
    for dbFile in dbFiles:
        dbName = dbFile[dbFile.rfind("/") + 1 : dbFile.rfind(".")]
        try:
            imported, notImported = db.addEntriesFromFile(dbName, dbFile)
            if notImported > 0:
                logging.warning(f"Warning: {notImported} entries from database '{dbName}' were not imported successfully")
            dbNames.append(dbName)
            logging.info(f"  Imported {imported} entries from database '{dbName}'")
        except IOError as e:
            logging.error(f"Cannot open database file '{dbName}' at '{dbFile}': {e}")
            continue

    # Optimize database for searching
    logging.info("Optimizing database for searching")
    db.optimizeDB()

    # Add all necessary annotation columns BEFORE searching
    logging.info("Adding annotation columns to dataframe")
    annotationColumns = []
    for dbName in dbNames:
        cols_to_add = [
            (f"DBs_{dbName}_count", 0),
            (f"DBs_{dbName}", ""),
            (f"DBs_RT_{dbName}_count", 0),
            (f"DBs_RT_{dbName}", ""),
        ]

        for col_name, default_value in cols_to_add:
            annotationColumns.append(col_name)
            if col_name not in results_df.columns:
                if isinstance(default_value, int):
                    results_df = results_df.with_columns(pl.lit(default_value, dtype=pl.Int64).alias(col_name))
                else:
                    results_df = results_df.with_columns(pl.lit(default_value, dtype=pl.Utf8).alias(col_name))

    # Create database search helper class
    class DBSearchHelper:
        def __init__(
            self,
            dbs,
            ppm,
            correctppmPosMode,
            correctppmNegMode,
            rtError,
            useRt,
            checkXnInHits,
            processedElement,
        ):
            self.ppm = ppm
            self.correctppmPosMode = correctppmPosMode
            self.correctppmNegMode = correctppmNegMode
            self.rtError = rtError
            self.useRt = useRt
            self.checkXnInHits = checkXnInHits
            self.processedElement = processedElement
            self.dbs = dbs
            self.fT = formulaTools()

        def searchForRow(self, row):
            """
            Search databases for a single row and return dict of hits per database.

            Returns:
                tuple: (hits_per_db, hit_objects)
                    - hits_per_db: Dictionary with formatted hits strings
                      Format: {dbName: {'hits': [...], 'hitsRT': [...]}}
                    - hit_objects: List of (hit, row_info) tuples for compound sheet
            """
            # Extract values from row
            mass = row["M"] if row["M"] is not None and row["M"] != "" and "," not in str(row["M"]) else None
            if mass is not None:
                try:
                    mass = float(mass)
                except:
                    mass = None

            mz = float(row["MZ"])
            polarity = row["Ionisation_Mode"]
            rt_min = float(row["RT"])
            charges = int(row["Charge"])

            xn = 0
            try:
                xn = int(row["Xn"])
            except:
                pass

            # Apply mass correction
            if mass is not None:
                mass = mass + mass * 1.0 * (self.correctppmPosMode if polarity == "+" else self.correctppmNegMode) / 1e6
            mz = mz + mz * 1.0 * (self.correctppmPosMode if polarity == "+" else self.correctppmNegMode) / 1e6

            # Search database
            hits_per_db = {}
            hit_objects = []

            for hit in self.dbs.searchDB(
                mass=mass,
                mz=mz,
                polarity=polarity,
                charges=charges,
                rt_min=rt_min if self.useRt else None,
                ppm=self.ppm,
                rt_error=self.rtError,
                checkXN=self.checkXnInHits,
                element=self.processedElement,
                Xn=xn,
                adducts=useAdducts,
            ):
                if hit.dbName not in hits_per_db:
                    hits_per_db[hit.dbName] = {"hits": [], "hitsRT": []}

                hit_str = f"(Name: {hit.name}, Type: {hit.hitType}, Num: {hit.num}, Formula: {hit.sumFormula}, RT: {hit.rt_min}, MassErrorPPM: {hit.matchErrorPPM:.5f}, MassErrorMass: {hit.matchErrorMass:.5f}, RTError: {float(hit.rt_min) - rt_min if hit.rt_min != '' else ''}, Additional information: {hit.additionalInfo})"
                hits_per_db[hit.dbName]["hits"].append(hit_str)

                if hit.rt_min != "":
                    try:
                        rtDelta = abs(float(hit.rt_min) - rt_min)
                        rt_hit_str = f"RT delta: {rtDelta:.2f} (Name: {hit.name} Type: {hit.hitType}, Num: {hit.num}, Formula: {hit.sumFormula}, RT: {hit.rt_min}, MassErrorPPM: {hit.matchErrorPPM:.5f}, MassErrorMass: {hit.matchErrorMass:.5f}, RTError: {float(hit.rt_min) - rt_min if hit.rt_min != '' else ''}, Additional information: {hit.additionalInfo})"
                        hits_per_db[hit.dbName]["hitsRT"].append((rtDelta, rt_hit_str))
                    except Exception as e:
                        logging.error(f"Error processing RT for database hit: {e}")

                # Store hit object with row information for compound-focused sheet
                row_info = {
                    "Feature_Num": row.get("Num"),
                    "Feature_RT": row.get("RT"),
                    "Feature_MZ": row.get("MZ"),
                    "Feature_Xn": row.get("Xn"),
                    "Feature_Ionisation_Mode": row.get("Ionisation_Mode"),
                    "Feature_Charge": row.get("Charge"),
                    "Feature_M": row.get("M"),
                }
                hit_objects.append((hit, row_info))

            return hits_per_db, hit_objects

    # Create search helper instance
    searcher = DBSearchHelper(
        db,
        ppm,
        correctppmPosMode,
        correctppmNegMode,
        rtError,
        useRt,
        checkXnInHits,
        processedElement,
    )

    # Process each row and update only those with hits
    total_rows = len(results_df)
    if pwMaxSet is not None:
        pwMaxSet(total_rows)

    logging.info(f"Searching database hits for {total_rows} metabolites")

    # Collect all hits for compound-focused sheet
    all_compound_hits = []

    for row_idx in range(total_rows):
        if pwValSet is not None and row_idx % 10 == 0:
            pwValSet(row_idx)

        # Get row as dict
        row = results_df.row(row_idx, named=True)

        # Search for database hits
        hits_per_db, hit_objects = searcher.searchForRow(row)

        # Update only if there are hits
        if hits_per_db:
            for dbName, hit_data in hits_per_db.items():
                # Update count and hits columns
                results_df[row_idx, f"DBs_{dbName}_count"] = len(hit_data["hits"])
                results_df[row_idx, f"DBs_{dbName}"] = '"%s"' % (";".join(hit_data["hits"]))

                # Update RT hits if available
                if len(hit_data["hitsRT"]) > 0:
                    sorted_hits = sorted(hit_data["hitsRT"], key=lambda x: x[0])
                    results_df[row_idx, f"DBs_RT_{dbName}_count"] = len(hit_data["hitsRT"])
                    results_df[row_idx, f"DBs_RT_{dbName}"] = '"%s"' % (";".join([b[1] for b in sorted_hits]))

        # Collect hit objects for compound sheet
        all_compound_hits.extend(hit_objects)

    # Save annotated results to new sheet
    logging.info(f"Saving annotated results to sheet: {new_sheet_name}")
    plDB.set_table(new_sheet_name, results_df)

    # Create compound-focused sheet
    if all_compound_hits:
        logging.info(f"Creating compound-focused sheet with {len(all_compound_hits)} database hits")
        compound_rows = []

        for hit, row_info in all_compound_hits:
            compound_row = {
                # Database entry information
                "DB_Name": hit.dbName,
                "DB_Num": hit.num,
                "DB_CompoundName": hit.name,
                "DB_SumFormula": hit.sumFormula,
                "DB_Mass": hit.mass,
                "DB_RT_min": hit.rt_min if hit.rt_min else None,
                "DB_MZ": hit.mz if hit.mz else None,
                "DB_Polarity": hit.polarity,
                "DB_AdditionalInfo": str(hit.additionalInfo),
                "HitType": hit.hitType,
                "MatchErrorPPM": hit.matchErrorPPM,
                "MatchErrorMass": hit.matchErrorMass,
                # Feature information where the hit was found
                "Feature_Num": row_info["Feature_Num"],
                "Feature_RT": row_info["Feature_RT"],
                "Feature_MZ": row_info["Feature_MZ"],
                "Feature_Xn": row_info["Feature_Xn"],
                "Feature_Ionisation_Mode": row_info["Feature_Ionisation_Mode"],
                "Feature_Charge": row_info["Feature_Charge"],
                "Feature_M": row_info["Feature_M"],
            }
            compound_rows.append(compound_row)

        compound_datatypes = {
            "DB_Name": pl.Utf8,
            "DB_Num": pl.Utf8,
            "DB_CompoundName": pl.Utf8,
            "DB_SumFormula": pl.Utf8,
            "DB_Mass": pl.Float64,
            "DB_RT_min": pl.Float64,
            "DB_MZ": pl.Float64,
            "DB_Polarity": pl.Utf8,
            "DB_AdditionalInfo": pl.Utf8,
            "HitType": pl.Utf8,
            "MatchErrorPPM": pl.Float64,
            "MatchErrorMass": pl.Float64,
            "Feature_Num": pl.Int64,
            "Feature_RT": pl.Float64,
            "Feature_MZ": pl.Float64,
            "Feature_Xn": pl.Int64,
            "Feature_Ionisation_Mode": pl.Utf8,
            "Feature_Charge": pl.Int64,
            "Feature_M": pl.Float64,
        }

        # Create dataframe for compound-focused sheet
        compound_df = pl.DataFrame(compound_rows, schema=compound_datatypes)

        # Sort by database name and compound name
        compound_df = compound_df.sort(["DB_Name", "DB_CompoundName", "Feature_Num"])

        # Save to compound-focused sheet
        compound_sheet_name = f"{new_sheet_name}_Compounds"
        logging.info(f"Saving compound-focused results to sheet: {compound_sheet_name}")
        plDB.set_table(compound_sheet_name, compound_df)
    else:
        logging.info("No database hits found, skipping compound-focused sheet")

    plDB.commit()
    plDB.close()

    logging.info(f"Database search annotation completed. Added {len(annotationColumns)} columns")
    return annotationColumns
