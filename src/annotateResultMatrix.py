from __future__ import print_function, division, absolute_import
import polars as pl
import logging
from .utils import add_sheet_to_excel, getDBSuffix, getDBFormat
from .PolarsDB import PolarsDB
from .resultsPostProcessing import searchDatabases
from .resultsPostProcessing import generateSumFormulas as sumFormulaGeneration
from .formulaTools import formulaTools

import json

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
        sheet_name: Name of the sheet to read from (e.g., "4_Reintegrated")
        new_sheet_name: Name of the sheet to write to (e.g., "6_Annotated")
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
            (f"DBs_{dbName}_count", ""),
            (f"DBs_{dbName}", ""),
            (f"DBs_RT_{dbName}_count", ""),
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

                hit_str = f"(Name: {hit.name}, Type: {hit.hitType}, Num: {hit.num}, Formula: {hit.sumFormula}, RT: {hit.rt_min}, MassErrorPPM: {hit.matchErrorPPM:.5f}, MassErrorMass: {hit.matchErrorMass:.5f}, Additional information: {hit.additionalInfo})"
                hits_per_db[hit.dbName]["hits"].append(hit_str)

                if hit.rt_min != None and hit.rt_min != "":
                    try:
                        rtDelta = abs(float(hit.rt_min) - rt_min)
                        rt_hit_str = f"RT delta: {rtDelta:.2f} (Name: {hit.name} Type: {hit.hitType}, Num: {hit.num}, Formula: {hit.sumFormula}, RT: {hit.rt_min}, MassErrorPPM: {hit.matchErrorPPM:.5f}, MassErrorMass: {hit.matchErrorMass:.5f}, RTError: {float(hit.rt_min) - rt_min if hit.rt_min != None and hit.rt_min != '' else ''}, Additional information: {hit.additionalInfo})"
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
                    "Feature_OGroup": row.get("OGroup"),
                    "Feature_Ion": row.get("Ion"),
                    "Feature_Loss": row.get("Loss"),
                    "Feature_Relative_peakarea_in_group": row.get("Relative_peakarea_in_group"),
                    "Feature_Average_peakarea": row.get("Average_peakarea"),
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
                results_df[row_idx, f"DBs_{dbName}_count"] = len(hit_data["hits"]) if hit_data["hits"] else ""
                # results_df[row_idx, f"DBs_{dbName}"] = '"%s"' % (";\n".join(hit_data["hits"])) if hit_data["hits"] else ""
                results_df[row_idx, f"DBs_{dbName}"] = json.dumps(hit_data["hits"], ensure_ascii=True) if hit_data["hits"] else ""

                # Update RT hits if available
                if len(hit_data["hitsRT"]) > 0:
                    sorted_hits = sorted(hit_data["hitsRT"], key=lambda x: x[0])
                    results_df[row_idx, f"DBs_RT_{dbName}_count"] = len(hit_data["hitsRT"]) if hit_data["hitsRT"] else ""
                    # results_df[row_idx, f"DBs_RT_{dbName}"] = '"%s"' % (";\n".join([b[1] for b in sorted_hits])) if hit_data["hitsRT"] else ""
                    results_df[row_idx, f"DBs_RT_{dbName}"] = json.dumps([b[1] for b in sorted_hits], ensure_ascii=True) if hit_data["hitsRT"] else ""

        # Collect hit objects for compound sheet
        all_compound_hits.extend(hit_objects)

    # Save annotated results to new sheet
    logging.info(f"Saving annotated results to sheet: {new_sheet_name}")
    plDB.set_table(new_sheet_name, results_df)

    # Create compound-focused sheet
    if all_compound_hits:
        logging.info(f"Creating compound-focused sheet with {len(all_compound_hits)} database hits")

        # Collect all unique additionalInfo keys across all hits so every row gets a column
        all_additional_keys = []
        seen_keys = set()
        for hit, _row_info in all_compound_hits:
            for k in hit.additionalInfo.keys():
                if k not in seen_keys:
                    all_additional_keys.append(k)
                    seen_keys.add(k)

        compound_rows = []
        for hit, row_info in all_compound_hits:
            compound_row = {
                # Database entry information
                "DB_Name": str(hit.dbName) if hit.dbName is not None else "",
                "DB_Num": str(hit.num) if hit.num is not None else "",
                "DB_CompoundName": str(hit.name) if hit.name is not None else "",
                "DB_SumFormula": str(hit.sumFormula) if hit.sumFormula is not None else "",
                "DB_Mass": float(hit.mass) if hit.mass is not None and hit.mass != "" else None,
                "DB_RT_min": float(hit.rt_min) if hit.rt_min is not None and hit.rt_min != "" else None,
                "DB_MZ": float(hit.mz) if hit.mz is not None and hit.mz != "" else None,
                "DB_Polarity": str(hit.polarity) if hit.polarity is not None else "",
                "HitType": str(hit.hitType) if hit.hitType is not None else "",
                "MatchErrorPPM": float(hit.matchErrorPPM) if hit.matchErrorPPM is not None else None,
                "MatchErrorMass": float(hit.matchErrorMass) if hit.matchErrorMass is not None else None,
                # Feature information where the hit was found
                "Feature_Num": row_info["Feature_Num"],
                "Feature_OGroup": row_info["Feature_OGroup"],
                "Feature_RT": row_info["Feature_RT"],
                "Feature_MZ": row_info["Feature_MZ"],
                "Feature_Xn": row_info["Feature_Xn"],
                "Feature_Ionisation_Mode": row_info["Feature_Ionisation_Mode"],
                "Feature_Charge": row_info["Feature_Charge"],
                "Feature_M": row_info["Feature_M"],
                "Feature_Ion": row_info["Feature_Ion"],
                "Feature_Loss": row_info["Feature_Loss"],
                "Feature_Relative_peakarea_in_group": row_info["Feature_Relative_peakarea_in_group"],
                "Feature_Average_peakarea": row_info["Feature_Average_peakarea"],
            }
            # Expand additionalInfo into individual columns with a DB_Info_ prefix
            for k in all_additional_keys:
                compound_row[f"DB_Info_{k}"] = str(hit.additionalInfo.get(k, ""))
            compound_rows.append(compound_row)

        # Build schema overrides only for the known numeric columns; let polars infer the rest
        schema_overrides = {
            "DB_Mass": pl.Float64,
            "DB_RT_min": pl.Float64,
            "DB_MZ": pl.Float64,
            "MatchErrorPPM": pl.Float64,
            "MatchErrorMass": pl.Float64,
            "Feature_RT": pl.Float64,
            "Feature_MZ": pl.Float64,
            "Feature_M": pl.Float64,
            "Feature_Relative_peakarea_in_group": pl.Float64,
            "Feature_Average_peakarea": pl.Float64,
        }

        # Create dataframe for compound-focused sheet
        compound_df = pl.DataFrame(compound_rows, schema_overrides=schema_overrides, infer_schema_length=len(compound_rows))

        # Sort by database name, compound name, then feature id
        sort_cols = [c for c in ["DB_Name", "DB_CompoundName", "Feature_Num"] if c in compound_df.columns]
        if sort_cols:
            compound_df = compound_df.sort(sort_cols)

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


def annotateWithSumFormulas(
    file,
    sheet_name,
    useAtoms,
    atomsRange,
    processedElement,
    useExactXn,
    ppm,
    ppmCorrectionPosMode,
    ppmCorrectionNegMode,
    useAdducts,
    pwMaxSet=None,
    pwValSet=None,
    nCores=1,
):
    """
    Annotate metabolites with generated sum formulas using PolarsDB.

    This function:
    1. Loads the results table from the specified sheet using PolarsDB
    2. Generates sum formulas for each metabolite based on element ranges
    3. Adds sum formula columns with different element combinations
    4. Overwrites the input sheet with annotated results

    Args:
        file: Path to the results file (PolarsDB format)
        sheet_name: Name of the sheet to read from and write to (e.g., "6_Annotated")
        useAtoms: List of atoms to use in sum formulas (e.g., ["C", "H", "O", "N", "S"])
        atomsRange: List of [min, max] ranges for each atom
        processedElement: Element to check in formulas (e.g., "C")
        useExactXn: How to check Xn ("Exact", "Minimum", "Don't use", or "PlusMinus_X")
        ppm: Mass accuracy in ppm
        ppmCorrectionPosMode: PPM correction for positive mode
        ppmCorrectionNegMode: PPM correction for negative mode
        useAdducts: List of adduct definitions [[name, mzoffset, polarity, charge, mCount], ...]
        pwMaxSet: Progress callback for max value
        pwValSet: Progress callback for current value
        nCores: Number of CPU cores to use

    Returns:
        List of sum formula column names added
    """
    logging.info(f"Starting sum formula generation using PolarsDB")
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

    # Use the new Polars-based sum formula generation directly
    logging.info(f"Generating sum formulas with atoms: {useAtoms}")
    results_df, sf_compound_hits = sumFormulaGeneration.annotatePolarsTableWithSumFormulas(
        results_df,
        useAtoms,
        atomsRange,
        processedElement,
        useExactXn,
        ppm=ppm,
        ppmCorrectionPosMode=ppmCorrectionPosMode,
        ppmCorrectionNegMode=ppmCorrectionNegMode,
        adducts=useAdducts,
        pwMaxSet=pwMaxSet,
        pwValSet=pwValSet,
        nCores=nCores,
    )

    # Get list of added columns
    smCol = "SFs"
    annotationColumns = [
        f"{smCol}_CHO",
        f"{smCol}_CHOS",
        f"{smCol}_CHOP",
        f"{smCol}_CHON",
        f"{smCol}_CHONP",
        f"{smCol}_CHOPS",
        f"{smCol}_CHONS",
        f"{smCol}_CHONPS",
        f"{smCol}_all",
        f"{smCol}_CHO_count",
        f"{smCol}_CHOS_count",
        f"{smCol}_CHOP_count",
        f"{smCol}_CHON_count",
        f"{smCol}_CHONP_count",
        f"{smCol}_CHOPS_count",
        f"{smCol}_CHONS_count",
        f"{smCol}_CHONPS_count",
        f"{smCol}_all_count",
    ]

    # Overwrite the input sheet with annotated results
    logging.info(f"Overwriting sheet '{sheet_name}' with sum formula annotations")
    plDB.set_table(sheet_name, results_df)

    # Create formula-focused transposed sheet (one row per formula hit, mirroring the _Compounds sheet)
    sf_sheet_name = f"{sheet_name}_SumFormulas"
    if sf_compound_hits:
        logging.info(f"Creating sum formula compound sheet with {len(sf_compound_hits)} hits")
        schema_overrides = {
            "MassErrorPPM": pl.Float64,
            "MassErrorMass": pl.Float64,
            "Feature_RT": pl.Float64,
            "Feature_MZ": pl.Float64,
            # Feature_M can hold comma-separated neutral masses so keep it as a string
            "Feature_M": pl.Utf8,
            "Feature_Relative_peakarea_in_group": pl.Float64,
            "Feature_Average_peakarea": pl.Float64,
        }
        sf_hits_df = pl.DataFrame(
            sf_compound_hits,
            schema_overrides=schema_overrides,
            infer_schema_length=len(sf_compound_hits),
        )
        sort_cols = [c for c in ["Element_Class", "SumFormula", "Feature_Num"] if c in sf_hits_df.columns]
        if sort_cols:
            sf_hits_df = sf_hits_df.sort(sort_cols)
        logging.info(f"Saving sum formula compound sheet to: {sf_sheet_name}")
        plDB.set_table(sf_sheet_name, sf_hits_df)
    else:
        logging.info("No sum formula hits found, skipping compound-focused sheet")

    plDB.commit()
    plDB.close()

    logging.info(f"Sum formula generation completed. Added {len(annotationColumns)} columns")
    return annotationColumns
