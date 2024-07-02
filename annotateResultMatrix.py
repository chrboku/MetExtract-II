import csv

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
    for groupName, groupProps in groups.iteritems():
        for rowName in groupProps["colsNames"]:
            for i in range(len(rowNames)):
                if rowNames[i] == rowName:
                    groupProps["cols"].append(i)


# re-arrange data and add statistic columns
def addStatsColumnToResults(metaFile, groups, toFile, outputOrder, commentStartingCharacter="#"):
    data = []
    comments = []

    mzInd=-1
    rtInd=-1

    with open(metaFile, "rb") as x:
        meta = csv.reader(x, delimiter="\t")

        rowNum = 0
        for line in meta:
            if line[0].startswith(commentStartingCharacter):
                comments.append(line[0].strip())
                continue

            if rowNum == 0:
                headers = line
                j=0
                for header in headers:
                    if header=="MZ":
                        mzInd=j
                    if header=="RT":
                        rtInd=j
                    j+=1
                matchRows(groups, line)
            else:
                data.append(line)
            rowNum += 1

    for groupName in outputOrder:
        headers.append(groupName)

    for row in data:
        for groupName in outputOrder:
            found = 0
            for pos in groups[groupName]["cols"]:
                try:
                    if row[pos].strip() != "":
                        found += 1
                except:
                    pass;
            row.append(str(found))

    data = sorted(data, key=lambda x: (float(x[rtInd]), float(x[mzInd])))   # sort results according to rt and mz
    data.insert(0, headers)

    with open(toFile, "w") as x:
        metaWriter = csv.writer(x, delimiter="\t")
        for row in data:
            metaWriter.writerow(row)
        for comment in comments:
            metaWriter.writerow([comment])


# omit those columns that have less results than required by the user input
def performGroupOmit(infile, groupStats, outfile, commentStartingCharacter="#"):
    data = []
    notUsed = []
    falsePositives = []
    headers = {}
    hrow = []
    comments = []

    with open(infile, "rb") as x:
        meta = csv.reader(x, delimiter="\t")

        rowNum = 0
        for line in meta:
            if line[0].startswith(commentStartingCharacter):
                comments.append(line[0].strip())
                continue

            if rowNum == 0:
                hrow = line
                for i in range(len(line)):
                    headers[line[i]] = i

            else:
                use = False
                allGomit = True
                isFalsePositive = False
                for gname, gmin, gomit, gremoveAsFalsePositive in groupStats:
                    if gomit:
                        use = use or (int(line[headers[gname]]) >= gmin)
                        allGomit=False
                    if gremoveAsFalsePositive:
                        isFalsePositive = isFalsePositive or int(line[headers[gname]]) > 0

                if isFalsePositive:
                    falsePositives.append(line)
                if use or allGomit:     ## either use it because it was found more than n times in a group or because no omit was used
                    data.append(line)
                else:
                    notUsed.append(line)
            rowNum += 1

    with open(outfile, "w") as x:
        metaWriter = csv.writer(x, delimiter="\t")
        metaWriter.writerow(hrow)
        for row in data:
            metaWriter.writerow(row)
        for comment in comments:
            metaWriter.writerow([comment])
        metaWriter.writerow(["## most likely true positives"])


    if len(notUsed)>0:
        with open(outfile.replace(".tsv", ".omitteds.tsv"), "w") as x:
            metaWriter = csv.writer(x, delimiter="\t")
            metaWriter.writerow(hrow)
            for row in notUsed:
                metaWriter.writerow(row)
            for comment in comments:
                metaWriter.writerow([comment])
            metaWriter.writerow(["## features that have not been detected in a sufficiently high number of samples (group parameters omit)"])


    if len(notUsed)>0:
        with open(outfile.replace(".tsv", ".falsePositives.tsv"), "w") as x:
            metaWriter = csv.writer(x, delimiter="\t")
            metaWriter.writerow(hrow)
            for row in falsePositives:
                metaWriter.writerow(row)
            for comment in comments:
                metaWriter.writerow([comment])
            metaWriter.writerow(["## features that have  been detected in at least one 'blank' sample"])

