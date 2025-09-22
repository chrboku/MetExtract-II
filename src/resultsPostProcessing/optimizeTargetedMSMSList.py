from .utils import Bunch
import csv


def getGreedyMSMSList():
    pass


if __name__ == "__main__":
    inF = "D:/180615_438_Labelboxpaper/results.tsv"
    samplesToUse = [
        "pos_12CKarurDON_13CKarurDON_Rep1.mzXML",
        "pos_12CKarurDON_13CP17x54DON_Rep2.mzXML",
        "pos_12CKarurMock_13CKarurMock_Rep2.mzXML",
    ]
    numCol = "Num"
    rtCol = "RT"  ## min
    mzCol = "MZ"

    maxSimMSMS = 5
    rtRange = 0.1  ## +/- min

    headers = {}
    data = []
    with open(inF, "r", encoding="utf-8") as fIn:
        reader = csv.reader(fIn, delimiter="\t")

        for rowi, row in enumerate(reader):
            if len(row) == 0:  # Skip empty lines
                continue
            if rowi == 0:
                for coli, header in enumerate(row):
                    headers[header] = coli

            elif not row[0].startswith("#"):
                a = Bunch(
                    num=row[headers[numCol]],
                    rt=float(row[headers[rtCol]]),
                    mz=float(row[headers[mzCol]]),
                )
                data.append(a)

            else:
                pass
