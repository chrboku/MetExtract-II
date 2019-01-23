from formulaTools import formulaTools


# class and functionality for calculating sum formulas for a given molecular mass
class SGRGenerator:
    def __init__(self, atoms=None):
        if atoms is None:
            self.atoms = {"C": [12.0, 4], "H": [1.007825, 1], "O": [15.994915, 6], "N": [14.003074, 5],
                          "P": [30.973763, 5], "S": [31.972072, 6], "Cl": [34.968853, 7], "Na": [22.98977],
                          "Fe": [55.934939, 2], "Cr": [51.94051, 6]};
        else:
            self.atoms = atoms

    #   Seven Golden Rules for chemical formula generation from molecular mass
    #   Kind and Fiehn (2006)
    #   DOI: 10.1186/1471-2105-8-105
    #
    #   Rule #1 (restriction of elements)                        Implemented
    #   Rule #2 (LEWIS and SENIOR chemical rules)                Implemented
    #   Rule #3 (isotope patterns)                               Not implemented; Not required; exact number of carbon atoms is available
    #   Rule #4 (hydrogen/carbon ratios)                         Implemented
    #   Rule #5 (element ratio of N, O, P, and S versus carbon)  Implemented
    #   Rule #6 (element ratio probabilities)                    Implemented
    #   Rule #7 (presence of trimethylsilylated compounds)       Not implemented; Not required when working with molecular mass
    #

    def secondGoldenRule(self, useAtoms, curAtoms):
        valences = 0
        atomsWithOddValences = 0
        maxVal = 0
        for i in range(0, len(useAtoms)):
            if curAtoms[i] != 0:
                valAtom = self.atoms[useAtoms[i]][1]

                valences = valences + valAtom * curAtoms[i]

                if bool(valAtom % 2):
                    atomsWithOddValences = atomsWithOddValences + curAtoms[i]

                maxVal = max(maxVal, valAtom)
        #print valences, atomsWithOddValences, maxVal, useAtoms, curAtoms, (not(bool(valences%2)) or not(bool(atomsWithOddValences%2))) and (valences>=(2*maxVal)) and (valences>=(2*(sum(curAtoms)-1)))
        return (not (bool(valences % 2)) or not (bool(atomsWithOddValences % 2))) and (valences >= (2 * maxVal)) and (
            valences >= (2 * (sum(curAtoms) - 1)))

    #Golden rule number 4 and 5
    def goldenRuleGeneric(self, useAtoms, atomsRange, i, curAtoms, fixed):
        lookup = {"H": [0.2, 3.1], "F": [0, 1.5], "Cl": [0, 0.8], "Br": [0, 0.8], "N": [0, 1.3], "O": [0, 1.2],
                  "P": [0, 0.3], "S": [0, 0.8], "Si": [0, 0.5]}
        if "C" in useAtoms:
            cCount = -1
            for j in range(0, len(useAtoms)):
                if useAtoms[j] == "C":
                    cCount = curAtoms[j]
            for i in range(0, len(useAtoms)):
                if useAtoms[i] in lookup.keys() and useAtoms[i] not in fixed:
                    atomsRange[i] = [int(cCount * lookup[useAtoms[i]][0]), int(cCount * lookup[useAtoms[i]][1] + 1)]
        return useAtoms, atomsRange, curAtoms

    def sixthGoldenRule(self, useAtoms, atomsRange, curAtoms, fixed):

        atoms = ["N", "O", "P", "S"]
        count = [-1, -1, -1, -1]
        pos = [-1, -1, -1, -1]
        for i in range(0, 4):
            for j in range(0, len(useAtoms)):
                if useAtoms[j] == atoms[i]:
                    count[i] = curAtoms[i]

        if "N" in useAtoms:
            if "O" in useAtoms:
                if "P" in useAtoms:
                    if "S" in useAtoms:
                        if min(count) > 1:
                            if pos[0] not in fixed:
                                atomsRange[pos[0]] = [atomsRange[pos[0]][0], 10]
                            if pos[1] not in fixed:
                                atomsRange[pos[1]] = [atomsRange[pos[1]][0], 20]
                            if pos[2] not in fixed:
                                atomsRange[pos[2]] = [atomsRange[pos[2]][0], 4]
                            if pos[3] not in fixed:
                                atomsRange[pos[3]] = [atomsRange[pos[3]][0], 3]
                    else:
                        if min(count[0:3]) > 3:
                            if pos[0] not in fixed:
                                atomsRange[pos[0]] = [atomsRange[pos[0]][0], 11]
                            if pos[1] not in fixed:
                                atomsRange[pos[1]] = [atomsRange[pos[1]][0], 22]
                            if pos[2] not in fixed:
                                atomsRange[pos[2]] = [atomsRange[pos[2]][0], 6]
                elif "S" in useAtoms:
                    if min(count[0], count[1], count[3]) > 6:
                        if pos[0] not in fixed:
                            atomsRange[pos[0]] = [atomsRange[pos[0]][0], 19]
                        if pos[1] not in fixed:
                            atomsRange[pos[1]] = [atomsRange[pos[1]][0], 14]
                        if pos[3] not in fixed:
                            atomsRange[pos[3]] = [atomsRange[pos[3]][0], 8]

            elif "P" in useAtoms:
                if "S" in useAtoms:
                    if min(count[2], count[3], count[0]) > 1:
                        if pos[2] not in fixed:
                            atomsRange[pos[2]] = [atomsRange[pos[2]][0], 3]
                        if pos[3] not in fixed:
                            atomsRange[pos[3]] = [atomsRange[pos[3]][0], 3]
                        if pos[0] not in fixed:
                            atomsRange[pos[0]] = [atomsRange[pos[0]][0], 4]
        elif "O" in useAtoms:
            if "P" in useAtoms:
                if "S" in useAtoms:
                    if min(count[1], count[2], count[3]) > 1:
                        if pos[1] not in fixed:
                            atomsRange[pos[1]] = [atomsRange[pos[1]][0], 14]
                        if pos[2] not in fixed:
                            atomsRange[pos[2]] = [atomsRange[pos[2]][0], 3]
                        if pos[3] not in fixed:
                            atomsRange[pos[3]] = [atomsRange[pos[3]][0], 3]

        return useAtoms, atomsRange

    def perm(self, molMass, massDiff, useAtoms=["C", "H", "O"], atomsRange=[[0, 10], [0, 20], [0, 6]], i=0,
             curAtoms=[0, 0, 0], curMass=0., fixed=[], useSevenGoldenRules=True, useSecondRule=True):
        ret = []

        for c in range(max(0, atomsRange[i][0]), atomsRange[i][1] + 1):
            #use Golden rule number 4 and 5
            if useSevenGoldenRules and useAtoms[i] == "C":
                useAtoms, atomsRange, curAtoms = self.goldenRuleGeneric(useAtoms, atomsRange, i, curAtoms, fixed)
                useAtoms, atomsRange = self.sixthGoldenRule(useAtoms, atomsRange, curAtoms, fixed)
                pass
            curAtoms[i] = c
            t = curMass + self.atoms[useAtoms[i]][0] * c

            if abs(t - molMass) < massDiff and (not useSevenGoldenRules or (self.secondGoldenRule(useAtoms, curAtoms) or not(useSecondRule))):
                ret.append("".join(
                    [str(useAtoms[j]) + str(curAtoms[j] > 1 and curAtoms[j] or "") for j in range(0, len(useAtoms)) if
                     curAtoms[j] != 0]))
            elif t > (molMass + massDiff):
                break
            else:
                if (i + 1) < len(useAtoms):
                    d = self.perm(molMass, massDiff, useAtoms, atomsRange, i + 1, curAtoms, t, fixed,
                                  useSevenGoldenRules=useSevenGoldenRules, useSecondRule=useSecondRule)
                    if len(d) > 0:
                        ret.extend(d)
        curAtoms[i] = 0
        return ret

    #Golden rule number 1
    def firstGoldenRule(self, molMass, useAtoms, atomsRange, fixed):
        assert (len(useAtoms) == len(atomsRange))

        swi = {'C': lambda y: y < 500 and 39 or (y < 1000 and 78 or 162),
               'H': lambda y: y < 500 and 72 or (y < 1000 and 126 or 162),
               'N': lambda y: y < 500 and 20 or (y < 1000 and 20 or 162),
               'O': lambda y: y < 500 and 20 or (y < 1000 and 27 or 162),
               'P': lambda y: y < 500 and 9 or (y < 1000 and 9 or 162),
               'S': lambda y: y < 500 and 10 or (y < 1000 and 14 or 162),
               'F': lambda y: y < 500 and 16 or (y < 1000 and 34 or 162),
               'Cl': lambda y: y < 500 and 10 or (y < 1000 and 12 or 162),
               'Br': lambda y: y < 500 and 4 or (y < 1000 and 8 or 162),
               'Si': lambda y: y < 500 and 8 or (y < 1000 and 14 or 162)}

        for i in range(0, len(useAtoms)):
            if useAtoms[i] in swi and useAtoms[i] not in fixed:
                atomsRange[i][1] = min(swi[useAtoms[i]](molMass), atomsRange[i][1])

        return atomsRange


    #make sure, C is always the first element in useAtoms
    def findFormulas(self, molMass, ppm=5., useAtoms=["C", "H", "O", "N"], atomsRange=[], fixed=[],
                     useSevenGoldenRules=True, useSecondRule=True):

        #use Golden rule number 1
        if len(atomsRange) == 0:
            atomsRange = [[0, 1000] for u in useAtoms]
        if len(atomsRange)!=len(useAtoms):
            raise Exception("useAtoms and atomsRange parameter are not consistent (have a different length)")
        cAtomsRange=[]
        for aR in atomsRange:
            if isinstance(aR, int):
                cAtomsRange.append([aR, aR])
            elif len(aR)==2 and isinstance(aR[0], int) and isinstance(aR[1], int):
                cAtomsRange.append([aR[0], aR[1]])
            else:
                raise Exception("Malformed atomsRange parameter")
        atomsRange=cAtomsRange

        if isinstance(fixed, str):
            ft=formulaTools()
            fixed=ft.parseFormula(fixed).keys()

        if useSevenGoldenRules:
            atomsRange = self.firstGoldenRule(molMass, useAtoms, atomsRange, fixed=fixed)

        res = self.perm(molMass, molMass * (ppm / 1000000.), useAtoms, atomsRange, curAtoms=[t[0] for t in atomsRange],
                        fixed=fixed, useSevenGoldenRules=useSevenGoldenRules, useSecondRule=useSecondRule)

        return res


if __name__=="__main__":
    sfg=SGRGenerator()

    from formulaTools import formulaTools
    ft=formulaTools()




    mz=79.96599


    adducts=[0]

    for add in adducts:
        print add
        m=mz-add

        formsCRes=sfg.findFormulas(m, useAtoms=["C", "N", "H", "O", "P", "S"], atomsRange=[(1,32), (0,500), (0,10000), (0,400), [0, 10], [0, 20]],
                                   useSevenGoldenRules=True, useSecondRule=True, ppm=5.)

        from utils import Bunch, printObjectsAsTable
        from formulaTools import formulaTools
        fT=formulaTools()

        bs=[]
        for f in formsCRes:
            elems=fT.parseFormula(f)
            bs.append(Bunch(formula=f, mass=fT.calcMolWeight(elems), diffPPM=(fT.calcMolWeight(elems)-m)*1000000./m, diffPPMAt300=(fT.calcMolWeight(elems)-m)*1000000./300))

        printObjectsAsTable(bs, attrs=["formula", "mass", "diffPPM", "diffPPMAt300"])



