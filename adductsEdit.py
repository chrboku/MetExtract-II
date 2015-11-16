
class ConfiguredAdduct():
    def __init__(self, name="", mzoffset=1.99705, charge=1, polarity="+", mCount=1, entryType="user"):
        self.name = name
        self.mzoffset = mzoffset
        self.charge = charge
        self.polarity = polarity
        self.mCount = mCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredAdduct %s %.5f %s %d" % (self.name, self.mzoffset, self.polarity, self.charge)


class ConfiguredElement():
    def __init__(self, name="", weight=12., numberValenzElectrons=1, minCount=1, maxCount=1, entryType="user"):
        self.name = name
        self.weight = weight
        self.numberValenzElectrons = numberValenzElectrons
        self.minCount = minCount
        self.maxCount = maxCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredElement %s %.5f %d %d %d" % (
            self.name, self.weight, self.numberValenzElectrons, self.minCount, self.maxCount)
