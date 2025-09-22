class ConfiguredHeteroAtom:
    def __init__(
        self,
        name="34S",
        mzOffset=1.995796,
        relativeAbundance=0.0443,
        minCount=0,
        maxCount=4,
        entryType="user",
    ):
        self.name = name
        self.mzOffset = mzOffset
        self.relativeAbundance = relativeAbundance
        self.minCount = minCount
        self.maxCount = maxCount
        self.entryType = entryType

    def __str__(self):
        return "ConfiguredHeteroAtom: " + self.name
