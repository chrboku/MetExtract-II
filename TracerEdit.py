from formulaTools import getIsotopeMass


class ConfiguredTracer:
    def __init__(
        self,
        name="",
        elementCount=15,
        isotopeA="12C",
        isotopeB="13C",
        enrichmentA=0.9893,
        enrichmentB=0.995,
        amountA=0.9,
        amountB=1.0,
        monoisotopicRatio=0,
        maxRelNegBias=30,
        maxRelPosBias=30,
        tracerType="user",
        id=-1,
        mzDelta=None,
    ):
        self.id = id

        self.name = name  # 0
        self.elementCount = elementCount  # 1
        self.isotopeA = isotopeA  # 2
        self.isotopeB = isotopeB  # 3
        self.enrichmentA = enrichmentA  # 4
        self.enrichmentB = enrichmentB  # 5
        self.amountA = amountA  # 6
        self.amountB = amountB  # 7
        self.monoisotopicRatio = monoisotopicRatio  # 8
        self.maxRelNegBias = maxRelNegBias  # 9
        self.maxRelPosBias = maxRelPosBias  # 10
        self.tracerType = tracerType  # 11

        if mzDelta is None:
            self.mzDelta = (
                getIsotopeMass(self.isotopeB)[0] - getIsotopeMass(self.isotopeA)[0]
            )
        else:
            self.mzDelta = mzDelta

    def __str__(self):
        return (
            "ConfiguredTracer: %s %s %s %d  enrichment: %.3f %.3f amount %.1f %.1f monoisotopicRatio %.3f bias: %.1f %.1f tracerType %s"
            % (
                self.name,
                self.isotopeA,
                self.isotopeB,
                self.elementCount,
                self.enrichmentA,
                self.enrichmentB,
                self.amountA,
                self.amountB,
                self.monoisotopicRatio,
                self.maxRelNegBias,
                self.maxRelPosBias,
                self.tracerType,
            )
        )
