"""Comprehensive tests for src/formulaTools.py."""

import pytest
from formulaTools import (
    ELECTRON_MASS,
    _parse_charge_suffix,
    formulaTools,
    getElementOfIsotope,
    getIsotopeMass,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def ft():
    return formulaTools()


# ---------------------------------------------------------------------------
# _parse_charge_suffix
# ---------------------------------------------------------------------------


class TestParseChargeSuffix:
    def test_neutral(self):
        assert _parse_charge_suffix("C6H12O6") == ("C6H12O6", 0)

    def test_single_plus(self):
        assert _parse_charge_suffix("C15H20O6+") == ("C15H20O6", 1)

    def test_double_plus(self):
        assert _parse_charge_suffix("C15H20O6++") == ("C15H20O6", 2)

    def test_triple_plus(self):
        assert _parse_charge_suffix("C15H20O6+++") == ("C15H20O6", 3)

    def test_numeric_plus(self):
        assert _parse_charge_suffix("C15H20O6+2") == ("C15H20O6", 2)

    def test_numeric_plus3(self):
        assert _parse_charge_suffix("C15H20O6+3") == ("C15H20O6", 3)

    def test_single_minus(self):
        assert _parse_charge_suffix("C15H20O6-") == ("C15H20O6", -1)

    def test_double_minus(self):
        assert _parse_charge_suffix("C15H20O6--") == ("C15H20O6", -2)

    def test_numeric_minus(self):
        assert _parse_charge_suffix("C15H20O6-2") == ("C15H20O6", -2)

    def test_inline_charge_in_group_not_supported(self):
        # Charge notation embedded inside parentheses (e.g. (C2H4-)3) is not a
        # suffix charge on the whole formula, so it is returned unchanged with charge 0.
        clean, charge = _parse_charge_suffix("(C2H4-)3NO2")
        assert charge == 0
        assert clean == "(C2H4-)3NO2"

    def test_whitespace_stripped(self):
        clean, charge = _parse_charge_suffix("  C6H12O6+  ")
        assert clean == "C6H12O6"
        assert charge == 1

    def test_bare_carbon_cation(self):
        assert _parse_charge_suffix("C+") == ("C", 1)

    def test_bare_carbon_anion(self):
        assert _parse_charge_suffix("C-") == ("C", -1)


# ---------------------------------------------------------------------------
# parseFormula  (basic element parsing)
# ---------------------------------------------------------------------------


class TestParseFormula:
    def test_glucose(self, ft):
        elems = ft.parseFormula("C6H12O6")
        assert elems == {"C": 6, "H": 12, "O": 6}

    def test_alanine(self, ft):
        elems = ft.parseFormula("C3H7NO2")
        assert elems["C"] == 3
        assert elems["H"] == 7
        assert elems["N"] == 1
        assert elems["O"] == 2

    def test_single_atom(self, ft):
        assert ft.parseFormula("C") == {"C": 1}

    def test_parentheses(self, ft):
        elems = ft.parseFormula("C2(H2O)3")
        assert elems["C"] == 2
        assert elems["H"] == 6
        assert elems["O"] == 3

    def test_isotope_notation(self, ft):
        elems = ft.parseFormula("[13C]C5H12O6")
        assert elems["C"] == 6
        assert elems["13C"] == 1

    def test_spaces_ignored(self, ft):
        assert ft.parseFormula("C 6 H 12 O 6") == ft.parseFormula("C6H12O6")

    def test_charge_suffix_stripped(self, ft):
        """parseFormula must silently strip charge notation."""
        assert ft.parseFormula("C15H20O6+") == ft.parseFormula("C15H20O6")
        assert ft.parseFormula("C15H20O6+2") == ft.parseFormula("C15H20O6")
        assert ft.parseFormula("C15H20O6-") == ft.parseFormula("C15H20O6")

    def test_multi_character_element(self, ft):
        elems = ft.parseFormula("NaCl")
        assert elems["Na"] == 1
        assert elems["Cl"] == 1


# ---------------------------------------------------------------------------
# parseFormulaWithCharge
# ---------------------------------------------------------------------------


class TestParseFormulaWithCharge:
    def test_neutral(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C6H12O6")
        assert elems == {"C": 6, "H": 12, "O": 6}
        assert charge == 0

    def test_single_plus(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C15H20O6+")
        assert charge == 1
        assert elems["C"] == 15

    def test_double_plus(self, ft):
        _, charge = ft.parseFormulaWithCharge("C15H20O6++")
        assert charge == 2

    def test_numeric_charge(self, ft):
        _, charge = ft.parseFormulaWithCharge("C15H20O6+2")
        assert charge == 2

    def test_single_minus(self, ft):
        _, charge = ft.parseFormulaWithCharge("C15H20O6-")
        assert charge == -1

    def test_numeric_minus(self, ft):
        _, charge = ft.parseFormulaWithCharge("C15H20O6-2")
        assert charge == -2

    def test_elements_unchanged_with_charge(self, ft):
        elems_charged, _ = ft.parseFormulaWithCharge("C6H12O6+")
        elems_neutral = ft.parseFormula("C6H12O6")
        assert elems_charged == elems_neutral


# ---------------------------------------------------------------------------
# get_formal_charge
# ---------------------------------------------------------------------------


class TestGetFormalCharge:
    @pytest.mark.parametrize(
        "formula, expected",
        [
            ("C6H12O6", 0),
            ("C15H20O6+", 1),
            ("C15H20O6++", 2),
            ("C15H20O6+2", 2),
            ("C15H20O6+++", 3),
            ("C15H20O6+3", 3),
            ("C15H20O6-", -1),
            ("C15H20O6--", -2),
            ("C15H20O6-2", -2),
            ("C-", -1),
            ("C+", 1),
        ],
    )
    def test_charge_values(self, ft, formula, expected):
        assert ft.get_formal_charge(formula) == expected


# ---------------------------------------------------------------------------
# calcMolWeight — neutral
# ---------------------------------------------------------------------------


class TestCalcMolWeightNeutral:
    def test_glucose_mass(self, ft):
        # Monoisotopic mass of C6H12O6 = 180.06339
        elems = ft.parseFormula("C6H12O6")
        mw = ft.calcMolWeight(elems)
        assert abs(mw - 180.06339) < 0.001

    def test_water_mass(self, ft):
        # H2O monoisotopic = 18.01056
        elems = ft.parseFormula("H2O")
        mw = ft.calcMolWeight(elems)
        assert abs(mw - 18.01056) < 0.001

    def test_alanine_mass(self, ft):
        # C3H7NO2 = 89.04768
        elems = ft.parseFormula("C3H7NO2")
        mw = ft.calcMolWeight(elems)
        assert abs(mw - 89.04768) < 0.001

    def test_neutral_charge_zero_unchanged(self, ft):
        elems = ft.parseFormula("C6H12O6")
        assert ft.calcMolWeight(elems) == ft.calcMolWeight(elems, charge=0)


# ---------------------------------------------------------------------------
# calcMolWeight — charged (electron mass adjustment)
# ---------------------------------------------------------------------------


class TestCalcMolWeightCharged:
    def test_cation_loses_one_electron(self, ft):
        """C+ should equal C_neutral - ELECTRON_MASS ≈ 12 - 0.000549."""
        elems = ft.parseFormula("C")
        neutral = ft.calcMolWeight(elems)
        cation = ft.calcMolWeight(elems, charge=1)
        assert abs(cation - (neutral - ELECTRON_MASS)) < 1e-9

    def test_anion_gains_one_electron(self, ft):
        """C- should equal C_neutral + ELECTRON_MASS."""
        elems = ft.parseFormula("C")
        neutral = ft.calcMolWeight(elems)
        anion = ft.calcMolWeight(elems, charge=-1)
        assert abs(anion - (neutral + ELECTRON_MASS)) < 1e-9

    def test_c_plus_approx_11_99945(self, ft):
        """C+ ≈ 12.0 - 0.000549 ≈ 11.99945."""
        elems = ft.parseFormula("C")
        cation_mass = ft.calcMolWeight(elems, charge=1)
        assert abs(cation_mass - 11.99945) < 0.0001

    def test_c_minus_approx_12_00055(self, ft):
        """C- ≈ 12.0 + 0.000549 ≈ 12.00055."""
        elems = ft.parseFormula("C")
        anion_mass = ft.calcMolWeight(elems, charge=-1)
        assert abs(anion_mass - 12.00055) < 0.0001

    def test_doubly_charged_cation(self, ft):
        elems = ft.parseFormula("C15H20O6")
        neutral = ft.calcMolWeight(elems)
        dication = ft.calcMolWeight(elems, charge=2)
        assert abs(dication - (neutral - 2 * ELECTRON_MASS)) < 1e-9

    def test_doubly_charged_anion(self, ft):
        elems = ft.parseFormula("C15H20O6")
        neutral = ft.calcMolWeight(elems)
        dianion = ft.calcMolWeight(elems, charge=-2)
        assert abs(dianion - (neutral + 2 * ELECTRON_MASS)) < 1e-9

    def test_roundtrip_charge_symmetry(self, ft):
        """calcMolWeight(e, +1) + calcMolWeight(e, -1) == 2 * calcMolWeight(e, 0)."""
        elems = ft.parseFormula("C6H12O6")
        neutral = ft.calcMolWeight(elems)
        assert abs(ft.calcMolWeight(elems, 1) + ft.calcMolWeight(elems, -1) - 2 * neutral) < 1e-12


# ---------------------------------------------------------------------------
# Charged formula end-to-end: parseFormulaWithCharge + calcMolWeight
# ---------------------------------------------------------------------------


class TestChargedFormulaEndToEnd:
    def test_c15h20o6_neutral_mass(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C15H20O6")
        assert charge == 0
        mw = ft.calcMolWeight(elems, charge)
        # neutral monoisotopic C15H20O6 = 15*12 + 20*1.007825 + 6*15.994915 ≈ 296.12599
        assert abs(mw - 296.12599) < 0.001

    def test_c15h20o6_singly_charged_cation(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C15H20O6+")
        assert charge == 1
        mw = ft.calcMolWeight(elems, charge)
        # [M]+ = neutral - 1 electron ≈ 296.12544
        assert abs(mw - (296.12599 - ELECTRON_MASS)) < 0.001

    def test_c15h20o6_doubly_charged_from_plusplus(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C15H20O6++")
        assert charge == 2
        mw = ft.calcMolWeight(elems, charge)
        # full ionic mass = neutral - 2 * ELECTRON_MASS
        expected = 296.12599 - 2 * ELECTRON_MASS
        assert abs(mw - expected) < 0.001

    def test_plus2_equals_plusplus(self, ft):
        elems_a, charge_a = ft.parseFormulaWithCharge("C15H20O6+2")
        elems_b, charge_b = ft.parseFormulaWithCharge("C15H20O6++")
        assert charge_a == charge_b == 2
        assert ft.calcMolWeight(elems_a, charge_a) == ft.calcMolWeight(elems_b, charge_b)

    def test_singly_charged_anion(self, ft):
        elems, charge = ft.parseFormulaWithCharge("C15H20O6-")
        assert charge == -1
        mw = ft.calcMolWeight(elems, charge)
        # [M]- = neutral + 1 electron
        neutral = ft.calcMolWeight(ft.parseFormula("C15H20O6"))
        assert abs(mw - (neutral + ELECTRON_MASS)) < 1e-9


# ---------------------------------------------------------------------------
# isIso / getElementFor
# ---------------------------------------------------------------------------


class TestIsotopeHelpers:
    def test_is_iso_true(self, ft):
        assert ft.isIso("13C") is True

    def test_is_iso_false(self, ft):
        assert ft.isIso("C") is False

    def test_is_iso_two_char_element(self, ft):
        assert ft.isIso("Na") is False

    def test_get_element_for_13c(self, ft):
        elem, num = ft.getElementFor("13C")
        assert elem == "C"
        assert num == "13"

    def test_get_element_for_37cl(self, ft):
        elem, num = ft.getElementFor("37Cl")
        assert elem == "Cl"
        assert num == "37"

    def test_get_element_for_requires_isotope(self, ft):
        with pytest.raises(Exception):
            ft.getElementFor("C")


# ---------------------------------------------------------------------------
# flatToString
# ---------------------------------------------------------------------------


class TestFlatToString:
    def test_glucose_roundtrip(self, ft):
        elems = ft.parseFormula("C6H12O6")
        s = ft.flatToString(elems)
        # Re-parse the string and compare element counts
        assert ft.parseFormula(s) == elems

    def test_single_atom_no_subscript(self, ft):
        elems = ft.parseFormula("C")
        assert ft.flatToString(elems) == "C"

    def test_accepts_string_input(self, ft):
        s = ft.flatToString("C6H12O6")
        assert "C6" in s
        assert "H12" in s

    def test_html_subscripts(self, ft):
        elems = ft.parseFormula("C6H12O6")
        s = ft.flatToString(elems, prettyPrintWithHTMLTags=True)
        assert "<sub>" in s
        assert "</sub>" in s

    def test_isotope_bracket_notation(self, ft):
        elems = ft.parseFormula("[13C]C5H12O6")
        s = ft.flatToString(elems)
        assert "[13C]" in s


# ---------------------------------------------------------------------------
# calcDifferenceBetweenSumFormulas
# ---------------------------------------------------------------------------


class TestCalcDifference:
    def test_glucose_minus_water(self, ft):
        diff = ft.calcDifferenceBetweenSumFormulas("H2O", "C6H12O6")
        # Parent C6H12O6 minus fragment H2O → loss H2O
        assert diff == {"C": 6, "H": 10, "O": 5}

    def test_identical_formulas(self, ft):
        diff = ft.calcDifferenceBetweenSumFormulas("C6H12O6", "C6H12O6")
        assert diff == {}

    def test_full_loss(self, ft):
        diff = ft.calcDifferenceBetweenSumFormulas("C6H12O6", "C6H12O6")
        assert len(diff) == 0


# ---------------------------------------------------------------------------
# getIsotopeMass (module-level helper)
# ---------------------------------------------------------------------------


class TestGetIsotopeMass:
    def test_13c_mass(self):
        mass, element = getIsotopeMass("13C")
        assert element == "Carbon"
        assert abs(mass - 13.00335) < 0.0001

    def test_37cl_mass(self):
        mass, element = getIsotopeMass("37Cl")
        assert element == "Chlorine"
        assert abs(mass - 36.96590) < 0.0001

    def test_unknown_returns_minus_one(self):
        mass, element = getIsotopeMass("999X")
        assert mass == -1


# ---------------------------------------------------------------------------
# getElementOfIsotope (module-level helper)
# ---------------------------------------------------------------------------


class TestGetElementOfIsotope:
    def test_13c(self):
        assert getElementOfIsotope("13C") == "C"

    def test_37cl(self):
        assert getElementOfIsotope("37Cl") == "Cl"


# ---------------------------------------------------------------------------
# ELECTRON_MASS constant sanity
# ---------------------------------------------------------------------------


class TestElectronMass:
    def test_value_in_range(self):
        # CODATA electron mass in Da: 5.485799e-4
        assert 5.48e-4 < ELECTRON_MASS < 5.49e-4

    def test_is_positive(self):
        assert ELECTRON_MASS > 0
