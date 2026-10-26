"""Exercise the real rule table and checker without optional plotting imports."""
import ast
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import re
import unittest

import numpy as np
import pandas as pd

from PublicScripts.Lipids import LipidFragPrediction as prediction
from PublicScripts.Lipids.LipidFunctions import (
    fix_adduct, match_frag_to_spec, parse_spec_string, spec_to_nl,
)


class TestEtherPCRules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        path = Path(__file__).resolve().parents[1] / "PublicScripts/Lipids/LipidChecker.py"
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        rule_text = next(ast.literal_eval(node.value) for node in tree.body
                         if isinstance(node, ast.Assign)
                         and any(isinstance(t, ast.Name) and t.id == "lipid_data"
                                 for t in node.targets))
        cls.rules = pd.read_csv(StringIO(rule_text), sep="\t", index_col=False)
        for column in cls.rules.select_dtypes(include="object"):
            cls.rules[column] = cls.rules[column].str.strip()
        namespace = dict(vars(prediction))
        namespace.update(pd=pd, np=np, re=re, lipid_df=cls.rules,
                         fix_adduct=fix_adduct, match_frag_to_spec=match_frag_to_spec,
                         parse_spec_string=parse_spec_string, spec_to_nl=spec_to_nl)
        names = {"get_hg_frag", "get_tail_frag", "get_match_result",
                 "check_disqualifying_fragments", "check_fragments"}
        functions = ast.Module(body=[node for node in tree.body
                                     if isinstance(node, ast.FunctionDef)
                                     and node.name in names], type_ignores=[])
        exec(compile(functions, str(path), "exec"), namespace)
        cls.checker = staticmethod(namespace["check_fragments"])

    def check_spectra(self, sample, reference=None):
        # Published library peaks for PC O-14:0_14:0, supplied independently
        # of the predictor so the test checks the neutral-loss mass arithmetic.
        def spectrum(peaks):
            return " ".join(f"{mz}:100" for mz in peaks)
        row = {"Ontology": "EtherPC", "Adduct type": "[M+CH3COO]-",
               "Metabolite name": "PC O-28:0|PC O-14:0_14:0",
               "Reference m/z": 722.53412, "Tolerance": 0.05,
               "T1": "O-14:0", "T2": "14:0",
               "MS/MS spectrum": spectrum(sample),
               "Ref Spec": spectrum(sample if reference is None else reference)}
        with redirect_stdout(StringIO()):
            result = self.checker(pd.DataFrame([row]))
        return bool(result.loc[0, "Overall Match Status"])

    def test_acid_loss_alternative(self):
        self.assertTrue(self.check_spectra([227.2015, 420.2880]))

    def test_ketene_loss_alternative(self):
        self.assertTrue(self.check_spectra([227.2015, 438.2986]))

    def test_old_two_peaks_are_insufficient(self):
        self.assertFalse(self.check_spectra([227.2015, 648.4968]))

    def test_complementary_products_need_acyl_fragment(self):
        self.assertFalse(self.check_spectra([420.2880, 438.2986]))

    def test_reference_must_support_the_same_alternative(self):
        self.assertFalse(self.check_spectra([227.2015, 420.2880],
                                          [227.2015, 438.2986]))

    def test_two_fragments_per_negative_rule(self):
        rules = self.rules[(self.rules["Class"] == "EtherPC")
                           & (self.rules["Adduct"] == "[M+CH3COO]")]
        self.assertEqual(len(rules), 2)
        self.assertTrue(all(len(rule.split(",")) == 2
                            for rule in rules["Required Frags"]))


if __name__ == "__main__":
    unittest.main()
