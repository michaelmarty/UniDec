"""LPE rules require supported fragments and retain positive-mode cautions."""
import ast
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import unittest

import pandas as pd
from PublicScripts.Lipids.LipidFunctions import fix_adduct
from tests import test_etherpc_rules


class TestLPERules(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        loader = test_etherpc_rules.TestEtherPCRules
        loader.setUpClass()
        cls.checker = staticmethod(loader.checker)
        path = Path(__file__).resolve().parents[1] / "PublicScripts/Lipids/LipidChecker.py"
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        main = next(n for n in tree.body if isinstance(n, ast.If)
                    and isinstance(n.test, ast.Compare)
                    and isinstance(n.test.left, ast.Name)
                    and n.test.left.id == "__name__")
        start = next(i for i, n in enumerate(main.body)
                     if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == "clean_adduct"
                             for t in n.targets))
        end = next(i for i, n in enumerate(main.body)
                   if isinstance(n, ast.Expr) and isinstance(n.value, ast.Call)
                   and isinstance(n.value.func, ast.Attribute)
                   and n.value.func.attr == "to_excel")
        cls.final_checks = compile(ast.Module(body=main.body[start:end], type_ignores=[]),
                                   str(path), "exec")

    def check(self, sample):
        rows = pd.DataFrame([{
            "Ontology": "LPE", "Adduct type": "[M+H]+",
            "Metabolite name": "LPE 18:1", "Reference m/z": 480.30847,
            "Tolerance": 0.05, "T1": "18:1",
            "MS/MS spectrum": sample, "Ref Spec": "339.2894:999",
        }])
        with redirect_stdout(StringIO()):
            result = self.checker(rows)
        namespace = {"df": result, "fix_adduct": fix_adduct}
        exec(self.final_checks, namespace)
        return namespace["df"].iloc[0]

    def test_headgroup_rule_is_accepted_with_chain_caution(self):
        row = self.check("339.2894:999")
        self.assertIn("HGNL(141): Yes", row["Matched Required Fragments"])
        self.assertTrue(row["Overall Match Status"])
        self.assertEqual(row["Validation Status"], "Accepted")
        self.assertIn("acyl chain unconfirmed", row["Identification Status"])

    def test_missing_headgroup_is_failed(self):
        row = self.check("308.296:999")
        self.assertFalse(row["Overall Match Status"])
        self.assertEqual(row["Validation Status"], "Failed")
        self.assertEqual(row["Identification Status"], "Failed")

    def check_negative(self, sample, reference):
        rows = pd.DataFrame([{
            "Ontology": "LPE", "Adduct type": "[M-H]-",
            "Metabolite name": "LPE 18:1", "Reference m/z": 478.29392,
            "Tolerance": 0.05, "T1": "18:1",
            "MS/MS spectrum": sample, "Ref Spec": reference,
        }])
        with redirect_stdout(StringIO()):
            return bool(self.checker(rows).iloc[0]["Overall Match Status"])

    def test_negative_196_alternative(self):
        self.assertTrue(self.check_negative("196.0376:100 281.248:999",
                                            "196.0376:100 281.248:999"))

    def test_negative_214_alternative(self):
        self.assertTrue(self.check_negative("214.0482:100 281.248:999",
                                            "214.0482:100 281.248:999"))

    def test_negative_headgroup_requires_chain_fragment(self):
        self.assertFalse(self.check_negative("196.0376:100 214.0482:999",
                                             "196.0376:100 214.0482:100 281.248:999"))

    def test_negative_chain_alone_fails(self):
        self.assertFalse(self.check_negative("281.248:999",
                                             "196.0376:100 214.0482:100 281.248:999"))

    def test_negative_alternative_requires_reference_support(self):
        self.assertFalse(self.check_negative("214.0482:100 281.248:999",
                                             "196.0376:100 281.248:999"))

    def test_figure5_excludes_provisional_and_unknown_labels(self):
        path = (Path(__file__).resolve().parents[3] / "Scripts/Annika"
                / "DDA Paper Figures/Figure6.py")
        tree = ast.parse(path.read_text(encoding="utf-8-sig"))
        start = next(i for i, n in enumerate(tree.body)
                     if isinstance(n, ast.Assign)
                     and any(isinstance(t, ast.Name) and t.id == "true_values"
                             for t in n.targets))
        end = next(i for i, n in enumerate(tree.body[start:], start)
                   if isinstance(n, ast.Assign)
                   and any(isinstance(t, ast.Subscript) and isinstance(t.slice, ast.Constant)
                           and t.slice.value == "ID Status" for t in n.targets))
        rows = pd.DataFrame({"Overall Match Status": [True, False, "Provisional", None, False],
                             "Validation Status": ["Accepted", "Failed", "Provisional", None,
                                                   "Provisional"]})
        namespace = {"df": rows}
        with redirect_stdout(StringIO()):
            exec(compile(ast.Module(body=tree.body[start:end], type_ignores=[]),
                         str(path), "exec"), namespace)
        self.assertEqual(namespace["df"]["Overall Match Status"].tolist(), [True, False])


if __name__ == "__main__":
    unittest.main()
