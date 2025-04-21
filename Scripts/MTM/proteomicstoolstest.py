import unittest
from Scripts.MTM.proteomicstools import buildSubscriptDict, generateList, sanitizeInput, buildMatches, evalMatchedPeaks
from unidec.IsoDec.match import MatchedPeak


class TestSubscriptDict(unittest.TestCase):
    #normal test case
    def test_build_subscript_regular(self):
        input_data = "AS[10.5]DG[23.3]"
        expected_dict = {'S2': 10.5, 'G4': 23.3}
        result = buildSubscriptDict(input_data)
        self.assertEqual(result, expected_dict)

    #invalid input characters within sequence brackets
    def test_build_subscript_dict_invalid(self):
        input_data = 'AS[SKDZ]DG'
        expected_dict = {}
        result = buildSubscriptDict(input_data)
        self.assertEqual(result, expected_dict)

class TestGenerateList(unittest.TestCase):
    #The actual mass isnt as important here so much as the structure of the list
    #lets just make sure the number and ion type are populating
    #'A', 'B', 'C' will load left to right
    def test_generate_list_abc(self):
        input_data = "AS[10.5]DG[23.3]"
        ion_type = 'c'
        result = generateList(input_data, ion_type)
        expected = ['C1', 'C2', 'C3']
        actual = []
        for i in result:
            res = i.split(":")
            actual.append(res[0])
        self.assertEqual(actual, expected)

    #X, X+1, Y, Y-1, and Z will load right to left
    def test_generate_list_xyz(self):
            input_data = "AS[10.5]DG[23.3]"
            ion_type = 'z'
            result = generateList(input_data, ion_type)
            expected = ['Z3', 'Z2', 'Z1']
            actual = []
            for i in result:
                res = i.split(":")
                actual.append(res[0])
            self.assertEqual(actual, expected)

        #test invalid ion type
    def test_generate_list_invalid_iontype(self):
        input_data = "AS[10.5]DG[23.3]"
        ion_type = '&'
        result = generateList(input_data, ion_type)
        expected = None
        self.assertEqual(result, expected)


class TestSanitizeInput(unittest.TestCase):
    def test_sanitize_input(self):
        input_data = "AS[10.5]DG[23.3]"
        expected = "ASDG"
        result = sanitizeInput(input_data)
        self.assertEqual(result, expected)

class TestBuildMatches(unittest.TestCase):
    #Testing that the matches will actually render from a list of fragment possibilities and a given mass
    def test_build_matches(self):
        expected = ['C2']
        width = 0.05
        ion_list = ['C1:278.6', 'C2:492.699999', 'C3: 748.9']
        results = buildMatches(492.7, ion_list, 0.05)
        self.assertEqual(expected, results)

class TestEvalMatchedPeaks(unittest.TestCase):

    def test_eval_matched_peaks(self):
        # standard testing, expecting a tuple of a list of index match positions,
        #then the number of unique matches
        seq = "AS[10.5]DG[23.3]"
        #seq will give ['C1:88.06366289', 'C2:185.59569129', 'C3:300.62263431']
        matched_peaks = [
            MatchedPeak(z=0, mz=0),
            MatchedPeak(z=0, mz=0),
            MatchedPeak(z=0, mz=0),
        ]
        monoiso_values = [88.06366290, 300.62263432, 145.0293]
        for peak, monoiso in zip(matched_peaks, monoiso_values):
            peak.monoiso = monoiso
            peak.monoisos.append(monoiso)
        expected = ([1,3],2)
        result = evalMatchedPeaks(matched_peaks, seq, 'c')
        self.assertEqual(expected,result)


if __name__ == '__main__':
    unittest.main()
