from pygoslin.parser.Parser import LipidParser
lipid_parser = LipidParser()

lipid_name = "PE 16:1-12:0"
lipid = lipid_parser.parse(lipid_name)

if lipid != None: print(lipid.get_lipid_string())

fas = lipid.lipid.fa

for f in fas:
    print(fas[f].num_carbon, fas[f].double_bonds)