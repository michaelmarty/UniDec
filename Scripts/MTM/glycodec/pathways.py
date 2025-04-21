from glypy.enzyme import make_n_glycan_pathway


glycosidases, glycosyltransferases, seeds = make_n_glycan_pathway()

for name, enz in glycosidases.items():
    identifying_information = enz.identifying_information
    if identifying_information is not None:
        description = identifying_information.ec_number, identifying_information.name
    else:
        description = ('', ')')
    print("%s \t %s \t %s" % ((name, ) + description))
print('\n')
for name, enz in glycosyltransferases.items():
    identifying_information = enz.identifying_information
    if identifying_information is not None:
        description = identifying_information.ec_number, identifying_information.name
    else:
        description = ('', ')')
    print("%s \t %s \t %s" % ((name, ) + description))
