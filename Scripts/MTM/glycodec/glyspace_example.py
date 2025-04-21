from glypy.io import glyspace
GlycoRDF = glyspace.GlycoRDF
NSGlycan = glyspace.NSGlycan
NSRDF = glyspace.NSRDF

example = glyspace.get("G08059ZN")
print('\n'.join(repr(example).splitlines()[:10]) + '\n...')
print("Has Primary ID %s" % example.has_primary_id)
seq = (filter(lambda x: 'iupac_extended' in x, example.has_glycosequence))[0]
print("Selected IUPAC Extended Sequence \"%s\"" % seq)
print(seq.has_sequence)

glycosidases = GlycoRDF[:NSRDF.type:NSGlycan.Glycosidase_Treatment]
for g in glycosidases:
    print(g)
