from glypy.io import iupac, glycoct, wurcs, linear_code

# GlyTouCan:G08059ZN
glycoct_text = '''RES
1b:b-dgal-HEX-1:5
2s:n-acetyl
3b:b-dglc-HEX-1:5
4s:n-acetyl
5s:sulfate
6b:b-dglc-HEX-1:5
7s:n-acetyl
8b:a-lgal-HEX-1:5|6:d
9b:b-dgal-HEX-1:5
10b:a-dgro-dgal-NON-2:6|1:a|2:keto|3:d
11s:n-acetyl
LIN
1:1d(2+1)2n
2:1o(3+1)3d
3:3d(2+1)4n
4:3o(6+1)5n
5:1o(6+1)6d
6:6d(2+1)7n
7:6o(3+1)8d
8:6o(4+1)9d
9:9o(3+2)10d
10:10d(5+1)11n
'''

reference = structure = glycoct.loads(glycoct_text)
print("%d Monosaccharides in Structure" % (len(structure), ))
iupac_text = iupac.dumps(structure)
print("IUPAC Text =\n%s" % iupac_text)
structure = iupac.loads(iupac_text)
print("Preserved? %r\n" % (reference == structure, ))
wurcs_text = wurcs.dumps(structure)
print("WURCS 2.0 Text =\n%s" % wurcs_text)
structure = wurcs.loads(wurcs_text)
print("Preserved? %r\n" % (reference == structure, ))
linear_code_text = linear_code.dumps(structure)
print("LinearCode Text = %s" % linear_code_text)
structure = linear_code.loads(linear_code_text)
print("Preserved? %r" % (reference == structure, ))
