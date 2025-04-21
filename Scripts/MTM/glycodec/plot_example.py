import matplotlib
matplotlib.use("agg")
from glypy import plot
from glypy.io import glycoct

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

structure = glycoct.loads(glycoct_text)
draw_tree, ax = plot.plot(structure, label=True)
ax.set_xlim(-4, 0.6)
ax.set_ylim(-0.8, 1.65)
ax.figure.savefig("G08059ZN.pdf", bbox_inches="tight")
