from glypy import monosaccharides
from glypy.algorithms import similarity
from glypy.io import iupac
from glypy.io.nomenclature import identity

Hex = monosaccharides.Hex
Gal = monosaccharides.Gal

print("Gal is Hex = %r" % identity.is_a(Gal, Hex))
matched, total = similarity.monosaccharide_similarity(Gal, Hex)
print("Distance = %d" % (total - matched))
print("Hex is Gal = %r" % identity.is_a(Hex, Gal))
matched, total = similarity.monosaccharide_similarity(Hex, Gal)
print("Distance = %d" % (total - matched))
print("\n")

print("Identify Specially Named Monosaccharides")
Bac = monosaccharides.Bac
print("%s = %s" % (iupac.dumps(Bac), identity.identify(Bac)))
BacNAc = Bac.clone()
BacNAc.add_substituent('n-acetyl', 3)
print("BacNAc is equal to Bac? %r" % (BacNAc == Bac, ))
print("BacNAc is a Bac? %r" % (identity.is_a(BacNAc, Bac, exact=False), ))
print("\n")

print("is-a test for deoxy Hexose")
dHex = iupac.loads("?-?-6-deoxy-Hex")
for spec in ("?-L-Fucp", "?-L-Rhap", "?-D-Quip", "?-D-NeupAc"):
    res = iupac.loads(spec)
    print("Is %s a dHex? %r" % (spec, identity.is_a(res, dHex)))
