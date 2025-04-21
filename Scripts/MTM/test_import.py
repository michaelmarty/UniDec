from unidec import modules as bp, engine
import unidec.tools as ud
import numpy as np

path = "C:\\Python\\UniDec3\\TestSpectra\\myo.txt"

eng = engine.UniDec()
eng.open_file(path, time_range=[0, 1])
#eng.process_data()
#eng.run_unidec()
eng.unidec_imports(everything=True)
eng.pick_peaks()
eng.dscore()

#plt.plot(eng.data.massdat[:,0], eng.data.massdat[:,1])
#plt.show()

for p in eng.pks.peaks:
    print(p.mass, p.height, p.dscore)

seq = "GLSDGEWQQVLNVWGKVEADIAGHGQEVLIRLFTGHPETLEKFDKFKHLKTEAEMKASEDLKKHGTVVLTALGGILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISDAIIHVLHSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"

expectedmass = bp.calc_pep_mass(seq)

print(expectedmass)

oligos = [[0, expectedmass, 0, 1, "P"], [0, 212, 0, 1, "Biotinylation"]]

oligomasslist, oligonames = ud.make_all_matches(oligos)

print(oligomasslist, oligonames)

matches = ud.match(eng.pks, oligomasslist, oligonames, tolerance=10)

m2 = np.transpose(matches)

for m in m2:
    if m[3] != '':
        print(m)
