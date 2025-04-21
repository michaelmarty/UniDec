import MDAnalysis as mda
import os


topdir = "D:\\VMD\\AqpZMutants\\R224A_CG"
os.chdir(topdir)

u = mda.Universe('step5_assembly.psf', 'step7.1_production.xtc')
#u = mda.Universe('step7.1_production.tpr', 'step7.1_production.xtc')

protein = u.select_atoms('segid PRO*')
protein.write('wrapped.gro')
print(len(protein))
for ts in u.trajectory:
    protein.unwrap(compound='residues')

protein.write('wrapped.xtc', frames='all')
protein.write('wrapped.gro')