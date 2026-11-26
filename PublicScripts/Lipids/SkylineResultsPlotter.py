"""Quantify Skyline lipid exports and compare replicate groups.

Products excludes named precursors and renamed rows at or above precursor
m/z (0.005 m/z tolerance, same charge when supplied). This also rejects
over-shifted labelled precursor/isotope rows in re-imported standard lists.
For PC/EtherPC +H, Products uses only the phosphocholine headgroup: m/z 184.0733
for both analytes and glycerol-d5 PC standards (0.005 m/z window). The glycerol
label does not shift this fragment. Previously exported standard transitions at
189.1047 require corrected Skyline transitions and re-extraction, not relabeling.
Acetate still uses the existing product sum. Selection uses measured Product Mz,
so it works with both original fragment names and re-imported "Ion [...]" names.
PC/EtherPC uses acetate when quantifiable in a replicate, otherwise +H, after
normalizing each adduct against its matching standard adduct. Quant Adduct
columns record the choice. Standards with >5-fold raw-area ranges within a group
are retained and flagged for review; different standard species use their own
catalogue concentrations. Analyte/standard area ratios >50 are also warnings,
not rejection thresholds. Neither is a calibrated assay acceptance limit.
Standards with missing/nonpositive areas cannot serve as denominators.
TIC normalization remains
optional and produces relative lipid composition rather than absolute abundance.
compare_pipeline keeps the union of input lipids without lipid/class exclusion
lists. Each run writes <input-stem>_quantification_notes.csv to record warnings
and unavailable quantification; every lipid remains in the comparison tables.
With missing_set2_as_zero enabled, input entries absent from set2 for a lipid
present in the control become zero detected signal. Quantification failures
remain unavailable. Tests using zeros explicitly assume those zero signals;
infinite fold changes appear as labelled boundary triangles in the volcano.
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib as mpl
import io
from PublicScripts.Lipids.LipidFunctions import simplify_class_name, get_color_for_class
from copy import deepcopy
from matplotlib.gridspec import GridSpec

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['errorbar.capsize'] = 3
mpl.rcParams['patch.force_edgecolor'] = True
mpl.rcParams['patch.facecolor'] = 'b'
mpl.rcParams['lines.markersize'] = 7
mpl.rcParams['font.size'] = 12

mpl.use("WxAgg")

israw = r"""
Class,Lipid,Name 2,Concentration (ug/mL),m/z (+),m/z (-),RT,Product Number,SMILES
CE,14:1 cholesteryl-d7 ester,IS SE 27:1(+[2]H7) 14:1,25,602.58,600.58,10,700220,[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCCC/C=C\CCCC)=O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
CE,16:1 cholesteryl-d7 ester,IS SE 27:1(+[2]H7) 16:1,50,630.61,628.61,11,700221,[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCCC/C=C\CCCCCC)=O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
CE,18:1 cholesteryl-d7 ester,IS SE 27:1(+[2]H7) 18:1,75,658.64,656.64,11.7,700222,[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCCC/C=C\CCCCCCCC)=O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
CE,20:3 cholesteryl-d7 ester,IS SE 27:1(+[2]H7) 20:3,50,682.64,680.64,12,700223,[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
CE,22:4 cholesteryl-d7 ester,IS SE 27:1(+[2]H7) 22:4,25,708.66,706.66,12,700226,[H][C@@]12[C@]([C@](CC[C@H](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
Cer,C24:1 Ceramide-d7 (d18:1-d7/24:1),IS Cer 18:1;2(+[2]H7)/24:1,75,655.67,653.67,12.4,860679,OC[C@]([H])(NC(CCCCCCCCCCCCC/C=C\CCCCCCCC)=O)[C@@](O)([H])/C=C/CCCCCCCCCCC(C(C([2H])([2H])[2H])([2H])[2H])([2H])[2H]
Cer,C22:1 Ceramide-d7 (d18:1-d7/22:1),IS Cer 18:1;2(+[2]H7)/22:1,50,627.63,625.63,11,860745,[H][C@](/C=C/CCCCCCCCCCC(C([2H])(C([2H])([2H])[2H])[2H])([2H])[2H])(O)[C@@]([H])(NC(CCCCCCCCCCC/C=C\CCCCCCCC)=O)CO
Cer,C20:1 Ceramide-d7 (d18:1-d7/20:1),IS Cer 18:1;2(+[2]H7)/20:1,25,599.6,597.6,10,860746,[H][C@](/C=C/CCCCCCCCCCC(C([2H])(C([2H])([2H])[2H])[2H])([2H])[2H])(O)[C@@]([H])(NC(CCCCCCCCC/C=C\CCCCCCCC)=O)CO
Cer,C18:1 Ceramide-d7 (d18:1-d7/18:1),IS Cer 18:1;2(+[2]H7)/18:1,50,571.57,569.57,9,860747,[H][C@](/C=C/CCCCCCCCCCC(C([2H])(C([2H])([2H])[2H])[2H])([2H])[2H])(O)[C@@]([H])(NC(CCCCCCC/C=C\CCCCCCCC)=O)CO
Cer,C16:1 Ceramide-d7 (d18:1-d7/16:1),IS Cer 18:1;2(+[2]H7)/16:1,75,543.54,541.54,8,860748,[H][C@](/C=C/CCCCCCCCCCC(C([2H])(C([2H])([2H])[2H])[2H])([2H])[2H])(O)[C@@]([H])(NC(CCCCCCC/C=C\CCCCCC)=O)CO
CL,18:2-18:2-18:2-18:2 CL-d5,IS CL(+[2]H5) 18:2_18:2_18:2_18:2,15,1454.9,1452.9963,10.6,791108,O=P([O-])(OC[C@@]([H])(OC(CCCCCCC/C=C\C/C=C\CCCCC)=O)COC(CCCCCCC/C=C\C/C=C\CCCCC)=O)OC([2H])([2H])C([2H])(O)C([2H])([2H])OP([O-])(OC[C@]([H])(OC(CCCCCCC/C=C\C/C=C\CCCCC)=O)COC(CCCCCCC/C=C\C/C=C\CCCCC)=O)=O.[NH4+].[NH4+]
DG,17:0-22:4 DG-d5,IS DG(+[2]H5) 17:0_22:4,25,664.59,662.59,10.9,800823,[2H][C@](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)(C(O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
DG,17:0-18:1 DG-d5,IS DG(+[2]H5) 17:0_18:1,75,614.57,612.57,11,800824,OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H]
DG,17:0-20:3 DG-d5,IS DG(+[2]H5) 17:0_20:3,50,638.57,636.57,10.7,800825,OC([C@]([2H])(OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H]
DG,17:0-16:1 DG-d5,IS DG(+[2]H5) 16:1_17:0,50,586.54,584.54,10.2,800826,OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H]
DG,17:0-14:1 DG-d5,IS DG(+[2]H5) 14:1_17:0,25,558.51,556.51,9,800827,OC([C@]([2H])(OC(CCCCCCC/C=C\CCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H]
LPA,18:1 Lyso PA-d5, IS LPA 18:1(+[2]H9),25,446.3, 444.3, 1.75, 793074, [2H]C([2H])([2H])C([2H])([2H])C([2H])([2H])CCCC/C=C\CCCCCCC(=O)OC[C@H](O)COP([O-])([O-])=O.[Na+]
LPC,17:0 Lyso PC-d5,IS LPC 17:0(+[2]H5),50,492.38,490.38,3,855679,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPC,19:0 Lyso PC-d5,IS LPC 19:0(+[2]H5),25,520.41,518.41,4,855778,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPC,15:0 Lyso PC-d5,IS LPC 15:0(+[2]H5),25,464.35,462.35,2,870309,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPE,15:0 Lyso PE-d5,IS LPE(+[2]H5) 15:0,25,422.3,420.3,2.5,856709,[O-]P(OCC[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPE,17:0 Lyso PE-d5,IS LPE(+[2]H5) 17:0,50,450.33,448.33,3.2,856710,[O-]P(OCC[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPE,19:0 Lyso PE-d5,IS LPE(+[2]H5) 19:0,25,478.36,476.36,4.3,856716,[O-]P(OCC[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPG,15:0 Lyso PG-d5,IS LPG(+[2]H5) 15:0,25,476.3,474.3,2,858123,[O-]P(OCC(O)CO)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPG,19:0 Lyso PG-d5,IS LPG(+[2]H5) 19:0,25,532.4,530.4,3,858129,[O-]P(OCC(O)CO)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPG,17:0 Lyso PG-d5,IS LPG(+[2]H5) 17:0,50,504.3,502.3,2.5,858130,[O-]P(OCC(O)CO)(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPI,19:0 Lyso PI-d5,IS LPI(+[2]H5) 19:0,25,620.4,618.4,2.9,850106,[2H][C@](O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCCCC)=O)([2H])[2H]
LPI,15:0 Lyso PI-d5,IS LPI(+[2]H5) 15:0,25,564.3,562.3,1.95,850107,[2H][C@](O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCC)=O)([2H])[2H]
LPI,17:0 Lyso PI-d5,IS LPI(+[2]H5) 17:0,50,592.4,590.4,2.3,850108,[2H][C@](O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
LPS,15:0 Lyso PS-d5,IS LPS(+[2]H5) 15:0,25,489.3,487.3,1.9,858146,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPS,19:0 Lyso PS-d5,IS LPS(+[2]H5) 19:0,25,545.4,543.4,3,858147,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
LPS,17:0 Lyso PS-d5,IS LPS(+[2]H5) 17:0,50,517.3,515.3,2.3,858148,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PA,15:0-18:1 PA-d7,IS PA 15:0_18:1(+[2]H7),50,668.5, 666.5, 5.6, 791642, CCCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])([O-])=O)OC(=O)CCCCCCC/C=C\CCCCCC([2H])([2H])C([2H])([2H])C([2H])([2H])[2H]
PC,17:0-22:4 PC-d5,IS PC(+[2]H5) 17:0_22:4,50,829.64,827.64,9.1,855678,[2H][C@](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)(C(OP(OCC[N+](C)(C)C)([O-])=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PC,17:0-20:3 PC-d5,IS PC(+[2]H5) 17:0_20:3,100,803.62,801.62,8.9,855680,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PC,17:0-18:1 PC-d5,IS PC(+[2]H5) 17:0_18:1,150,779.62,777.62,9.2,855681,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PC,17:0-16:1 PC-d5,IS PC(+[2]H5) 16:1_17:0,100,751.59,749.59,8.5,855682,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PC,17:0-14:1 PC-d5,IS PC(+[2]H5) 14:1_17:0,50,723.56,721.56,7.5,855683,[O-]P(OCC[N+](C)(C)C)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PE,17:0-22:4 PE-d5,IS PE(+[2]H5) 17:0-22:4,25,787.59,785.59,9.3,856717,[2H][C@](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)(C(OP(OCC[NH3+])([O-])=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PE,17:0-20:3 PE-d5,IS PE(+[2]H5) 17:0-20:3,50,761.58,759.58,9.1,856718,[O-]P(OCC[NH3+])(OC([C@]([2H])(OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PE,17:0-18:1 PE-d5,IS PE(+[2]H5) 17:0-18:1,75,737.58,735.58,9.6,856719,[O-]P(OCC[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PE,17:0-16:1 PE-d5,IS PE(+[2]H5) 16:1-17:0,50,709.55,707.55,8.7,856720,[O-]P(OCC[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PE,17:0-14:1 PE-d5,IS PE(+[2]H5) 14:1-17:0,25,681.52,679.52,7.7,856721,[O-]P(OCC[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PG,17:0-22:4 PG-d5,IS PG(+[2]H5) 17:0-22:4,25,818.6,816.6,7.6,858131,[2H][C@](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)(C(OP(OCC(O)CO)([O-])=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PG,17:0-20:3 PG-d5,IS PG(+[2]H5) 17:0-20:3,50,792.6,790.6,7.4,858132,[O-]P(OCC(O)CO)(OC([C@]([2H])(OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PG,17:0-18:1 PG-d5,IS PG(+[2]H5) 17:0-18:1,75,768.5,766.5,7.8,858133,[O-]P(OCC(O)CO)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PG,17:0-16:1 PG-d5,IS PG(+[2]H5) 16:1-17:0,50,740.5,738.5,7.1,858134,[O-]P(OCC(O)CO)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PG,17:0-14:1 PG-d5,IS PG(+[2]H5) 14:1-17:0,25,712.5,710.5,6.3,858135,[O-]P(OCC(O)CO)(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PI,17:0-14:1 PI-d5,IS PI(+[2]H5) 14:1-17:0,25,800.5,798.5,6.1,850109,[2H][C@](OC(CCCCCCC/C=C\CCCC)=O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PI,17:0-16:1 PI-d5,IS PI(+[2]H5) 16:1-17:0,50,828.6,826.6,6.9,850110,[2H][C@](OC(CCCCCCC/C=C\CCCCCC)=O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PI,17:0-18:1 PI-d5,IS PI(+[2]H5) 17:0-18:1,75,856.6,854.6,7.6,850111,[2H][C@](OC(CCCCCCC/C=C\CCCCCCCC)=O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PI,17:0-20:3 PI-d5,IS PI(+[2]H5) 17:0-20:3,50,880.6,878.6,7.3,850112,[2H][C@](OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PI,17:0-22:4 PI-d5,IS PI(+[2]H5) 17:0-22:4,25,906.6,904.6,7.5,850118,[2H][C@](OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)(C(OP([O-])(O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)=O)([2H])[2H])C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H]
PS,17:0-22:4 PS-d5,IS PS(+[2]H5) 17:0-22:4,25,831.6,829.6,7.6,858149,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(OC(CCCCC/C=C\C/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PS,17:0-20:3 PS-d5,IS PS(+[2]H5) 17:0-20:3,50,805.6,803.6,7.4,858150,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(OC(CCCCCC/C=C\C/C=C\C/C=C\CCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PS,17:0-18:1 PS-d5,IS PS(+[2]H5) 17:0-18:1,75,781.6,779.6,7.7,858151,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PS,17:0-16:1 PS-d5,IS PS(+[2]H5) 16:1-17:0,50,753.5,751.5,7,858152,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
PS,17:0-14:1 PS-d5,IS PS(+[2]H5) 14:1-17:0,25,725.5,723.5,6.2,858153,[O-]P(OC[C@](C([O-])=O)([H])[NH3+])(OC([C@]([2H])(OC(CCCCCCC/C=C\CCCC)=O)C(OC(CCCCCCCCCCCCCCCC)=O)([2H])[2H])([2H])[2H])=O
SM,18:1 SM (d18:1/18:1)-d9,IS SM(+[2]H9) 18:1;2/18:1,50,738.64,736.64,8,860740,[H][C@](/C=C/CCCCCCCCCCCCC)(O)[C@@]([H])(NC(CCCCCCC/C=C\CCCCCCCC)=O)COP([O-])(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])=O
SM,16:1 SM (d18:1/16:1)-d9,IS SM(+[2]H9) 18:1;2/16:1,75,710.61,708.61,7.2,860741,[H][C@](/C=C/CCCCCCCCCCCCC)(O)[C@@]([H])(NC(CCCCCCC/C=C\CCCCCC)=O)COP([O-])(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])=O
SM,20:1 SM (d18:1/20:1)-d9,IS SM(+[2]H9) 18:1;2/20:1,25,766.67,764.67,8.9,860742,[H][C@](/C=C/CCCCCCCCCCCCC)(O)[C@@]([H])(NC(CCCCCCCCC/C=C\CCCCCCCC)=O)COP([O-])(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])=O
SM,22:1 SM (d18:1/22:1)-d9,IS SM(+[2]H9) 18:1;2/22:1,50,794.7,792.7,9.8,860743,[H][C@](/C=C/CCCCCCCCCCCCC)(O)[C@@]([H])(NC(CCCCCCCCCCC/C=C\CCCCCCCC)=O)COP([O-])(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])=O
SM,24:1 SM (d18:1/24:1)-d9,IS SM(+[2]H9) 18:1;2/24:1,75,822.73,820.73,10,860744,[H][C@](/C=C/CCCCCCCCCCCCC)(O)[C@@]([H])(NC(CCCCCCCCCCCCC/C=C\CCCCCCCC)=O)COP([O-])(OCC[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H])=O
ST,d7-cholesterol,IS ST 27:1;1(+[2]H7),10,376.5,374.5,7,700041,[H][C@@]12[C@]([C@](CC[C@H](O)C3)(C)C3=CC2)([H])CC[C@@]4(C)[C@@]1([H])CC[C@]4([H])[C@]([H])(C)CCCC([2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]
TG,14:0-13:0-14:0 TG-d5,IS TG(+[2]H5) 13:0-14:0-14:0,25,714.66,712.66,15.6,860906,CCCCCCCCCCCCCC(OC([2H])([2H])C([2H])(OC(CCCCCCCCCCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCC)=O)=O
TG,14:0-15:1-14:0 TG-d5,IS TG(+[2]H5) 14:0-14:0-15:1,50,740.67,738.67,15.6,860907,CCCCCCCCCCCCCC(OC([2H])([2H])C([2H])(OC(CCCCCCCC/C=C\CCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCC)=O)=O
TG,14:0-17:1-14:0 TG-d5,IS TG(+[2]H5) 14:0-14:0-17:1,75,768.71,766.71,15.9,860908,CCCCCCCCCCCCCC(OC([2H])([2H])C([2H])(OC(CCCCCCCC/C=C\CCCCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCC)=O)=O
TG,16:0-17:1-16:0 TG-d5,IS TG(+[2]H5) 16:0-16:0-17:1,125,824.77,822.77,16.5,860909,[2H]C(C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCCC/C=C\CCCCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O
TG,16:0-15:1-16:0 TG-d5,IS TG(+[2]H5) 15:1-16:0-16:0,100,796.74,794.74,16.3,860910,[2H]C(C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCCC/C=C\CCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O
TG,16:0-19:2-16:0 TG-d5,IS TG(+[2]H5) 16:0-16:0-19:2,100,850.78,848.78,16.6,860911,[2H]C(C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O)(OC(CCCCCCCC/C=C\C/C=C\CCCCC)=O)C([2H])([2H])OC(CCCCCCCCCCCCCCC)=O
TG,18:1-19:2-18:1 TG-d5,IS TG(+[2]H5) 18:1-18:1-19:2,50,902.81,900.81,16.6,860912,[2H]C(C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O)(OC(CCCCCCCC/C=C\C/C=C\CCCCC)=O)C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O
TG,18:1-21:2-18:1 TG-d5,IS TG(+[2]H5) 18:1-18:1-21:2,25,930.85,928.85,16.8,860913,[2H]C(C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O)(OC(CCCCCCCCCC/C=C\C/C=C\CCCCC)=O)C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O
TG,18:1-17:1-18:1 TG-d5,IS TG(+[2]H5) 17:1-18:1-18:1,75,876.8,874.8,16.5,860914,[2H]C(C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O)(OC(CCCCCCCC/C=C\CCCCCC)=O)C([2H])([2H])OC(CCCCCCC/C=C\CCCCCCCC)=O
"""

# Load the inline CSV into a DataFrame for quick use
default_isdf = pd.read_csv(io.StringIO(israw))
# Rename "Lipid" column to "Name 1" and "Name 2" to "Lipid"
default_isdf = default_isdf.rename(columns={"Lipid": "Name 1", "Name 2": "Lipid"})
# Drop "IS " prefix from "Lipid" column
default_isdf["Lipid"] = default_isdf["Lipid"].str.strip().str.replace(r"^IS\s+", "", regex=True)
# Replace - in "Lipid" column with _
default_isdf["Lipid"] = default_isdf["Lipid"].str.replace("-", "_")
# Replace TAG with TG and DAG with DG in "Lipid" column
default_isdf["Lipid"] = default_isdf["Lipid"].str.replace("TAG", "TG")
default_isdf["Lipid"] = default_isdf["Lipid"].str.replace("DAG", "DG")
default_isdf["Lipid"] = default_isdf["Lipid"].str.replace("SHexCer", "HexCer")
default_isdf["Lipid"] = default_isdf["Lipid"].str.replace("Hex2Cer", "HexCer")
print(default_isdf.to_string())

# Drop last two columns
# default_isdf = default_isdf.iloc[:, :-2]
# Calculated mass from m/z avg
default_isdf["Mass"] = (default_isdf["m/z (+)"] + default_isdf["m/z (-)"]) / 2.
# Calculate conc in uM from ug/mL using mass
default_isdf["Concentration (uM)"] = default_isdf["Concentration (ug/mL)"] / default_isdf["Mass"]

default_is_mapper = {"GM3": "PI", "PEtOH": "PG", "SHexCer": "Cer", "EtherPC": "PC", "HexCer": "Cer", "Hex2Cer": "Cer"}
area_match_is_classes = {"PC", "EtherPC", "EtherPE","PE", "PI", "SM"}
preferred_quant_adducts = {"PC": "[M+CH3COO]", "EtherPC": "[M+CH3COO]"}
# QC compares replicate areas within experimental groups, never between groups.
IS_AREA_VARIATION_WARNING_FOLD = 5.0
PC_HEADGROUP_MZ = 184.0733
PC_HEADGROUP_TOLERANCE = 0.005
QUANTIFICATION_NOTE_COLUMNS = ["Input", "Molecule", "Class", "Adduct", "Reason", "Stage",
                   "Rows removed", "Row unit", "Replicates"]
# print(default_isdf.to_string())
# exit()

def safedivide(a, b):
    """
    Divide numpy arrays b from a. Define a/0=0 to avoid zero division error.
    :param a: Numerator
    :param b: Denominator
    :return: a/b
    """
    c = deepcopy(b)
    c[b != 0] = a[b != 0] / b[b != 0]
    return c

def get_row_dict(id, mdf):
    row = mdf[mdf["Alignment ID"] == id]
    if len(row) == 0:
        return None
    return row.iloc[0].to_dict()

def translate_to_skyline(vdf, mdf):
    """
    Translates a DataFrame from MS-DIAL format to Skyline format. Takes dfs from PeakValue and PeakMaster files.
    """
    rows = []
    for i, row in vdf.iterrows():
        peakid = row["ID"]
        rowdict = get_row_dict(peakid, mdf)

        name = rowdict["Metabolite name"] if rowdict is not None else f"Unknown_{peakid}"
        listname = rowdict["Ontology"] if rowdict is not None else "Unknown"
        fragmentname = "None"
        area = row["Area"]
        replicate = row["File"]

        rows.append({"Molecule List Name": listname,
                          "Molecule": name,
                          "Fragment Name": fragmentname,
                          "Area": area,
                          "Replicate Name": replicate})

    # Create a new DataFrame for Skyline
    sdf = pd.DataFrame(rows)
    return sdf

def clean_classes(df):
    # Switch any DAG to DG and TAG to TG in "Molecule List Name"
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("DAG", "DG")
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("TAG", "TG")
    # Switch "SE 27:1 to CE"
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("SE 27:1", "CE")
    # Switch any lipids with "SE" in the name to CE
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("SE", "CE")
    # Switch name ST 27:1;1 to ST
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("ST 27:1;1", "ST")
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("SHexCer", "HexCer")
    df["Molecule List Name"] = df["Molecule List Name"].str.replace("Hex2Cer", "HexCer")
    return df

def clean_names(names, split_side=1):
    # Split at the | character and take the split_side part (0 for left, 1 for right)
    # If "|" is not present, keep the original name
    names = names.str.split("|").str[split_side].fillna(names)
    # names = names.str.split("|").str[split_side]
    # Replace any - with _ in "Molecule" column
    names = names.str.replace("-", "_")
    return names

def fix_name_order(names):
    fixednames = []
    for name in names:
        if not isinstance(name, str):
            print(name)
            fixednames.append(name)
            continue
        if "_" in name:
            p1 = name.split(" ")[0]
            p2 = name.split(" ")[1]
            # If "P" or "O" in p2 just add original name and continue
            if "P" in p2 or "O" in p2:
                fixednames.append(name)
                continue
            # For tails in p2, split at the _ and then order them in ascending order and join with a _
            tails = p2.split("_")
            tails = sorted(tails, key=lambda x: int(x.split(":")[0]))
            p2 = "_".join(tails)
            fixednames.append(f"{p1} {p2}")
        else:
            fixednames.append(name)
    return fixednames

def split_names(df, mol_col="Molecule"):
    df = df.copy()
    firstnames = []
    secondnames = []
    for i, row in df.iterrows():
        name = row[mol_col]
        if "|" in name:
            first, second = name.split("|")
        else:
            first, second = name, name
        firstnames.append(first)
        secondnames.append(second)

    df["First Name"] = firstnames
    df["Second Name"] = secondnames
    return df


def filter_transition_mode(df, mode):
    """Return the transitions used by the requested quantification mode."""
    if mode == "Products":
        precursor = df["Fragment Ion"].str.contains("precursor", case=False, na=False)
        # Re-imported Skyline projects can rename precursor rows to "Ion [...]".
        # Same-charge products at/above the precursor include unfragmented ions
        # and over-shifted labelled isotope rows, not genuine neutral-loss products.
        if {"Precursor Mz", "Product Mz"}.issubset(df.columns):
            precursor_mz = pd.to_numeric(df["Precursor Mz"], errors="coerce")
            product_mz = pd.to_numeric(df["Product Mz"], errors="coerce")
            unfragmented = (np.isfinite(precursor_mz) & np.isfinite(product_mz)
                            & (precursor_mz > 0) & (product_mz >= precursor_mz - 0.005))
            if {"Precursor Charge", "Product Charge"}.issubset(df.columns):
                unfragmented &= (pd.to_numeric(df["Precursor Charge"], errors="coerce")
                                 == pd.to_numeric(df["Product Charge"], errors="coerce"))
            precursor |= unfragmented
        products = df[~precursor]
        adduct_col = "Precursor Adduct" if "Precursor Adduct" in products.columns else "Adduct"
        if {"Molecule", "Molecule List Name", adduct_col}.issubset(products.columns):
            pc_h = (products["Molecule List Name"].isin(["PC", "EtherPC"])
                    & products[adduct_col].isin(["[M+H]", "[M5H2+H]"]))
            if pc_h.any():
                if "Product Mz" not in products:
                    raise ValueError("PC +H headgroup quantification requires Product Mz")
                # UltimateSPLASH PC-d5 labels the glycerol, not phosphocholine.
                # Analytes and these standards share the unshifted headgroup.
                matched = np.isclose(pd.to_numeric(products["Product Mz"], errors="coerce"),
                                     PC_HEADGROUP_MZ, rtol=0, atol=PC_HEADGROUP_TOLERANCE)
                products = products[~pc_h | matched]
        return products
    if mode == "Precursors":
        return df[df["Fragment Ion"].str.contains("precursor", case=False, na=False)]
    if mode == "Tails":
        return df[df["Fragment Ion"].str.startswith("T", na=False)]
    if mode == "Heads":
        return df[df["Fragment Ion"].str.startswith("H", na=False)]
    if mode == "All":
        return df
    raise ValueError(f"Unknown transition mode: {mode}")


def sum_transitions(df, mode="Products", drop_IS=True, normalize_IS=True, normalize_TIC=True,
                    conc_col="Concentration (uM)", customisdf="Default", is_replicate_groups=None,
                    exclusion_audit=None):
    adduct_col="Precursor Adduct" if "Precursor Adduct" in df.columns else "Adduct"
    mol_col = "Molecule"
    class_col = "Molecule List Name"
    rep_col = "Replicate Name"

    # Drop any with nan in rep_col, mol_col, adduct_col, or class_col
    df = df.dropna(subset=[rep_col, mol_col, adduct_col, class_col]).copy()

    molecules = df[mol_col].unique()
    replicates = df[rep_col].unique()

    df = clean_classes(df)

    # Fix adduct names
    # Swap [M7H2+H] to [M+H] and such
    df = df.replace({adduct_col: {"[M7H2+H]": "[M+H]", "[M7H2-H]": "[M-H]",
                                  "[M5H2+H]":"[M+H]", "[M5H2-H]":"[M-H]",
                                  "[M5H2+CH3COO]": "[M+CH3COO]", "[M5H2-CH3COO]": "[M-CH3COO]",
                                  "[M9H2+H]": "[M+H]", "[M9H2-H]": "[M-H]",
                                  "[M7H2+H-H2O]": "[M+H-H2O]", "[M5H2+NH4]": "[M+NH4]",
                                  "[M9H2+CH3COO]": "[M+CH3COO]", "[M9H2-CH3COO]": "[M-CH3COO]"}})

    rtcol = "RT" if "RT" in df.columns else "Retention Time"
    # Precompute average RT per molecule
    rt_df = df.groupby(mol_col)[rtcol].mean().reset_index()

    mzcol = "Precursor Mz"
    mz_df = df.groupby(mol_col)[mzcol].mean().reset_index()

    # Filter once, consistently for every lipid and project, before summing.
    # Retain the existing RT/mass metadata calculation above.
    df = filter_transition_mode(df, mode)

    newrows = []
    for m in molecules:
        # Get all adducts for the molecule
        subdf = df[df[mol_col]==m]
        adducts = subdf[adduct_col].unique()
        for a in adducts:
            newrow = {"Molecule": m, "Molecule List Name": df[df[mol_col] == m][class_col].values[0],
                      "MZ": mz_df[mz_df[mol_col] == m][mzcol].values[0],
                      "RT": rt_df[rt_df[mol_col] == m][rtcol].values[0],
                      "Adduct": a}
            for r in replicates:
                mask = (df[mol_col] == m) & (df[rep_col] == r) & (df[adduct_col] == a)
                subset = df[mask]

                # if m == "IS LPE(+[2]H5) 17:0":
                #     print(f"DEBUG: Subset for {m} {a} {r}:")
                #     print(subset)

                total_area = subset["Area"].sum(min_count=1)
                newrow[r] = total_area
            if normalize_IS:
                # If any IS values have a zero in any replicate, skip this IS compound
                if any(not np.isfinite(newrow[r]) or newrow[r] <= 0 for r in replicates) and m.startswith("IS"):
                    print(f"WARNING: Removing IS compound {m} with adduct {a} because it is 0")
                    if exclusion_audit is not None:
                        exclusion_audit.append({"Molecule": m, "Class": newrow[class_col], "Adduct": a,
                            "Reason": "Internal standard has zero or missing product area in at least one replicate.",
                            "Stage": "IS QC", "Rows removed": 1, "Row unit": "Lipid-adduct estimates",
                            "Replicates": "; ".join(r for r in replicates if not np.isfinite(newrow[r]) or newrow[r] <= 0)})
                    # print(newrow)
                    continue
            newrows.append(newrow)

    outdf = pd.DataFrame(newrows)

    if normalize_IS:
        print("Normalizing IS")
        outdf = normalize_is(outdf, replicates, conc_col, isdf=customisdf,
                             replicate_groups=is_replicate_groups, exclusion_audit=exclusion_audit)
        if len(outdf) == 0:
            raise ValueError("No data left after IS normalization. Check if IS compounds are present.")
        # Drop species that are 0 in all replicates after IS normalization
        outdf = outdf[~(outdf[replicates].sum(axis=1) == 0)]

    if drop_IS:
        print("Dropping IS compounds")
        # Drop any rows where Molecule starts with "IS_"
        outdf = outdf[~outdf[mol_col].str.startswith("IS")]

    if normalize_TIC:
        print("Normalizing TIC")
        outdf = normalize_tic(outdf, replicates)

    # Set any values less than minval to minval to avoid issues with log transformation later
    # for r in replicates:
    #     outdf[r] = outdf[r].apply(lambda x: max(x, minval))

    # Calculate CV for each molecule across all replicates
    # Drop all where mean is 0 to avoid division by zero
    outdf["Mean"] = outdf[replicates].mean(axis=1)
    outdf = outdf[outdf["Mean"] > 0].copy()
    outdf["CV"] = outdf[replicates].std(axis=1, ddof=1) / outdf[replicates].mean(axis=1)

    return outdf, replicates


def normalize_tic(df, replicates):
    normdf = df.copy()
    for r in replicates:
        total = normdf[r].sum()
        normdf[r] = normdf[r] / total
    return normdf

def make_trim_mask(values, trim_fraction):
    n = len(values)
    if n == 0:
        return np.zeros(0, dtype=bool)
    if trim_fraction <= 0 or n < 3:
        return np.ones(n, dtype=bool)
    ntrim = int(np.floor(n * trim_fraction))
    if ntrim == 0:
        return np.ones(n, dtype=bool)
    if 2 * ntrim >= n:
        ntrim = max((n - 1) // 2, 0)
    order = np.argsort(values)
    keep = np.zeros(n, dtype=bool)
    keep[order[ntrim:n - ntrim]] = True
    return keep

def choose_reference_sample(matrix, totals):
    scaled = matrix.divide(totals.replace(0, np.nan), axis=1)
    uq = scaled.replace(0, np.nan).quantile(0.75, axis=0)
    uq = uq.replace([np.inf, -np.inf], np.nan).dropna()
    if len(uq) > 0:
        return (uq - uq.median()).abs().idxmin()
    return (totals[totals > 0] - totals[totals > 0].median()).abs().idxmin()

def calc_tmm_factor(obs, ref, obs_total, ref_total, logratio_trim=0.30, sum_trim=0.05):
    mask = np.isfinite(obs) & np.isfinite(ref) & (obs > 0) & (ref > 0)
    if np.sum(mask) < 3:
        return 1.0

    obs = obs[mask]
    ref = ref[mask]

    norm_obs = obs / obs_total
    norm_ref = ref / ref_total
    m_vals = np.log2(norm_obs / norm_ref)
    a_vals = 0.5 * np.log2(norm_obs * norm_ref)

    keep = make_trim_mask(m_vals, logratio_trim) & make_trim_mask(a_vals, sum_trim)
    if np.sum(keep) == 0:
        keep = np.ones(len(m_vals), dtype=bool)

    obs_keep = obs[keep]
    ref_keep = ref[keep]
    m_keep = m_vals[keep]

    with np.errstate(divide='ignore', invalid='ignore'):
        inv_var = ((obs_total - obs_keep) / (obs_total * obs_keep)) + ((ref_total - ref_keep) / (ref_total * ref_keep))
        weights = np.where(np.isfinite(inv_var) & (inv_var > 0), 1.0 / inv_var, np.nan)

    valid = np.isfinite(weights) & np.isfinite(m_keep)
    if np.sum(valid) == 0:
        mean_m = np.mean(m_keep)
    else:
        mean_m = np.average(m_keep[valid], weights=weights[valid])

    factor = float(2 ** mean_m)
    if not np.isfinite(factor) or factor <= 0:
        return 1.0
    return factor


def normalize_trimmed_mean_mvalues(df, set1, set2):
    """
    Normalize replicate columns using trimmed mean of M-values (TMM).

    This follows the standard TMM idea of estimating a per-sample scaling factor
    against a reference sample by trimming extreme M-values and A-values. The
    replicate columns from set1 and set2 are scaled while metadata columns are
    preserved unchanged.

    :param df: Wide dataframe containing replicate intensity columns.
    :param set1: Replicate column names for group 1.
    :param set2: Replicate column names for group 2.
    :return: Copy of df with normalized replicate columns.
    """
    replicates = list(dict.fromkeys(list(set1) + list(set2)))
    normdf = df.copy()

    data = normdf[replicates].apply(pd.to_numeric, errors="coerce").fillna(0.0).astype(float)
    lib_sizes = data.sum(axis=0).astype(float)

    ref_col = choose_reference_sample(data, lib_sizes)
    ref_total = float(lib_sizes[ref_col])
    ref_values = data[ref_col].values.astype(float)

    factors = pd.Series(1.0, index=replicates, dtype=float)
    for r in replicates:
        obs_total = float(lib_sizes[r])
        if obs_total <= 0:
            print(f"WARNING: Replicate '{r}' has zero total intensity. Leaving as zeros after TMM normalization.")
            factors[r] = 1.0
            continue
        if r == ref_col:
            factors[r] = 1.0
            continue
        factors[r] = calc_tmm_factor(data[r].values.astype(float), ref_values, obs_total, ref_total)

    valid_factors = factors[(factors > 0) & np.isfinite(factors)]
    if len(valid_factors) > 0:
        factors = factors / np.exp(np.mean(np.log(valid_factors)))

    effective_lib_sizes = lib_sizes * factors
    valid_effective = effective_lib_sizes[(effective_lib_sizes > 0) & np.isfinite(effective_lib_sizes)]
    if len(valid_effective) == 0:
        print("WARNING: No valid effective library sizes were produced by TMM normalization. Returning unmodified data.")
        normdf[replicates] = data
        return normdf

    scale_target = float(np.exp(np.mean(np.log(valid_effective))))
    normalized = data.copy()
    for r in replicates:
        denom = effective_lib_sizes[r]
        if not np.isfinite(denom) or denom <= 0:
            normalized[r] = 0.0
        else:
            normalized[r] = data[r] / denom * scale_target

    normdf[replicates] = normalized
    return normdf


def get_is_conc_from_isdf(isrow, isdf, conc_col):
    try:
        isname = isrow["Molecule"].values[0].strip()
        isclass = isrow["Molecule List Name"].values[0]
    except Exception as e:
        isname = isrow["Molecule"].strip()
        isclass = isrow["Molecule List Name"].strip()

    # Remove only leading 'IS ' (with space), keep rest intact
    if isname.startswith("IS "):
        isname = isname[3:].strip()
    # Skyline and the standard table may use '-' and '_' interchangeably in lipid tails.
    normalized_isname = isname.replace("-", "_")
    normalized_isdf_names = (isdf["Lipid"].astype(str).str.strip()
                             .str.replace(r"^IS\s+", "", regex=True).str.replace("-", "_", regex=False))
    # Find matches
    isdf_row = isdf[(normalized_isdf_names == normalized_isname) & (isdf["Class"] == isclass)]
    # If none or multiple matches, print warning and return None
    if len(isdf_row) == 0:
        print(f"WARNING: No matching IS concentration found for '{isname}' of class '{isclass}'. Quantification unavailable.")
        return None
    if len(isdf_row) > 1:
        print(f"WARNING: Multiple matching IS found in IS dataframe for {isname}")
        return None
    conc = isdf_row[conc_col].values[0]
    return conc


def get_area_matched_is(row, isrows, replicates):
    """Select one internal standard whose mean peak area is closest to the analyte."""
    analyte_values = pd.to_numeric(row[replicates], errors="coerce")
    analyte_area = analyte_values.mean()
    if not np.isfinite(analyte_area) or analyte_area <= 0:
        return None

    standard_values = isrows[replicates].apply(pd.to_numeric, errors="coerce")
    standard_areas = standard_values.mean(axis=1)
    valid = np.isfinite(standard_areas) & (standard_areas > 0)
    if not valid.any():
        return None

    # Compare fold differences so standards above and below the analyte are treated symmetrically.
    area_distances = np.abs(np.log(standard_areas[valid] / analyte_area))
    return isrows.loc[[area_distances.idxmin()]]


def calc_is(df, isrows, replicates, conc_col, isdf=default_isdf, match_by_area=False):
    # Copy just the replicate columns
    normdata = df[replicates].copy()
    cls = df["Molecule List Name"].values[0] if "Molecule List Name" in df.columns else "Unknown"
    adduct = df["Adduct"].values[0] if "Adduct" in df.columns else "Unknown"
    if len(df) == 0:
        print("WARNING: No data rows found. Skipping IS calculation.")
        return normdata
    if len(isrows) == 0:
        print(f"WARNING: No IS rows found for class {cls} and adduct {adduct}. Excluding this adduct from quantification.")
        return normdata * np.nan

    # Single IS row
    if len(isrows) == 1:
        for r in replicates:
            isval = float(isrows[r].values[0])
            if isdf is not None:
                conc = get_is_conc_from_isdf(isrows, isdf, conc_col)
                if conc is not None:
                    isval /= conc
                    # print(f"INFO: Normalizing by IS value {isval} for class {cls} and adduct {adduct} in replicate {r}.")
                else:
                    print(f"WARNING: No valid concentration found for IS row in class {cls} and adduct {adduct}. Excluding this adduct from quantification.")
                    normdata[r] = np.nan
                    continue
            if isval == 0:
                print(f"WARNING: IS value is zero for class {cls} and adduct {adduct} in replicate {r}. Excluding this adduct from quantification.")
                normdata[r] = np.nan
                continue
            normdata[r] = normdata[r].astype(float) / isval

    # Area matching uses one standard per analyte across all replicates so that the
    # selected standard does not change between experimental groups. Other classes match by RT.
    else:
        for i, row in df.iterrows():
            rt = row.get("RT", None)
            area_matched_isrow = get_area_matched_is(row, isrows, replicates) if match_by_area else None
            if area_matched_isrow is None and rt is None:
                print(f"WARNING: No RT found for row {i} but multiple standards are found. Excluding this adduct from quantification.")
                normdata.loc[i, replicates] = np.nan
                continue

            for r in replicates:
                possible_rows = isrows[isrows[r] > 0]
                if area_matched_isrow is not None and float(area_matched_isrow[r].iloc[0]) > 0:
                    closest_isrow = area_matched_isrow
                else:
                    closest_isrow = possible_rows.iloc[(possible_rows["RT"] - rt).abs().argsort()[:1]]
                if len(closest_isrow) == 0:
                    print(f"WARNING: No IS rows with nonzero value found for class {cls} and adduct {adduct} in replicate {r}. Excluding this adduct from quantification.")
                    normdata.at[i, r] = np.nan
                    continue
                isval = float(closest_isrow[r].values[0])
                if isdf is not None:
                    conc = get_is_conc_from_isdf(closest_isrow, isdf, conc_col)
                    if conc is not None:
                        rawarea = normdata.at[i, r]
                        intratio = rawarea / isval
                        if intratio > 50:
                            if match_by_area:
                                standard_name = closest_isrow["Molecule"].iloc[0]
                                print(f"WARNING: Analyte/internal standard peak-area ratio of {intratio} for "
                                      f"{row['Molecule']} ({adduct}, {r}) exceeds 50. Using {standard_name}, "
                                      f"the closest-area standard; retaining concentration-corrected quantification. "
                                      f"Review the assay's validated range.")
                            else:
                                print(f"WARNING: High internal standard ratio of {intratio} for row {i} in class {cls} and adduct {adduct} in replicate {r}. This may indicate an issue with the IS normalization.")
                        isval /= conc
                        #
                        # corrected = float(rawarea) / isval
                        # if corrected > 1:
                        #     print(row.to_string())
                        #     print(closest_isrow.to_string())
                        #     print("ISval:", isval, "Raw area:", rawarea, "Corrected area:", corrected, "Int Ratio:", intratio)
                    else:
                        print(f"WARNING: No valid concentration found for IS row in class {cls} and adduct {adduct}. Excluding this adduct from quantification.")
                        normdata.at[i, r] = np.nan
                        continue
                if isval == 0:
                    print(f"WARNING: IS value is zero for class {cls} and adduct {adduct} in replicate {r}. Excluding this adduct from quantification.")
                    normdata.at[i, r] = np.nan
                    continue
                normdata.at[i, r] = float(normdata.at[i, r]) / isval

    return normdata

def get_is(df, class_name, adduct, is_mapper_dict=None):
    if is_mapper_dict is not None and class_name in is_mapper_dict:
        mapped_class = is_mapper_dict[class_name]
        print(f"Mapping class '{class_name}' to '{mapped_class}' for IS lookup.")
        class_name = mapped_class
    if "Ether" in class_name:
        # print(f"Class '{class_name}' contains 'Ether'. Attempting to find IS for 'non-Ether' version of the class.")
        class_name = class_name.replace("Ether", "").strip()
    subdf = df[(df["Molecule List Name"] == class_name) & (df["Adduct"] == adduct)]
    isrows = subdf[subdf["Molecule"].str.startswith("IS ")]
    return isrows

def dereplicate_adducts(normdf, replicates):
    # Average valid adduct estimates for each lipid without allowing a failed adduct to dilute the result.
    outrows = []
    lipids = normdf["Molecule"].unique()
    for lipid in lipids:
        subdf = normdf[normdf["Molecule"] == lipid]
        adducts = subdf["Adduct"].unique()
        newrow = subdf.iloc[0].copy()
        if len(adducts) > 1:
            newrow["Adduct"] = "+".join(adducts)
        preferred = preferred_quant_adducts.get(newrow["Molecule List Name"])
        used_adducts = []
        for r in replicates:
            if preferred is not None:
                # Each adduct was already normalized to its own matching IS.
                # Choose acetate per replicate when usable; otherwise use +H.
                newrow[r] = np.nan
                newrow[f"Quant Adduct {r}"] = ""
                for adduct in [preferred, "[M+H]"]:
                    values = pd.to_numeric(subdf.loc[subdf["Adduct"] == adduct, r], errors="coerce")
                    values = values[np.isfinite(values) & (values > 0)]
                    if len(values):
                        newrow[r] = float(values.iloc[0])
                        newrow[f"Quant Adduct {r}"] = adduct
                        used_adducts.append(adduct)
                        break
                if not np.isfinite(newrow[r]):
                    # Retain a measured zero only if neither adduct is usable.
                    values = pd.to_numeric(subdf[r], errors="coerce")
                    if np.isfinite(values).all() and values.eq(0).all():
                        newrow[r] = 0.
                    else:
                        print(f"WARNING: No valid acetate or +H quantification for {lipid} "
                              f"in {r}; retaining a missing estimate, not zero.")
                continue
            valid_values = pd.to_numeric(subdf[r], errors="coerce").dropna()
            newrow[r] = valid_values.mean() if len(valid_values) > 0 else np.nan
        if preferred is not None:
            newrow["Adduct"] = "+".join(dict.fromkeys(used_adducts))
        outrows.append(newrow)
    outdf = pd.DataFrame(outrows)
    print(len(normdf), "rows before dereplication, ", len(outdf), "rows after dereplication.")
    return outdf


def check_standard_variation(df, replicate_groups, exclusion_audit=None):
    """Flag within-standard area variation without rejecting valid standards.

    Different concentrations between standard species are expected and handled
    by calc_is. Nonpositive/nonfinite denominators remain unusable.
    """
    if replicate_groups is None:
        return df
    keep = pd.Series(True, index=df.index)
    for index, row in df.loc[df["Molecule"].str.startswith("IS ")].iterrows():
        for group in replicate_groups:
            if len(group) < 2:
                continue
            values = pd.to_numeric(row[group], errors="coerce").to_numpy(dtype=float)
            invalid = not np.isfinite(values).all() or (values <= 0).any()
            ratio = np.inf if invalid else values.max() / values.min()
            if invalid:
                keep.at[index] = False
                reason = "Internal standard has nonpositive/nonfinite area within comparison group."
                print(f"WARNING: Removing IS {row['Molecule']} ({row['Adduct']}): {reason}")
            elif ratio > IS_AREA_VARIATION_WARNING_FOLD:
                reason = (f"IS area varies {ratio:.2f}-fold within comparison group; "
                          "retained for concentration-corrected normalization. Review peaks and standard loading.")
                print(f"WARNING: {row['Molecule']} ({row['Adduct']}): {reason}")
            else:
                continue
            if exclusion_audit is not None:
                exclusion_audit.append({"Molecule": row["Molecule"], "Class": row["Molecule List Name"],
                        "Adduct": row["Adduct"], "Reason": reason,
                        "Stage": "IS QC" if invalid else "IS QC warning", "Rows removed": int(invalid),
                        "Row unit": "Lipid-adduct estimates",
                        "Replicates": "; ".join(group)})
            if invalid:
                break
    return df.loc[keep].copy()


def normalize_is(df, replicates, conc_col, isdf="Default", dereplicate=True, replicate_groups=None,
                 exclusion_audit=None):
    if type(isdf) == str:
        if isdf == "Default":
            isdf = default_isdf
        else:
            raise ValueError(f"Invalid value for isdf: {isdf}. Must be 'Default' or a custom IS dataframe.")

    normdf = check_standard_variation(df.copy(), replicate_groups, exclusion_audit=exclusion_audit)
    for class_name, preferred_adduct in preferred_quant_adducts.items():
        class_mask = normdf["Molecule List Name"] == class_name
        excluded_mask = class_mask & ~normdf["Adduct"].isin([preferred_adduct, "[M+H]"])
        excluded_count = excluded_mask.sum()
        if excluded_count > 0:
            print(f"Using {preferred_adduct} with +H fallback for {class_name} quantification; "
                  f"excluding {excluded_count} rows with other adducts.")
            normdf = normdf[~excluded_mask].copy()

    classes = normdf["Molecule List Name"].unique()
    for r in replicates:
        normdf[r] = normdf[r].astype(float)

    for c in classes:
        mask1 = normdf["Molecule List Name"] == c
        subdf_class = normdf[mask1]
        adducts = subdf_class["Adduct"].unique()
        for a in adducts:
            mask2 = subdf_class["Adduct"] == a
            subdf = subdf_class[mask2]
            # Remove ones with IS in the name from subdf to avoid self-normalization
            mask3 = ~subdf["Molecule"].str.startswith("IS ")
            subdf = subdf[mask3]
            if len(subdf)==0:
                continue
            # Get IS rows for this class
            isrow = get_is(normdf, c, a, is_mapper_dict=default_is_mapper)
            # Calculate normalized data for this class and replicate
            normdata = calc_is(subdf, isrow, replicates, conc_col, isdf=isdf,
                               match_by_area=c in area_match_is_classes)
            # Merge back into normdf
            mask = (mask1) & (mask2) & (mask3)
            normdf.loc[mask, replicates] = normdata

    # Repeat with IS lipids
    for c in classes:
        mask1 = normdf["Molecule List Name"] == c
        subdf_class = normdf[mask1]
        adducts = subdf_class["Adduct"].unique()
        for a in adducts:
            mask2 = subdf_class["Adduct"] == a
            subdf = subdf_class[mask2]
            mask3 = subdf["Molecule"].str.startswith("IS ")
            subdf = subdf[mask3]
            # Get IS rows for this class
            isrow = get_is(normdf, c, a, is_mapper_dict=default_is_mapper)
            if len(isrow) == 0 or len(subdf)==0:
                continue
            # Calculate normalized data for this class and replicate
            normdata = calc_is(subdf, isrow, replicates, conc_col, isdf=isdf)
            # Merge back into normdf
            mask = (mask1) & (mask2) & (mask3)
            normdf.loc[mask, replicates] = normdata

    if dereplicate:
        normdf = dereplicate_adducts(normdf, replicates)

    return normdf


def detection_table(quantified, source, replicates, drop_IS=True):
    """Restore the input catalogue and distinguish missing rows from failed quantification.

    source is captured BEFORE transition filtering. An input row that fails
    quantification must never become a zero measurement.
    """
    source = source.copy()
    clean_classes(source)
    keys = ["Molecule", "Molecule List Name"]
    catalogue = source[keys].drop_duplicates()
    if drop_IS:
        catalogue = catalogue.loc[~catalogue["Molecule"].str.startswith("IS")]
    out = catalogue.merge(quantified, on=keys, how="left", validate="one_to_one")
    for r in replicates:
        present = pd.MultiIndex.from_frame(source.loc[source["Replicate Name"].eq(r), keys].drop_duplicates())
        has_rows = pd.MultiIndex.from_frame(out[keys]).isin(present)
        out.loc[~has_rows, r] = np.nan
        out[f"Detection {r}"] = np.where(
            ~has_rows, "Not measured",
            np.where(np.isfinite(out[r]) & out[r].gt(0), "Quantified",
            np.where(out[r].eq(0), "Not detected", "Quantification failed")))
    return out


def comparison_stats(df, set1, set2, paired=False, bh_correction=True, include_nondetections=False):
    """Keep every lipid and optionally test complete measurements containing zeros.

    Tests including non-detections assume zero signal, not known biological
    absence. Preserve infinite log2 fold changes rather than adding a pseudocount.
    """
    out = df.copy()
    replicates = set1 + set2
    out["Mean"] = out[replicates].mean(axis=1)
    for label, group in [("Set1", set1), ("Set2", set2)]:
        values = out[group]
        out[f"{label}_Detected"] = (np.isfinite(values) & values.gt(0)).sum(axis=1)
        out[f"{label}_Available"] = np.isfinite(values).sum(axis=1)
        complete = np.isfinite(values).all(axis=1)
        out[f"{label}_Mean"] = values.mean(axis=1).where(complete)
        out[f"{label}_CV"] = (values.std(axis=1) / values.mean(axis=1)).where(complete)
    d1, d2 = out["Set1_Detected"], out["Set2_Detected"]
    a1, a2 = out["Set1_Available"], out["Set2_Available"]
    out["Detection change"] = np.select(
        [(d1 > 0) & (d2 == 0) & (a2 == len(set2)),
         (d2 > 0) & (d1 == 0) & (a1 == len(set1)),
         (d1 > 0) & (d2 == 0), (d2 > 0) & (d1 == 0), (d1 > 0) & (d2 > 0)],
        ["Set1 only; Set2 ND", "Set2 only; Set1 ND", "Set1 detected; Set2 unavailable",
         "Set2 detected; Set1 unavailable", "Detected in both"], default="No quantified detection")
    eligible = np.isfinite(out[replicates]).all(axis=1) & out[replicates].gt(0).all(axis=1)
    if include_nondetections:
        eligible = (np.isfinite(out[replicates]).all(axis=1)
                    & out[replicates].ge(0).all(axis=1) & out[replicates].gt(0).any(axis=1))
        if paired:
            # Paired ratios require a positive denominator in every control replicate.
            eligible &= out[set1].gt(0).all(axis=1)
    if min(len(set1), len(set2)) < 2:
        eligible[:] = False
    out["Comparison status"] = "Not tested: ND or unavailable measurements"
    out["Test basis"] = np.where(out[replicates].eq(0).any(axis=1),
                                "Includes zero non-detections", "Positive measurements")
    metrics = ["FC", "Log2FC", "StdDevFC", "CI95_FC", "p-value", "Raw p-value", "Log10p"]
    summarycols = [f"Ratio_{i + 1}" for i in range(len(set1))] if paired else ["Set1_Mean", "Set2_Mean"]
    for col in metrics + ([c for c in summarycols] if paired else []):
        out[col] = np.nan
    if eligible.any():
        tested, _ = stats_calcs(out.loc[eligible].copy(), set1, set2, paired, bh_correction,
                               allow_zeros=include_nondetections)
        out.loc[eligible, "Comparison status"] = "Not tested: undefined t-test"
        out.loc[tested.index, "Comparison status"] = "Tested"
        for col in metrics + (summarycols if paired else []):
            out.loc[tested.index, col] = tested[col]
    return out, summarycols


def stats_calcs(df, replicates_set1, replicates_set2, paired=False, bh_correction=True, allow_zeros=False):

    df["Mean"] = df[replicates_set1 + replicates_set2].mean(axis=1)
    set1_means = df[replicates_set1].mean(axis=1)
    set2_means = df[replicates_set2].mean(axis=1)
    df["Set1_Mean"] = set1_means
    df["Set2_Mean"] = set2_means
    df["Set1_CV"] = df[replicates_set1].std(axis=1, ddof=1) / set1_means
    df["Set2_CV"] = df[replicates_set2].std(axis=1, ddof=1) / set2_means

    summarycols = []
    if paired:
        if len(replicates_set1) != len(replicates_set2):
            raise ValueError("For paired t-test, the number of replicates in both sets must be equal.")
        set1_values = np.array(df[replicates_set1], dtype=np.float64)
        set2_values = np.array(df[replicates_set2], dtype=np.float64)
        ratios = set2_values / set1_values
        # Put ratios in their own columns
        for i in range(ratios.shape[1]):
            df[f"Ratio_{i + 1}"] = ratios[:, i]
            summarycols.append(f"Ratio_{i + 1}")
        ratiomean = ratios.mean(axis=1)
        df["FC"] = ratiomean
        with np.errstate(divide="ignore", invalid="ignore"):
            df["Log2FC"] = np.log2(np.array(ratiomean))
        df["StdDevFC"] = np.std(ratios, axis=1, ddof=1)
        df["CI95_FC"] = stats.t.ppf(0.975, df=len(replicates_set1) - 1) * (df["StdDevFC"] / np.sqrt(len(replicates_set1)))

    else:
        summarycols.append("Set1_Mean")
        summarycols.append("Set2_Mean")
        with np.errstate(divide="ignore", invalid="ignore"):
            df["FC"] = set2_means / set1_means
            df["Log2FC"] = np.log2(set2_means) - np.log2(set1_means)
        # Propagate standard deviation
        set1_std = df[replicates_set1].std(axis=1, ddof=1)
        set2_std = df[replicates_set2].std(axis=1, ddof=1)
        # Make sure to handle cases where set1_means or set2_means are zero to avoid division by zero
        df["StdDevFC"] = df["FC"] * np.sqrt(safedivide(set1_std , set1_means) ** 2 + safedivide(set2_std , set2_means) ** 2)
        df["CI95_FC"] = stats.t.ppf(0.975, df=len(replicates_set1) + len(replicates_set2) - 2) * (
                    df["StdDevFC"] / np.sqrt(len(replicates_set1) + len(replicates_set2)))

    # Drop any with a fold change of 0
    if not allow_zeros:
        df = df[df["FC"] != 0].copy()

    p_values = []
    for index, row in df.iterrows():
        set1_values = row[replicates_set1].values
        set2_values = row[replicates_set2].values

        set1_values = np.array(set1_values, dtype=np.float64)
        set2_values = np.array(set2_values, dtype=np.float64)

        if index == 0:
            print("N Set1:", len(set1_values), "N Set2:", len(set2_values))

        if paired:
            _, p = stats.ttest_rel(set1_values, set2_values)
        else:
            _, p = stats.ttest_ind(set1_values, set2_values)
        p_values.append(p)
    df["p-value"] = p_values

    # Drop any that are inf or nan in Log2FC or Log10p
    if allow_zeros:
        df = df.loc[df["FC"].notna() & np.isfinite(df["p-value"])].copy()
        # Ratio uncertainty is undefined when a group mean is zero.
        df.loc[df["Set1_Mean"].eq(0) | df["Set2_Mean"].eq(0), ["StdDevFC", "CI95_FC"]] = np.nan
    else:
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=["FC", "p-value"])
    df["Raw p-value"] = df["p-value"]

    if bh_correction and len(df):
        corrected = multipletests(df["p-value"], method='fdr_bh')
        df["p-value"] = corrected[1]

    df["Log10p"] = -np.log10(df["p-value"] + 1e-10)

    return df, summarycols


def class_calcs(df, replicate_set1, replicate_set2, class_col="Molecule List Name", normalize=True, paired=True,
                bh_correction=True, retain_unavailable=False, include_nondetections=False):
    classes = df[class_col].unique()
    replicates = replicate_set1 + replicate_set2

    results = []
    for cls in classes:
        newrow = {class_col: cls}
        for r in replicates:
            mask = df[class_col] == cls
            class_sum = df[mask][r].sum(min_count=1)
            newrow[r] = class_sum
        results.append(newrow)
    outdf = pd.DataFrame(results)

    if normalize:
        outdf = normalize_tic(outdf, replicates)

    calculate = comparison_stats if retain_unavailable else stats_calcs
    extra = {"include_nondetections": include_nondetections} if retain_unavailable else {}
    outdf, summarycols = calculate(outdf, replicate_set1, replicate_set2, paired=paired,
                                  bh_correction=bh_correction, **extra)
    # print(outdf)
    return outdf, summarycols


def make_volcano_plot(df, log2fc_col="Log2FC", neglogp_col="Log10p", title="Volcano Plot", ax=None, xlims=[-8, 8],
                      fold_range=1., show_detection_summary=True, show_group_label=True):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)
    # Color by Molecule List Name
    ax = plt.gca()
    if show_detection_summary and "Detection change" in df:
        changes = df["Detection change"].value_counts()
        notes = [f"{name}: {count}" for name, count in changes.items() if name != "Detected in both"]
        if notes:
            ax.text(0.5, -0.24, "Detection changes\n" + "\n".join(notes),
                    transform=ax.transAxes, va="top", ha="center", fontsize=8)
    group_label = ""
    if show_group_label and len(df) and "Set1 name" in df:
        group_label = f"\nSet2 / Set1: {df['Set2 name'].iloc[0]} / {df['Set1 name'].iloc[0]}"
    # Infinite FC values retain their statistical y-coordinate and get distinct
    # boundary markers. Clipping changes display positions only.
    df = df.loc[df[log2fc_col].notna() & np.isfinite(df[neglogp_col])].copy()
    classes = df["Molecule List Name"].unique()
    # colors = plt.cm.get_cmap('tab10', len(classes))
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    # sizes = df["Mean"]
    # # Scale sizes to between 20 and 50
    # sizes = 20 + (sizes - sizes.min()) / (sizes.max() - sizes.min()) * (50 - 20)

    # Clip log2fc values to xlims for better visualization
    plot_x = df[log2fc_col].clip(lower=xlims[0]+0.15, upper=xlims[1]-0.15)
    finite = np.isfinite(df[log2fc_col])
    plt.scatter(plot_x[finite], df.loc[finite, neglogp_col], alpha=0.7,
                c=df.loc[finite, "Molecule List Name"].map(get_color_for_class))
    for mask, marker, label in [(np.isneginf(df[log2fc_col]), "<", "Set2 ND (log2 FC = -infinity)"),
                                (np.isposinf(df[log2fc_col]), ">", "Set1 ND (log2 FC = +infinity)")]:
        if mask.any():
            plt.scatter(plot_x[mask], df.loc[mask, neglogp_col], marker=marker, s=45,
                        c=df.loc[mask, "Molecule List Name"].map(get_color_for_class),
                        edgecolors="black", linewidths=0.4, alpha=0.7, label=label)
    if (~finite).any():
        plt.legend(loc="best", fontsize=8, frameon=False)
    plt.title(title)
    plt.xlabel("Log\u2082 Fold Change" + group_label, fontsize=9)
    ylabel = "-Log\u2081\u2080 adjusted p-value" if df.attrs.get("bh_correction", False) else "-Log\u2081\u2080 p-value"
    plt.ylabel(ylabel)
    plt.axhline(y=-np.log10(0.05), color='r', linestyle='--')
    plt.axvline(x=fold_range, color='g', linestyle='--')
    plt.axvline(x=-fold_range, color='g', linestyle='--')
    plt.xlim(xlims)

    # Add cursor hover to show molecule names
    annot = plt.annotate("", xy=(0,0), xytext=(20,20),textcoords="offset points",
                         bbox=dict(boxstyle="round", fc="w"),
                         arrowprops=dict(arrowstyle="->"))
    annot.set_visible(False)
    def update_annot(ind):
        pos = sc.get_offsets()[ind["ind"][0]]
        annot.xy = pos
        text = "\n".join([df["Molecule"].iloc[n] for n in ind["ind"]])
        annot.set_text(text)
        annot.get_bbox_patch().set_facecolor("lightyellow")
        annot.get_bbox_patch().set_alpha(0.9)
    def hover(event):
        vis = annot.get_visible()
        if event.inaxes == ax:
            cont, ind = sc.contains(event)
            if cont:
                update_annot(ind)
                annot.set_visible(True)
                plt.draw()
            else:
                if vis:
                    annot.set_visible(False)
                    plt.draw()

    sc = plt.scatter(plot_x, df[neglogp_col], alpha=0)  # Invisible scatter for hover
    plt.gcf().canvas.mpl_connect("motion_notify_event", hover)




def make_class_plots(classdf, ax=None, fclimit=None):
    if ax is None:
        plt.figure(figsize=(10, 6))
    else:
        plt.sca(ax)
    classdf = classdf.loc[np.isfinite(classdf["FC"]) & np.isfinite(classdf["CI95_FC"])]
    if classdf.empty:
        plt.text(0.5, 0.5, "No testable class fold changes", ha="center", transform=plt.gca().transAxes)
        plt.title("Class Fold Changes")
        return
    # classes = classdf["Molecule List Name"].unique()
    # class_color_map = {cls: colors(i) for i, cls in enumerate(classes)}
    colors = [get_color_for_class(cls) for cls in classdf["Molecule List Name"]]

    # Labels are class name plus p-value
    labels = [f"{cls}\n(p={pval:.2f})" for cls, pval in zip(classdf["Molecule List Name"], classdf["p-value"])]

    if len(labels) > 10:
        rotation=90
    else:
        rotation=0

    plt.bar(classdf["Molecule List Name"], classdf["FC"], yerr=classdf["CI95_FC"], color=colors,
            label=labels, capsize=5)
    plt.xticks(rotation=rotation)
    # line at y=1
    plt.axhline(y=1, color='r', linestyle='--')
    plt.ylabel("Fold Change")
    plt.title("Class Fold Changes")
    plt.ylim(0, fclimit)
    if len(labels) <= 10:
        plt.legend()


def pca_plot(df, replicates_set1, replicates_set2, title="PCA Plot", ax=None):
    replicates = replicates_set1 + replicates_set2
    df = df.loc[np.isfinite(df[replicates]).all(axis=1)]
    if len(df) < 2:
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "PCA unavailable: fewer than two complete lipids", ha="center", transform=ax.transAxes)
        return
    data = df[replicates].values
    data = StandardScaler().fit_transform(data.T)

    pca = PCA(n_components=2)
    principalComponents = pca.fit_transform(data)
    pc_df = pd.DataFrame(data=principalComponents, columns=['PC1', 'PC2'])
    pc_df['Group'] = ['Set1'] * len(replicates_set1) + ['Set2'] * len(replicates_set2)

    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    colors = {'Set1': 'b', 'Set2': 'r'}
    for group in pc_df['Group'].unique():
        indicesToKeep = pc_df['Group'] == group
        plt.scatter(pc_df.loc[indicesToKeep, 'PC1'], pc_df.loc[indicesToKeep, 'PC2'], c=colors[group], s=50)
    plt.title(title)
    plt.xlabel(f'PC1 - {pca.explained_variance_ratio_[0] * 100:.2f}%')
    plt.ylabel(f'PC2 - {pca.explained_variance_ratio_[1] * 100:.2f}%')
    plt.legend(pc_df['Group'].unique())

def heatmap_plot(df, replicates, title="Heatmap", ax=None, lindex=0, colormap='viridis', normtype="minmax",
                 scale="linear", class_rect=False, fontsize=5, color_classes=False, colorbar=None):
    data = df[replicates].to_numpy(dtype=float, copy=True)
    data[~np.isfinite(data)] = np.nan
    if not np.isfinite(data).any():
        if ax is None:
            _, ax = plt.subplots()
        ax.text(0.5, 0.5, "No defined values", ha="center", transform=ax.transAxes)
        ax.set_title(title)
        return
    vmax = 1-1e-12
    vmin = 1e-12
    dmax = np.nanmax(data)
    dmin = np.nanmin(data)
    if normtype=="minmax":
        # Normalize data to 0-1 for better color mapping
        data = (data - dmin) / (dmax - dmin if dmax != dmin else 1.)
    elif normtype=="zscore":
        data = StandardScaler().fit_transform(data)
    elif normtype=="zerocentered":
        data = (data / (2 * (np.nanmax(np.abs(data)) or 1.))) + 0.5
        extreme = np.amax(np.abs([dmin, dmax]))
        dmin = -extreme
        dmax = extreme
    elif normtype=="onecentered":
        data = ((data - 1) / (2 * (np.nanmax(np.abs(data - 1)) or 1.))) + 0.5
        extreme = np.amax(np.abs([dmin, dmax]))
        dmin = -extreme
        dmax = extreme
    elif normtype=="zerocentered_clip":
        data = data + 0.5
        data[data < vmin] = vmin
        data[data > vmax] = vmax
    elif normtype=="onecentered_clip":
        data = data - 0.5
        data[data < vmin] = vmin
        data[data > vmax] = vmax

    print("Min max of heatmap data:", np.min(data), np.max(data))
    if ax is None:
        plt.figure(figsize=(10, 8))
    else:
        plt.sca(ax)

    cmap = plt.get_cmap(colormap).copy()
    cmap.set_bad("lightgray")
    plt.imshow(np.ma.masked_invalid(data), aspect='auto', cmap=cmap, norm=scale, vmax=vmax, vmin=vmin)
    # Labels are last 3 characters of replicates
    labels = [r[lindex:] for r in replicates]
    plt.xticks(ticks=np.arange(len(replicates)), labels=labels,
               rotation=30 if len(replicates) > 4 else 0, fontsize=9)
    # If color classes, color ytick labels by class with a colored dot ahead of the label
    ylabels = df['Molecule']
    if color_classes:
        classes = df["Molecule List Name"].values
        for i, cls in enumerate(classes):
            plt.text(-0.6, i, "●", color=get_color_for_class(cls), fontsize=fontsize*1.5, va='center', ha='right')
        # Add 3 spaces to the end of each label to avoid overlap with the colored dot
        ylabels = [label + "     " for label in ylabels]

    # Retain every data row while avoiding hundreds of overlapping tick labels.
    tickrows = np.arange(0, len(df), max(1, int(np.ceil(len(df) / 25.))))
    plt.yticks(ticks=tickrows, labels=np.asarray(ylabels)[tickrows], fontsize=fontsize)
    def describe_cell(x, y):
        row, col = int(round(y)), int(round(x))
        if 0 <= row < len(df) and 0 <= col < len(replicates):
            return f"{df['Molecule'].iloc[row]} | {replicates[col]}: {df[replicates[col]].iloc[row]}"
        return ""
    plt.gca().format_coord = describe_cell

    if colorbar is not None:
        cb = plt.colorbar(label=colorbar, pad=0.02, aspect=75)
        # Label ticks on color bar with dmin, 0, dmax
        cb.set_ticks([vmin, 0.5, vmax])
        cb.set_ticklabels([f"{dmin:.1f}", "0", f"{dmax:.1f}"])

    if class_rect:
        # Draw colored rectangles around each class
        classes = df["Molecule List Name"].values
        class_changes = [0]
        for i in range(1, len(classes)):
            if classes[i] != classes[i - 1]:
                class_changes.append(i)
        class_changes.append(len(classes))
        for i in range(len(class_changes) - 1):
            start = class_changes[i] - 0.5
            end = class_changes[i + 1] - 0.5
            rect = mpl.patches.Rectangle(( -0.5, start), len(replicates), end - start,
                                         linewidth=1, edgecolor=get_color_for_class(classes[class_changes[i]]),
                                         facecolor='none')
            plt.gca().add_patch(rect)


    plt.title(title + "\nGray = unavailable / undefined")

def pie_chart(df, set, ax=None, title="Class Distribution", otherthresh=0.035):
    class_sums = df.groupby("Molecule List Name")[set].sum()
    if ax is None:
        plt.figure(figsize=(8, 8))
    else:
        plt.sca(ax)

    # For sums that are less than 1% of total, group into "Other"
    total = class_sums.sum()
    if not np.isfinite(total) or total <= 0:
        plt.text(0.5, 0.5, "No quantified signal", ha="center", transform=plt.gca().transAxes)
        plt.title(title)
        return
    b1 = class_sums / total >= otherthresh
    class_sums_filtered = class_sums[b1]
    colors = [get_color_for_class(cls) for cls in class_sums_filtered.index]
    other_sum = total - class_sums_filtered.sum()
    if other_sum > 0:
        class_sums_filtered["Other"] = other_sum
        colors.append("#CCCCCC")


    plt.pie(class_sums_filtered, labels=class_sums_filtered.index, colors=colors, autopct='%1.0f%%', startangle=140)
    plt.title(title)

def plot_cv_hist(normdf, ax=None, nbins=50):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)
    classes = normdf["Molecule List Name"].unique()
    bins = np.linspace(0, 1, nbins)
    # Add maxcv to the end of bins to capture any CVs above 1
    bins = np.append(bins, normdf["CV"].max())
    hists = []
    for cls in classes:
        subset = normdf[normdf["Molecule List Name"] == cls]
        histdata = subset["CV"].dropna()
        hist, _ = np.histogram(histdata, bins=bins)
        hists.append(hist)
    hists = np.array(hists)
    bottom = np.zeros(hists.shape[1])
    for i, cls in enumerate(classes):
        plt.bar(bins[:-1], hists[i], width=bins[1] - bins[0], bottom=bottom, color=get_color_for_class(cls),
                edgecolor='black', label=cls)
        bottom += hists[i]

    plt.title("CV Distribution")
    plt.xlabel("Coefficient of Variation (CV)")
    plt.ylabel("Count")
    # plt.axvline(x=0.2, color='r', linestyle='--', label='CV = 0.2')
    plt.legend()

def cv_pipeline(filepath, results_set=None, mode="Products", drop_IS=True, normalize_IS=True, normalize_TIC=True, plot_results=True,
                figsavepath=None, nbins=50, cleanup=True, write_output=False):
    df = pd.read_csv(filepath)
    # Drop all but PC and EtherPC classes to debug
    # df = df[df["Molecule List Name"].isin(["PC", "EtherPC"])]

    if cleanup:
        # Drop replicate name or Area nan
        df = df.dropna(subset=["Replicate Name", "Area"])

    normdf, replicates = sum_transitions(df, mode=mode, drop_IS=drop_IS, normalize_IS=normalize_IS,
                                         normalize_TIC=normalize_TIC)

    # Sort NormDF by molecule list name, then by molecule name, then by adduct
    normdf = normdf.sort_values(by=["Molecule List Name", "Molecule", "Adduct"], ascending=[True, True, True]).reset_index(drop=True)

    if results_set is not None:
        # Take only replicates that are in results_set
        replicates = [r for r in replicates if r in results_set]

    if write_output:
        # Write normdf to xlsx
        outfile = os.path.splitext(filepath)[0] + "_cv_output.xlsx"
        normdf.to_excel(outfile, index=False)
        print(f"CV output written to {outfile}")
        print(normdf.to_string())

    if plot_results:
        # Plot CV distribution
        plt.figure(figsize=(8, 6))
        ncol=2
        nrow=1
        ax1 = plt.subplot(nrow, ncol, 1)
        # normdf = simplify_class_df(normdf, "Molecule List Name")
        plot_cv_hist(normdf, ax=ax1, nbins=nbins)
        ax2 = plt.subplot(nrow, ncol, 2)
        # Drop classes with CV greater than 1 for pie chart to avoid clutter
        normdf = normdf[normdf["CV"] < 1]
        pie_chart(normdf, "Mean", ax=ax2)
        plt.tight_layout()
        if figsavepath is not None:
            plt.savefig(figsavepath)
        plt.show()


def filter_normdf(normdf, fold_range=1.):
    normdfsorted = normdf.sort_values(by=["Log2FC"], ascending=[False])

    # Filter normdfsorted to only include rows where -log10p > 1.3 (p<0.05) and abs(log2FC) > 1
    normdfsorted = normdfsorted[
        (normdfsorted["Log10p"] > -np.log10(0.05)) & (normdfsorted["Log2FC"].abs() > fold_range)]
    if "Comparison status" in normdf:
        # Keep detection changes even when the statistical cutoff is not met.
        retain = normdf["Comparison status"].ne("Tested") | normdf["Detection change"].ne("Detected in both")
        normdfsorted = pd.concat([normdfsorted, normdf.loc[retain]]).drop_duplicates(subset=["Molecule"])
    if len(normdfsorted) == 0:
        print("No significant molecules found for heatmap after applying cutoffs.")
        # normdfsorted = normdf.sort_values(by=["Molecule List Name", "Mean"], ascending=[True, False])
        normdfsorted = normdf.sort_values(by=["Log2FC"], ascending=[False])
    return normdfsorted

def select_internal_standards(df, overrides, replicates, exclusion_audit=None):
    """Apply manually reviewed class-specific IS choices identically to both groups."""
    keep = pd.Series(True, index=df.index)
    for lipid_class, standard in (overrides or {}).items():
        if not isinstance(standard, str) or not standard.startswith("IS "):
            raise ValueError(f"{lipid_class}: specify the exact internal-standard Molecule name")
        class_rows = df["Molecule List Name"].eq(lipid_class)
        selected = df.loc[class_rows & df["Molecule"].eq(standard)]
        missing = set(replicates).difference(selected["Replicate Name"].dropna())
        if missing:
            raise ValueError(f"Selected {lipid_class} standard {standard} is missing product rows in: "
                             + ", ".join(sorted(missing)))
        excluded = class_rows & df["Molecule"].str.startswith("IS ") & df["Molecule"].ne(standard)
        if exclusion_audit is not None:
            for name, group in df.loc[excluded].groupby("Molecule", sort=False):
                exclusion_audit.append({"Molecule": name, "Class": lipid_class, "Adduct": "All",
                    "Reason": f"Reviewed IS selection: use {standard} in both comparison groups.",
                    "Stage": "IS selection", "Rows removed": len(group), "Row unit": "Input result rows",
                    "Replicates": "; ".join(sorted(group["Replicate Name"].unique()))})
        keep &= ~excluded
        print(f"Using reviewed {lipid_class} standard {standard} in both comparison groups.")
    return df.loc[keep].copy()


def compare_pipeline(file1, set1, set2, mode="Products", drop_IS=True, normalize_IS=True,
                  normalize_TIC=True, norm_tmm=False,
                     paired=True, bh_correction=True, plot_results=True, heatplots=True,
                  piecharts=True, figsavepath=None, heatmap_cutoffs=True, write_output=True, fold_range=1.5,
                     otherthresh=0.035, set1_name="Set 1", set2_name="Set 2",
                     show_plots=True,
                     is_standard_overrides=None, file2=None, missing_set2_as_zero=True,
                     cholesterol_file=None):
    """Compare replicate groups from one or two Skyline results CSVs.

    set1/set2 can each be {"file1": [replicate names], "file2": [replicate names]}.
    Either group can select samples from either or both files. A plain list
    searches both CSVs and requires each name to occur in only one file.
    File prefixes distinguish samples with identical names in different CSVs.
    set1_name/set2_name are the freely editable group labels for plots/exports.
    Set show_plots=False to build the dashboard and defer display until a later
    plt.show(), allowing a separate summary figure to be saved first.
    is_standard_overrides maps a class to one reviewed IS Molecule name to use
    in both groups. Other standards in that class are excluded and audited.
    missing_set2_as_zero treats absent input rows in set2 as zero signal when
    the lipid is present in set1. Quantification failures remain visible as missing estimates.
    Tests then include zeros, with this assumption recorded in Test basis.
    cholesterol_file optionally adds precursor transitions from a separate CSV
    before normalization. Its Replicate Name values must match selected samples;
    names shared by selected samples in both input files are ambiguous and rejected.
    """
    file1 = os.path.abspath(os.fspath(file1))
    paths = {"file1": file1}
    if file2 is not None:
        paths["file2"] = os.path.abspath(os.fspath(file2))
    sources = {label: pd.read_csv(path) for label, path in paths.items()}
    available = {label: set(frame["Replicate Name"].dropna()) for label, frame in sources.items()}
    selections = {label: {} for label in paths}
    used_samples = set()
    resolved_groups = []
    sample_selection = []
    for group_number, requested, group_name in [(1, set1, set1_name), (2, set2, set2_name)]:
        samples = []
        if isinstance(requested, dict):
            for label, names in requested.items():
                if label not in paths:
                    if not names:
                        continue
                    raise ValueError(f"No input configured for {label}")
                if not isinstance(names, (list, tuple)):
                    raise ValueError(f"Set{group_number} {label}: provide a list of replicate names")
                samples.extend((label, name) for name in names)
        else:
            if not isinstance(requested, (list, tuple)):
                raise ValueError(f"Set{group_number}: provide a replicate list or file-selection dictionary")
            for name in requested:
                matches = [label for label in paths if name in available[label]]
                if not matches:
                    raise ValueError(f"Comparison replicates not found in CSV: {name}")
                if len(matches) > 1:
                    raise ValueError(f"Replicate {name} occurs in both CSVs; select file1/file2 explicitly")
                samples.append((matches[0], name))
        if not samples:
            raise ValueError(f"Set{group_number} must select at least one sample")
        resolved = []
        for label, name in samples:
            if name not in available[label]:
                raise ValueError(f"Comparison replicates not found in CSV: {name} ({paths[label]})")
            if (label, name) in used_samples:
                raise ValueError(f"Sample {label}: {name} is selected more than once")
            used_samples.add((label, name))
            alias = f"{label}::{name}" if len(paths) > 1 else name
            selections[label][name] = alias
            resolved.append(alias)
            sample_selection.append({"Group": f"Set{group_number}", "Group label": group_name,
                                     "File": label, "Input": paths[label],
                                     "Original replicate": name, "Comparison replicate": alias})
        resolved_groups.append(resolved)
    set1, set2 = resolved_groups
    frames = []
    exclusion_audit = []
    for label, selected in selections.items():
        if not selected:
            continue
        input_file = paths[label]
        frame = sources[label].loc[sources[label]["Replicate Name"].isin(selected)].copy()
        frame["Original Replicate Name"] = frame["Replicate Name"]
        frame["Replicate Name"] = frame["Replicate Name"].map(selected)
        frame["Input Source"] = input_file
        frames.append(frame)
    df = pd.concat(frames, ignore_index=True)
    input_label = "; ".join(paths[label] for label, selected in selections.items() if selected)
    source_df = df.copy()
    df = filter_transition_mode(df, mode)
    if cholesterol_file is not None:
        cholesterol_file = os.path.abspath(os.fspath(cholesterol_file))
        cholesterol = pd.read_csv(cholesterol_file)
        aliases = {}
        for selected in selections.values():
            for name, alias in selected.items():
                if name in aliases and cholesterol["Replicate Name"].eq(name).any():
                    raise ValueError(f"Cholesterol replicate {name} matches selected samples in both CSVs")
                aliases[name] = alias
        cholesterol = cholesterol.loc[cholesterol["Replicate Name"].isin(aliases)].copy()
        if cholesterol.empty:
            raise ValueError("No cholesterol replicates match the selected samples")
        cholesterol["Original Replicate Name"] = cholesterol["Replicate Name"]
        cholesterol["Replicate Name"] = cholesterol["Replicate Name"].map(aliases)
        cholesterol["Input Source"] = cholesterol_file
        # Some APCI exports put the d7 precursor under the analyte's name.
        # Separate it before sum_transitions removes isotope labels from adducts.
        d7_standard = (cholesterol["Molecule"].eq("ST 27:1;1")
                       & cholesterol["Precursor Adduct"].eq("[M7H2+H-H2O]"))
        cholesterol.loc[d7_standard, "Molecule"] = "IS ST 27:1;1(+[2]H7)"
        source_df = pd.concat([source_df, cholesterol], ignore_index=True)
        cholesterol = filter_transition_mode(cholesterol, "Precursors")
        if cholesterol.empty:
            raise ValueError("No cholesterol precursor transitions remain for the selected samples")
        df = pd.concat([df, cholesterol], ignore_index=True)
        input_label += "; " + cholesterol_file

    notes_path = os.path.splitext(file1)[0] + "_quantification_notes.csv"
    if write_output:
        # Create the notes file before normalization, even if quantification later fails.
        pd.DataFrame(exclusion_audit, columns=QUANTIFICATION_NOTE_COLUMNS).to_csv(notes_path, index=False)
    if df.empty:
        raise ValueError("No result rows remain for the selected transition mode")
    missing = set(set1 + set2).difference(df["Replicate Name"].dropna())
    if missing:
        raise ValueError("Comparison replicates not found after selection: " + ", ".join(sorted(missing)))

    # Limit to molecule list name of pc for testing
    # df = df[df["Molecule List Name"] == "PC"]

    if is_standard_overrides:
        if not normalize_IS:
            raise ValueError("is_standard_overrides requires normalize_IS=True")
        df = select_internal_standards(df, is_standard_overrides, set1 + set2, exclusion_audit)
    normdf, replicates = sum_transitions(df, mode="All", drop_IS=drop_IS, normalize_IS=normalize_IS,
                                         normalize_TIC=normalize_TIC, is_replicate_groups=[set1, set2],
                                         exclusion_audit=exclusion_audit)
    if norm_tmm:
        print("TMM Normalization")
        normdf = normalize_trimmed_mean_mvalues(normdf, set1, set2)

    # Keep missing quantification estimates and per-replicate PC adduct choices
    # reviewable even when downstream statistics exclude incomplete molecules.
    quantified = detection_table(normdf, source_df,
                                 set1 + set2, drop_IS=drop_IS)
    quantified["Set2_AssignedZeros"] = 0
    if missing_set2_as_zero:
        control_present = quantified[[f"Detection {r}" for r in set1]].ne("Not measured").any(axis=1)
        for r in set2:
            absent = control_present & quantified[f"Detection {r}"].eq("Not measured")
            quantified.loc[absent, r] = 0.
            quantified.loc[absent, f"Detection {r}"] = "Absent in set2; assigned zero"
            quantified.loc[absent, "Set2_AssignedZeros"] += 1
    classdf, classsummarycols = class_calcs(quantified, set1, set2, paired=paired,
                                           bh_correction=bh_correction, retain_unavailable=True,
                                           include_nondetections=missing_set2_as_zero)
    normdf, summarycols = comparison_stats(quantified, set1, set2, paired=paired,
                                           bh_correction=bh_correction, include_nondetections=missing_set2_as_zero)
    normdf.attrs["bh_correction"] = bh_correction
    normdf["Major difference"] = normdf["Log2FC"].abs().ge(fold_range)
    normdf["Significant difference"] = normdf["Major difference"] & normdf["p-value"].lt(0.05)
    normdf["Set1 name"] = set1_name
    normdf["Set2 name"] = set2_name
    if is_standard_overrides:
        normdf["Configured IS"] = normdf["Molecule List Name"].map(is_standard_overrides)
    print(len(normdf), "molecules after stats calculations.")
    unavailable = quantified.loc[~quantified["Molecule"].isin(normdf.loc[normdf["Comparison status"].eq("Tested"), "Molecule"])]
    for _, row in unavailable.iterrows():
        missing_replicates = [r for r in set1 + set2 if not np.isfinite(row[r])]
        reason = ("Retained without t-test: missing quantification estimates."
                  if missing_replicates else "Retained without t-test: non-detections or undefined t-test; see Detection Changes.")
        exclusion_audit.append({"Input": input_label, "Molecule": row["Molecule"],
                               "Class": row["Molecule List Name"], "Adduct": row["Adduct"],
                               "Reason": reason, "Stage": "Comparison statistics only",
                               "Rows removed": 0, "Row unit": "Retained molecule; not t-tested",
                               "Replicates": "; ".join(missing_replicates)})

    normdf = split_names(normdf)
    print(normdf.to_string())

    if write_output:
        detection_changes = normdf.loc[normdf["Comparison status"].ne("Tested")
                                       | normdf["Detection change"].ne("Detected in both")
                                       | normdf["Test basis"].eq("Includes zero non-detections")]
        detection_changes.to_csv(os.path.splitext(file1)[0] + "_detection_changes.csv", index=False)
        for record in exclusion_audit:
            record.setdefault("Input", input_label)
        notes = pd.DataFrame(exclusion_audit, columns=QUANTIFICATION_NOTE_COLUMNS)
        notes.to_csv(notes_path, index=False)
        print(f"Saved quantification notes to {notes_path}")
        try:
            # Write to two sheets in an Excel file
            with pd.ExcelWriter(os.path.splitext(file1)[0] + "_stats.xlsx") as writer:
                pd.DataFrame(sample_selection).to_excel(writer, sheet_name="Sample Selection", index=False)
                classdf.to_excel(writer, sheet_name="Class Stats", index=False)
                normdf.to_excel(writer, sheet_name="Molecule Stats", index=False)
                detection_changes.to_excel(writer, sheet_name="Detection Changes", index=False)
                quantified.to_excel(writer, sheet_name="Quantified Molecules", index=False)
                notes.to_excel(writer, sheet_name="Quantification Notes", index=False)
                df.to_excel(writer, sheet_name="Raw Data", index=False)
        except:
            print("ERROR writing to Excel file. Close file and try again!")

    # if write_output:
    #     # # Write classdf to CSV
    #     classdf.to_csv(os.path.splitext(filepath)[0] + "_class_stats.csv", index=False)
    #     # # Write normdf to CSV
    #     normdf.to_csv(os.path.splitext(filepath)[0] + "_molecule_stats.csv", index=False)

    if plot_results:

        nrow = 1
        if heatplots:
            nrow +=1
        if piecharts:
            nrow +=1
        ncol = 3

        fig = plt.figure(figsize=(18, 10))
        gs = GridSpec(nrow, ncol, figure=fig)

        # ax1 = plt.subplot(nrow, ncol, 1)
        ax1 = fig.add_subplot(gs[0, 0])
        make_class_plots(classdf, ax=ax1)

        # ax2 = plt.subplot(nrow, ncol, 2)
        ax2 = fig.add_subplot(gs[0, 1])
        make_volcano_plot(normdf, title="Volcano Plot", ax=ax2, fold_range=fold_range, show_detection_summary=False)

        # ax3 = plt.subplot(nrow, ncol, 3)
        ax3 = fig.add_subplot(gs[0, 2])
        pca_plot(normdf, set1, set2, title="PCA Plot", ax=ax3)

        i=3
        if heatplots:
            # ax4 = plt.subplot(nrow, ncol, i+1)
            ax4 = fig.add_subplot(gs[1, 0])
            # Sort normdf by class and then by Mean
            if heatmap_cutoffs:
                normdfsorted = filter_normdf(normdf, fold_range=fold_range)
            else:
                normdfsorted = normdf.sort_values(by=["Molecule List Name", "Mean"], ascending=[True, False])

            heatmap_plot(normdfsorted, replicates, title="Heatmap", ax=ax4, normtype="minmax", scale="linear")

            # ax5 = plt.subplot(nrow, ncol, i+2)
            ax5 = fig.add_subplot(gs[1, 1])
            if paired:
                titleval = "Paired Ratio Heatmap"
                normtype = "onecentered"
                cmap = "RdBu_r"
            else:
                titleval = "Set Means Heatmap"
                normtype = "minmax"
                cmap = "viridis"
            heatmap_plot(normdfsorted, summarycols, title=titleval, ax=ax5, lindex=0, colormap=cmap, normtype=normtype)

            ax6 = fig.add_subplot(gs[1:, 2])
            cols = ["Log2FC", "StdDevFC"]
            heatmap_plot(normdfsorted, cols, title="Fold Change Heatmap", ax=ax6, lindex=0, colormap="RdBu_r",
                         normtype="zerocentered", fontsize=8)
            i += 3

        if piecharts:
            # ax7 = plt.subplot(nrow, ncol, i+1)
            ax7 = fig.add_subplot(gs[-1, 0])
            pie_chart(classdf, "Set1_Mean", ax=ax7, title="Class Distribution " + set1_name, otherthresh=otherthresh)
            # ax8 = plt.subplot(nrow, ncol, i+2)
            ax8 = fig.add_subplot(gs[-1, 1])
            pie_chart(classdf, "Set2_Mean", ax=ax8, title="Class Distribution " + set2_name, otherthresh=otherthresh)

        plt.tight_layout()

        if figsavepath is not None:
            plt.savefig(figsavepath)

        if show_plots:
            plt.show()

    return normdf, classdf

def is_analysis(file, mode="Products", isdf=None):
    df = pd.read_csv(file)
    isrows = df[df["Molecule"].str.startswith("IS ")]
    isrows, replicates = sum_transitions(isrows, mode=mode, drop_IS=False, normalize_IS=False, normalize_TIC=False)

    if isdf is None:
        isdf = default_isdf

    # Calculate response factors for each IS compound in each replicate
    for i, row in isrows.iterrows():
        isname = row["Molecule"].strip().replace("IS ", "")
        isclass = row["Molecule List Name"]
        conc = get_is_conc_from_isdf(row, isdf, "Concentration (uM)")
        # Create new columns of rep + _"RF" for response factor
        if conc is not None:
            for r in replicates:
                response_factor = row[r] / conc if conc != 0 else 0
                isrows.at[i, r + "_RF"] = response_factor
        else:
            print(f"IS Compound: {isname}, Class: {isclass}, Concentration: Unknown. Response Factor: Unknown (no matching IS found in IS dataframe)")

    # Mean response factor across replicates for each IS compound
    rf_cols = [r + "_RF" for r in replicates]
    isrows["Mean_RF"] = isrows[rf_cols].mean(axis=1)
    # CV rt
    isrows["CV_RF"] = isrows[rf_cols].std(axis=1, ddof=1) / isrows[rf_cols].mean(axis=1)
    print(isrows.to_string())

def lipid_bar_chart_compare(df, classes=["PC"], ax=None, fontsize=8, s1name="Set 1", s2name="Set 2", log2fclimit=0):
    if ax is None:
        plt.figure(figsize=(8, 6))
    else:
        plt.sca(ax)

    subdf = df[df["Molecule List Name"].isin(classes)]
    # Sort by molecule name
    subdf = subdf.sort_values(by=["Molecule"], ascending=[True])
    # Means are Set1_Mean and Set2_Mean columns
    # Errors are Set1_CV and Set2_CV columns multiplied by the means to get standard deviation
    x = np.arange(len(subdf))
    labels = subdf["Molecule"]
    # Drop before "|" in labels for better readability
    labels = clean_names(labels, split_side=1)
    labels = fix_name_order(labels)
    width = 0.35
    plt.bar(x - width/2, subdf["Set1_Mean"], width, yerr=subdf["Set1_CV"] * subdf["Set1_Mean"], label=s1name, capsize=5)
    plt.bar(x + width/2, subdf["Set2_Mean"], width, yerr=subdf["Set2_CV"] * subdf["Set2_Mean"], label=s2name, capsize=5)
    bmax = max((subdf["Set1_Mean"] + subdf["Set1_CV"] * subdf["Set1_Mean"]).max(), (subdf["Set2_Mean"] + subdf["Set2_CV"] * subdf["Set2_Mean"]).max())
    # Add significance stars based on p-value with thresholds of 0.05, 0.01, and 0.001
    offset = bmax * 0.02
    for i, pval in enumerate(subdf["p-value"]):
        log2fc = subdf["Log2FC"].iloc[i]
        if abs(log2fc) < log2fclimit:
            continue

        if pval < 0.001:
            plt.text(x[i], max(subdf["Set1_Mean"].iloc[i] + subdf["Set1_CV"].iloc[i] * subdf["Set1_Mean"].iloc[i],
                               subdf["Set2_Mean"].iloc[i] + subdf["Set2_CV"].iloc[i] * subdf["Set2_Mean"].iloc[i]) + offset,
                     "***", ha='center', va='bottom', fontsize=fontsize)
        elif pval < 0.01:
            plt.text(x[i], max(subdf["Set1_Mean"].iloc[i] + subdf["Set1_CV"].iloc[i] * subdf["Set1_Mean"].iloc[i],
                               subdf["Set2_Mean"].iloc[i] + subdf["Set2_CV"].iloc[i] * subdf["Set2_Mean"].iloc[i]) + offset,
                     "**", ha='center', va='bottom', fontsize=fontsize)
        elif pval < 0.05:
            plt.text(x[i], max(subdf["Set1_Mean"].iloc[i] + subdf["Set1_CV"].iloc[i] * subdf["Set1_Mean"].iloc[i],
                               subdf["Set2_Mean"].iloc[i] + subdf["Set2_CV"].iloc[i] * subdf["Set2_Mean"].iloc[i]) + offset,
                     "*", ha='center', va='bottom', fontsize=fontsize)

    plt.xticks(x, labels, rotation=90, fontsize=fontsize)
    plt.legend()
    plt.ylabel("Mol %")
    plt.title("Lipid Comparison for " + ", ".join(classes))



if __name__ == "__main__":
    file1 = r"Z:\Group Share\Annika\Stellar\EC\B2 Discs\Detergent trials\B2_vs_HEK_results.csv"
    file2 = None
    output_dir = Path(r"Z:\Group Share\Annika\Stellar\EC\B2 Discs\Comparison Results")
    output_name = "HEK_B2_withcholesterol_volcano_pies"
    set1_name = "HEK Extract"
    set2_name = "B2 Overexpressed Extract"
    cholesterol_file = r"Z:\Group Share\Annika\Stellar\EC\B2 Discs\APCI\Chol_TransResults_B2Extr.csv"  # Set to a cholesterol Skyline CSV with matching replicate names.
    # cholesterol_file = r"Z:\Group Share\Annika\Stellar\FAM\B2 Discs\Cholesterol\Molecule Transition Results_cholesterol.csv"
    # file1 = r"C:\Users\marty\Downloads\Extr_US_MoleculeTransResults.csv"

    # set1 = ["noEDTA_E_01", "noEDTA_E_02", "noEDTA_E_03"]
    # set2 = ["EDTA_E_01", "EDTA_E_02", "EDTA_E_03"]

    # set1 = ["noEDTA_FT_01", "noEDTA_FT_02", "noEDTA_FT_03"]
    # set2 = ["noEDTA_E_01", "noEDTA_E_02", "noEDTA_E_03"]
    #
    # set1 = ["EDTA_FT_01", "EDTA_FT_02", "EDTA_FT_03"]
    # set2 = ["EDTA_E_01", "EDTA_E_02", "EDTA_E_03"]

    # Select exact Replicate Name values from either or both files for each group.
    # Leave a file's list empty to omit its samples from that group.
    set1 = {"file1": ["HEK_01", "HEK_02", "HEK_03"]}
    set2 = {"file1": ["B2_01", "B2_02", "B2_03"]}
    # Count missing analyte rows as zero detected signal; preserve measured low values.
    missing_set2_as_zero = True
    fold_range = np.log2(2)  # Twofold effect-size cutoff, shared by plots and export flags.
    # Only this LPI standard passed manual review in EDTA; use it in B2 too.
    # is_standard_overrides = {"LPI": "IS LPI(+[2]H5) 17:0"}


    # is_analysis(file1)
    # exit()

    # cv_pipeline(file1, results_set=set2["file1"], mode="Products", drop_IS=False, normalize_IS=False, normalize_TIC=False,
    #             write_output=True, plot_results=True)
    #
    # exit()

    normdf, classdf = compare_pipeline(file1, set1, set2, mode="Products", file2=file2,
                                       cholesterol_file=cholesterol_file,
                                       missing_set2_as_zero=missing_set2_as_zero,
                                       drop_IS=True, normalize_IS=True, norm_tmm=False,
                  normalize_TIC=True, paired=False, bh_correction=True, plot_results=True, show_plots=False, write_output=True,
                     fold_range=fold_range, otherthresh=0.03, set1_name=set1_name, set2_name=set2_name)



    # Figure with pie 1 volcano pie 2 in a row
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))
    pie_chart(classdf, "Set1_Mean", ax=ax1, title=set1_name, otherthresh=0.015)
    make_volcano_plot(normdf, title="Volcano Plot", ax=ax2, fold_range=fold_range, xlims=[-6, 6],
                      show_detection_summary=False, show_group_label=False)
    pie_chart(classdf, "Set2_Mean", ax=ax3, title=set2_name, otherthresh=0.015)
    fig.tight_layout()
    output_dir.mkdir(parents=True, exist_ok=True)
    for extension in ["png", "pdf"]:
        output_path = output_dir / f"{output_name}.{extension}"
        fig.savefig(output_path, dpi=600)
        print(f"Saved figure to {output_path}")
    # Display both the full dashboard and the saved three-panel summary.
    plt.show()
