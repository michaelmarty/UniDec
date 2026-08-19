import unidec.tools as ud
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib_venn import venn2
from matplotlib_venn import venn3
import matplotlib as mpl
from matplotlib import rcParams

from UniDecImporter.ImporterFactory import ImporterFactory

# Set backend to Agg
mpl.use('WxAgg')
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'
if False:
    shared = 2828
    isodec_unique = 1303
    topfd_unique = 866

    venn = venn2(subsets=(isodec_unique, topfd_unique, shared), set_labels=("IsoDec", "TopFD"))
    venn.get_patch_by_id('10').set_color('red')
    venn.get_patch_by_id('01').set_color('blue')
    venn.get_patch_by_id('11').set_color('purple')

    label10 = venn.get_label_by_id('10')
    label10.set_position((-0.5, 0))
    label10.set_size(20)
    # plt.plot([-0.63,-0.55], [0,0], color="black", lw=2)


    label01 = venn.get_label_by_id('01')
    label01.set_position((0.5, 0))
    label01.set_size(20)
    # plt.plot([0.55,0.63], [0,0], color="black", lw=2)

    label11 = venn.get_label_by_id('11')
    label11.set_position((0, 0))
    label11.set_size(20)

    label_id = venn.get_label_by_id("A")
    label_id.set_size(20)
    label_id.set_position((-0.25,-0.55))

    label_tfd = venn.get_label_by_id("B")
    label_tfd.set_size(20)
    label_tfd.set_position((0.25,-0.55))
    plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\isodec_TopFD_prsms_noscanz.png", dpi=300, transparent=True)
    plt.savefig("Z:\\Group Share\\JGP\\JPST001885\\CA_UVPD_forID_manuscript\\Figures\\isodec_TopFD_prsms_noscanz.pdf", dpi=300, transparent=True)
    plt.show()

if True:
    unique_tfd = 1946
    unique_fd = 902
    tfd_fd = 810
    unique_iso = 868
    tfd_iso = 699
    fd_iso = 640
    all_three = 4565

    venn = venn3(subsets=(unique_tfd,unique_fd,tfd_fd,unique_iso,tfd_iso,fd_iso,all_three), set_labels=("TFD", "FLD", "ID-8"))
    venn.get_patch_by_id('100').set_color('red')
    venn.get_patch_by_id('010').set_color('yellow')
    venn.get_patch_by_id('001').set_color('blue')
    venn.get_patch_by_id('110').set_color('orange')
    venn.get_patch_by_id('101').set_color('purple')
    venn.get_patch_by_id('011').set_color('green')
    venn.get_patch_by_id('111').set_color('gray')


    venn.get_patch_by_id('100').set_edgecolor('black')
    venn.get_patch_by_id('010').set_edgecolor('black')
    venn.get_patch_by_id('001').set_edgecolor('black')
    venn.get_patch_by_id('110').set_edgecolor('black')
    venn.get_patch_by_id('101').set_edgecolor('black')
    venn.get_patch_by_id('011').set_edgecolor('black')
    venn.get_patch_by_id('111').set_edgecolor('black')

    unique_alpha = 0.8
    shared_alpha = 0.8

    venn.get_patch_by_id('100').set_alpha(unique_alpha)
    venn.get_patch_by_id('010').set_alpha(unique_alpha)
    venn.get_patch_by_id('001').set_alpha(unique_alpha)

    venn.get_patch_by_id('110').set_alpha(shared_alpha)
    venn.get_patch_by_id('101').set_alpha(shared_alpha)
    venn.get_patch_by_id('011').set_alpha(shared_alpha)
    venn.get_patch_by_id('111').set_alpha(shared_alpha)


    fontsize = 18

    venn.get_label_by_id('100').set_size(fontsize)
    venn.get_label_by_id('010').set_size(fontsize)
    venn.get_label_by_id('001').set_size(fontsize)
    venn.get_label_by_id('110').set_size(fontsize)
    venn.get_label_by_id('101').set_size(fontsize)
    venn.get_label_by_id('011').set_size(fontsize)
    venn.get_label_by_id('111').set_size(fontsize)

    # venn.get_label_by_id('100').set_fontweight('bold')
    # venn.get_label_by_id('010').set_fontweight('bold')
    # venn.get_label_by_id('001').set_fontweight('bold')
    # venn.get_label_by_id('110').set_fontweight('bold')
    # venn.get_label_by_id('101').set_fontweight('bold')
    # venn.get_label_by_id('011').set_fontweight('bold')
    # venn.get_label_by_id('111').set_fontweight('bold')

    # venn.get_label_by_id('110').set_position((-0.18,0.38))
    # venn.get_label_by_id('010').set_position((0.35,0.35))
    # venn.get_label_by_id('011').set_position((0.29,-0.30))
    # venn.get_label_by_id('001').set_position((0,-0.52))
    # venn.get_label_by_id('100').set_position((-0.71,0.0))
    # venn.get_label_by_id('101').set_position((-0.50,-0.45))
    #
    # plt.plot([-0.47, -0.35], [-0.40,-0.30], color='k', lw=2)
    # plt.plot([-0.55, -0.62], [0,0], color='k', lw=2)

    labelsize = 24
    venn.get_label_by_id('A').set_size(labelsize)
    venn.get_label_by_id('B').set_size(labelsize)
    venn.get_label_by_id('C').set_size(labelsize)

    plt.tight_layout()
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\toppic_venn3_scanoverlap.png",
                dpi=300, transparent=True)
    plt.savefig("C:\\Users\\MartyLabPC\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures\\toppic_venn3_scanoverlap.pdf",
                dpi=300, transparent=True)
    plt.show()


