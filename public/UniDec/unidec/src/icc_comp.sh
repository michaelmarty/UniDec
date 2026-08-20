#icc unidec.c -o unideclinux -mAVX -I/cm/shared/apps/hdf5_18/1.8.16/include -L/cm/shared/apps/hdf5_18/1.8.16/lib -lm -lhdf5_hl -lhdf5 -lz -mkl=parallel -qopenmp -openmp-link static -fast -O3 -std=c99 -w

icc UniDec.c -o unideclinux -mAVX -I/cm/shared/apps/hdf5_18/1.8.16/include -L/cm/shared/apps/hdf5_18/1.8.16/lib -lm -lhdf5_hl -lhdf5 -lz -mkl=parallel -qopenmp -O3 -std=c99 -w
