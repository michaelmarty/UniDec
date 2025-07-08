//
// Created by Michael Marty on 7/8/2025.
//

#ifndef UDMAIN_H
#define UDMAIN_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "hdf5.h"

#include "udmain.h"
#include "udcore.h"
#include "udtools.h"
#include "udio.h"
#include "h5io.h"
// #include "UD_IO.h"
#include "UD_score.h"
#include "udstruct.h"

Decon ExitToBlank(const Config config, Decon decon);
Decon MainDeconvolution(const Config config, const Input inp, const int silent, const int verbose);
void RunAutotune(Config *config, Input *inp, Decon *decon);
int run_unidec(int argc, char *argv[], Config config);

#endif //UDMAIN_H
