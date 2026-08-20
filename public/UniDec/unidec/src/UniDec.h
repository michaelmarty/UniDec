/*
 * UniDec.h
 *
 *  Created on: 1 Jul 2014
 *      Author: Michael.Marty
 */

/*
#ifndef UNIDEC_H_
#define UNIDEC_H_
#endif*/ /* UNIDEC_H_ */

#ifndef UNIDEC_HEADER
#define UNIDEC_HEADER

#define __STDC_WANT_LIB_EXT1__ 1
#include <stdio.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>

#include "udmain.h"
#include "UniDecIM_Main.h"
#include "MetaUniDec_Main.h"
#include "UniDecCD_Main.h"
#include "UD_conv.h"

Config ImportConfig(int argc, char * argv[], Config config);

#endif
