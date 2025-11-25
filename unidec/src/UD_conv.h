/*
* UD_conv.h
*
*  Created on : 11/06/22
* Author : Michael.Marty
*/

//
// 
// Copyright 2022 University of Arizona
//
//

#ifndef UD_CONV_H
#define UD_CONV_H
#include "udcore.h"
#include "udtools.h"
#include "udio.h"
#include "udstruct.h"

void create_sim_spec(Config config, Input inp, const float* blur, const float mass, float* simint);
int conv_main(int argc, char* argv[], Config config);

#endif
