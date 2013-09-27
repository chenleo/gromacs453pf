/*
 *                This source code is part of
 *
 *                 G   R   O   M   A   C   S
 *
 * Copyright (c) 1991-2000, University of Groningen, The Netherlands.
 * Copyright (c) 2001-2009, The GROMACS Development Team
 *
 * Gromacs is a library for molecular simulation and trajectory analysis,
 * written by Erik Lindahl, David van der Spoel, Berk Hess, and others - for
 * a full list of developers and information, check out http://www.gromacs.org
 *
 * This program is free software; you can redistribute it and/or modify it under 
 * the terms of the GNU Lesser General Public License as published by the Free 
 * Software Foundation; either version 2 of the License, or (at your option) any 
 * later version.
 * As a special exception, you may use this file as part of a free software
 * library without restriction.  Specifically, if other files instantiate
 * templates or use macros or inline functions from this file, or you compile
 * this file and link it with other files to produce an executable, this
 * file does not by itself cause the resulting executable to be covered by
 * the GNU Lesser General Public License.  
 *
 * In plain-speak: do not worry about classes/macros/templates either - only
 * changes to the library have to be LGPL, not an application linking with it.
 *
 * To help fund GROMACS development, we humbly ask that you cite
 * the papers people have written on it - you can find them on the website!
 */
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "types/nrnb.h"
#include "pf_nb_kernel_c.h"
#include "../nb_kerneltype.h"


/* Include standard kernel headers in local directory */
#include "pf_nb_kernel010.h"
#include "pf_nb_kernel020.h"
#include "pf_nb_kernel030.h"
#include "pf_nb_kernel100.h"
#include "pf_nb_kernel101.h"
#include "pf_nb_kernel102.h"
#include "pf_nb_kernel103.h"
#include "pf_nb_kernel104.h"
#include "pf_nb_kernel110.h"
#include "pf_nb_kernel111.h"
#include "pf_nb_kernel112.h"
#include "pf_nb_kernel113.h"
#include "pf_nb_kernel114.h"
#include "pf_nb_kernel120.h"
#include "pf_nb_kernel121.h"
#include "pf_nb_kernel122.h"
#include "pf_nb_kernel123.h"
#include "pf_nb_kernel124.h"
#include "pf_nb_kernel130.h"
#include "pf_nb_kernel131.h"
#include "pf_nb_kernel132.h"
#include "pf_nb_kernel133.h"
#include "pf_nb_kernel134.h"
#include "pf_nb_kernel200.h"
#include "pf_nb_kernel201.h"
#include "pf_nb_kernel202.h"
#include "pf_nb_kernel203.h"
#include "pf_nb_kernel204.h"
#include "pf_nb_kernel210.h"
#include "pf_nb_kernel211.h"
#include "pf_nb_kernel212.h"
#include "pf_nb_kernel213.h"
#include "pf_nb_kernel214.h"
#include "pf_nb_kernel220.h"
#include "pf_nb_kernel221.h"
#include "pf_nb_kernel222.h"
#include "pf_nb_kernel223.h"
#include "pf_nb_kernel224.h"
#include "pf_nb_kernel230.h"
#include "pf_nb_kernel231.h"
#include "pf_nb_kernel232.h"
#include "pf_nb_kernel233.h"
#include "pf_nb_kernel234.h"
#include "pf_nb_kernel300.h"
#include "pf_nb_kernel301.h"
#include "pf_nb_kernel302.h"
#include "pf_nb_kernel303.h"
#include "pf_nb_kernel304.h"
#include "pf_nb_kernel310.h"
#include "pf_nb_kernel311.h"
#include "pf_nb_kernel312.h"
#include "pf_nb_kernel313.h"
#include "pf_nb_kernel314.h"
#include "pf_nb_kernel320.h"
#include "pf_nb_kernel321.h"
#include "pf_nb_kernel322.h"
#include "pf_nb_kernel323.h"
#include "pf_nb_kernel324.h"
#include "pf_nb_kernel330.h"
#include "pf_nb_kernel331.h"
#include "pf_nb_kernel332.h"
#include "pf_nb_kernel333.h"
#include "pf_nb_kernel334.h"
#include "pf_nb_kernel400.h"
#include "pf_nb_kernel410.h"
#include "pf_nb_kernel420.h"
#include "pf_nb_kernel430.h"


static pf_nb_kernel_t *
pf_kernellist[eNR_NBKERNEL_NR] = 
{
    pf_nb_kernel010,
    pf_nb_kernel020,
    pf_nb_kernel030,
    pf_nb_kernel100,
    pf_nb_kernel101,
    pf_nb_kernel102,
    pf_nb_kernel103,
    pf_nb_kernel104,
    pf_nb_kernel110,
    pf_nb_kernel111,
    pf_nb_kernel112,
    pf_nb_kernel113,
    pf_nb_kernel114,
    pf_nb_kernel120,
    pf_nb_kernel121,
    pf_nb_kernel122,
    pf_nb_kernel123,
    pf_nb_kernel124,
    pf_nb_kernel130,
    pf_nb_kernel131,
    pf_nb_kernel132,
    pf_nb_kernel133,
    pf_nb_kernel134,
    pf_nb_kernel200,
    pf_nb_kernel201,
    pf_nb_kernel202,
    pf_nb_kernel203,
    pf_nb_kernel204,
    pf_nb_kernel210,
    pf_nb_kernel211,
    pf_nb_kernel212,
    pf_nb_kernel213,
    pf_nb_kernel214,
    pf_nb_kernel220,
    pf_nb_kernel221,
    pf_nb_kernel222,
    pf_nb_kernel223,
    pf_nb_kernel224,
    pf_nb_kernel230,
    pf_nb_kernel231,
    pf_nb_kernel232,
    pf_nb_kernel233,
    pf_nb_kernel234,
    pf_nb_kernel300,
    pf_nb_kernel301,
    pf_nb_kernel302,
    pf_nb_kernel303,
    pf_nb_kernel304,
    pf_nb_kernel310,
    pf_nb_kernel311,
    pf_nb_kernel312,
    pf_nb_kernel313,
    pf_nb_kernel314,
    pf_nb_kernel320,
    pf_nb_kernel321,
    pf_nb_kernel322,
    pf_nb_kernel323,
    pf_nb_kernel324,
    pf_nb_kernel330,
    pf_nb_kernel331,
    pf_nb_kernel332,
    pf_nb_kernel333,
    pf_nb_kernel334,
    pf_nb_kernel400,
    pf_nb_kernel410,
    pf_nb_kernel430,
	pf_nb_kernel010nf,
    pf_nb_kernel020nf,
    pf_nb_kernel030nf,
    pf_nb_kernel100nf,
    pf_nb_kernel101nf,
    pf_nb_kernel102nf,
    pf_nb_kernel103nf,
    pf_nb_kernel104nf,
    pf_nb_kernel110nf,
    pf_nb_kernel111nf,
    pf_nb_kernel112nf,
    pf_nb_kernel113nf,
    pf_nb_kernel114nf,
    pf_nb_kernel120nf,
    pf_nb_kernel121nf,
    pf_nb_kernel122nf,
    pf_nb_kernel123nf,
    pf_nb_kernel124nf,
    pf_nb_kernel130nf,
    pf_nb_kernel131nf,
    pf_nb_kernel132nf,
    pf_nb_kernel133nf,
    pf_nb_kernel134nf,
    pf_nb_kernel200nf,
    pf_nb_kernel201nf,
    pf_nb_kernel202nf,
    pf_nb_kernel203nf,
    pf_nb_kernel204nf,
    pf_nb_kernel210nf,
    pf_nb_kernel211nf,
    pf_nb_kernel212nf,
    pf_nb_kernel213nf,
    pf_nb_kernel214nf,
    pf_nb_kernel220nf,
    pf_nb_kernel221nf,
    pf_nb_kernel222nf,
    pf_nb_kernel223nf,
    pf_nb_kernel224nf,
    pf_nb_kernel230nf,
    pf_nb_kernel231nf,
    pf_nb_kernel232nf,
    pf_nb_kernel233nf,
    pf_nb_kernel234nf,
    pf_nb_kernel300nf,
    pf_nb_kernel301nf,
    pf_nb_kernel302nf,
    pf_nb_kernel303nf,
    pf_nb_kernel304nf,
    pf_nb_kernel310nf,
    pf_nb_kernel311nf,
    pf_nb_kernel312nf,
    pf_nb_kernel313nf,
    pf_nb_kernel314nf,
    pf_nb_kernel320nf,
    pf_nb_kernel321nf,
    pf_nb_kernel322nf,
    pf_nb_kernel323nf,
    pf_nb_kernel324nf,
    pf_nb_kernel330nf,
    pf_nb_kernel331nf,
    pf_nb_kernel332nf,
    pf_nb_kernel333nf,
    pf_nb_kernel334nf,
    pf_nb_kernel400nf,
    pf_nb_kernel410nf,
    pf_nb_kernel430nf,
};


void
pf_nb_kernel_setup(FILE *log,pf_nb_kernel_t **list)
{
  int i;
  pf_nb_kernel_t *p;

  if(NULL != log)
    fprintf(log,"Configuring pairwise force C nonbonded kernels...\n");

  for(i=0;i<eNR_NBKERNEL_NR;i++)
  {
    p = pf_kernellist[i];
    if(p!=NULL)
      list[i] = p;
  }
}    
