/*
 * Top level functions to add interactions.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_array_h
#define pf_array_h

#include "types/pf_array.h"
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_array_scalar.h"

void pf_atom_add_bonded(t_pf_global *pf_global, atom_id i, atom_id j, int type, rvec force);
void pf_atom_add_nonbonded_single(t_pf_global *pf_global, atom_id i, atom_id j, int type, real force, real dx, real dy, real dz);
void pf_atom_add_nonbonded(t_pf_global *pf_global, atom_id i, atom_id j, real pf_coul, real pf_lj, real dx, real dy, real dz);

#endif  /* pf_array_h */
