/*
 * Functions which act on types/pf_array structures.
 *
 * Copyright Bogdan Costescu 2010-2012
 */

#ifndef pf_utils_h
#define pf_utils_h

#include "types/pf_array.h"

t_pf_atom_id_list *pf_atom_id_list_alloc(int len);
void pf_atom_id_list_free(t_pf_atom_id_list *p);
void pf_atoms_scalar_real_divide(t_pf_atoms *atoms, real divisor);
void pf_atom_summed_merge_to_scalar(t_pf_atom_summed *src, t_pf_atom_scalar *dst, const rvec *x, int Vector2Scalar);
void pf_atoms_summed_merge_to_scalar(t_pf_atoms *atoms, const rvec *x, int Vector2Scalar);
real pf_vector2signedscalar(const rvec v, const rvec xi, const rvec xj, int Vector2Scalar);
atom_id *pf_init_sys2pf(int syslen);
void pf_fill_sys2pf(atom_id *sys2pf, atom_id *len, t_pf_atom_id_list *p);
char *pf_make_sys_in_group(int syslen, t_pf_atom_id_list *p);
void pf_fill_atom2residue(t_pf_global *pf_global, gmx_mtop_t *top_global);
t_pf_atom_id_list *pf_group2atoms(int len, atom_id *list);
t_pf_atom_id_list *pf_groupatoms2residues(t_pf_atom_id_list *atoms, t_pf_global *pf_global);
void pf_read_group(t_pf_global *pf_global, const char *ndxfile, char *groupname, char **sys_in_g);
void pf_check_sys_in_g(t_pf_global *pf_global);
int pf_interactions_type_str2val(char *typestr);
char *pf_interactions_type_val2str(int type);
void pf_atoms_alloc(int OnePair, t_pf_atoms *atoms, int syslen, char *name);
void pf_atoms_init(int OnePair, t_pf_atoms *atoms);
void pf_atoms_scalar_alloc(t_pf_atoms *atoms, int syslen, char *name);
void pf_atoms_scalar_init(t_pf_atoms *atoms);
void pf_atoms_and_residues_init(t_pf_global *pf_global);
void pf_output_data_per_frame_init(t_pf_global *pf_global);
t_pf_global *pf_init(FILE *fplog, int nfile, const t_filenm fnm[], gmx_mtop_t *top_global);

#endif  /* pf_utils_h */
