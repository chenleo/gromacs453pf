/*
 * Functions for per-atom interactions.
 *
 * Copyright Bogdan Costescu 2011-2013
 */

void pf_per_atom_real_set(t_pf_per_atom_real *per_atom_real, real val);
void pf_per_atom_real_init(t_pf_per_atom_real **per_atom_real, atom_id len, real val);
void pf_per_atom_real_int_set(t_pf_per_atom_real_int *per_atom_real_int, real val_real, int val_int);
void pf_per_atom_real_int_init(t_pf_per_atom_real_int **per_atom_real_int, atom_id len, real val_real, int val_int);
void pf_per_atom_real_write_frame(FILE *f, real *force, atom_id len, gmx_bool no_end_zeros);

void pf_per_atom_sum(t_pf_per_atom_real *per_atom_real, t_pf_atom_summed *atoms, atom_id atoms_len, const rvec *x, int Vector2Scalar);
void pf_per_atom_average(t_pf_per_atom_real_int *per_atom_average, t_pf_atom_summed *atoms, atom_id atoms_len, const rvec *x, int Vector2Scalar);
void pf_per_atom_minmax(t_pf_per_atom_real *per_atom_real, t_pf_atom_summed *atoms, atom_id atoms_len, gmx_bool findmax, const rvec *x, int Vector2Scalar);
