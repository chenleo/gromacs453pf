/*
 * Arrays of pairwise forces. This is a compromise between memory required to
 * store the pairwise forces and the CPU time needed to access them (insert,
 * lookup).
 * For each atom i there is a struct containing several arrays; there are
 * separate arrays for each interaction type, so there is no need to keep an
 * extra 'type' field. Each interaction takes space for one atom index and a
 * real vector containing the value of the force.
 * The arrays are grown as needed, through memory reallocations, similar to the
 * nonbonded lists; to avoid too many reallocations (memory management can be
 * expensive especially with certain MPI implementations which register memory
 * to achieve zero copy), they are always done by increasing the current size
 * with a factor > 1. However, this factor should be kept small enough to
 * avoid allocating too much memory which will then remain unused.
 * For each array there is a counter for the number of items currently in use
 * and one for the number of items for which space was allocated. If they are
 * equal and there is a request for storing another interaction, a reallocation
 * will happen.
 * Lookup is done through a linear search in the array. This uses CPU time, but
 * saves memory - organizing the items in a binary tree would require 2
 * pointers per interaction and allocations cannot be done in bulk anymore.
 * Insertion can be done directly at the end if it's known that the jjnr does
 * not already exist in the list; otherwise a linear search is needed.
 *
 * For the case of summed up forces per pair, a different set of data structures
 * are needed; the main reason is the extra data present in the summed up
 * interaction structure (the type) - it makes no sense to use a single structure
 * and ignore this field for the detailed storage of interactions, as this will
 * increase significantly the storage requirements. Furthermore, when storing
 * the atoms in group1, it makes no sense to keep all those NULL pointers
 * towards lists of interactions.
 * 
 * NOTE: all structures are defined for atoms, but are then used for both atoms and residues;
 * it made no sense to define new ones as they are functionally the same
 * 
 * Copyright Bogdan Costescu 2010-2013
 */

#ifndef pf_types_array_h
#define pf_types_array_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdio.h>

#include "simple.h"
#include "pf_array_detailed.h"
#include "pf_array_summed.h"
#include "pf_array_scalar.h"
#include "pf_per_atom.h"

/* type for a list of atom_id's, containing both the list itself and the
 * length of the list
 */
typedef struct {
  atom_id *list;
  int len;
} t_pf_atom_id_list;

/* TODO: the bonded interaction type  might be possible to be reduced to a 
 * single real vector (=eliminate jjnr) as the atoms interacting are known
 * at simulation start and do not change (also their order in the bond list
 * doesn't change). This would however require an extra step of looking up
 * the atom indeces when writing out the force matrix. This would also require
 * a change in the way the atom structure is accessed, it makes no sense to
 * keep an array of t_pf_interaction items, but only an array of reals 
 * representing forces.
 */
/* type to represent all atoms */

typedef struct {
  atom_id *sys2pf;		/* indexing table: real atom nr. to index in the pf array; this has length equal to the total nr. of atoms in system */
  atom_id len;			/* nr. of atoms from both groups, gives the length of the atoms lists below */
  /* only one of atoms and atom_summed lists will be initialized */
  t_pf_atom_detailed *detailed;	/* list of atoms */
  t_pf_atom_summed *summed;	/* list of atoms */
  t_pf_atom_scalar *scalar;     /* list of atoms */
} t_pf_atoms;

typedef struct {
  int period;       /* nr. of steps to average before writing; if 1 (default), no averaging is done; if 0, averaging is done over all steps so only one frame is written at the end */
  int steps;        /* counter for current step, incremented for every call of pf_save_and_write_scalar_averages(); when it reaches time_averages_steps, data is written */
  rvec *com;                    /* averaged residue COM coordinates over steps, needed for COM calculations; only initialized when pf_global->ResidueBased is non-zero */
} t_pf_time_averages;

typedef struct {
  gmx_bool bInitialized;        /* TRUE if pairwise forces should be written out, FALSE otherwise; if FALSE, many of the following structure members will not be initialized */
  /* AtomBased & ResidueBased below are not boolean, to allow expressing preferences not only for storing in memory but also about the format they are written in:
   * value | value in fi file | meaning
   * 0     | no               | no storing (default)
   * 1     | scalar           | store as vectors, write as scalars
   * 2     | vector           | store as vectors, write as vectors
   * 3     | compat_bin       | store as vectors, write in compatibility mode (signed scalars) in binary
   * 4     | compat_ascii     | store as vectors, write in compatibility mode (signed scalars) in ascii
   * Also see the constants defined below.
   */
  int AtomBased;                /* independent of ResidueBased */
  int ResidueBased;             /* independent of AtomBased */
  /* OnePair defines the way the interactions between the same pair of atoms are stored:
   * value | value in fi file | meaning
   * 0     | detailed         | each interaction is stored separately; it's possible to have the same pair in several interaction lists (default)
   * 1     | summed           | each interaction is stored once
   */
  int OnePair;
  /* Vector2Scalar defines the way the force vector to scalar conversion is done:
   * value | value in fi file | meaning
   * 0     | norm             | takes the norm of the vector
   * 1     | projection       | takes the projection of the force on the direction of the 2 atoms
   * As the conversion affects all scalar writing modes (PF_FILE_OUT*), this is kept as a separate
   * setting rather than creating separates modes for norm and projection.
   */
  int Vector2Scalar;
  int syslen_atoms;             /* total nr. of atoms in the system; this is a local copy to avoid passing too many variables down the function call stack */
  int syslen_residues;          /* maximum of residue nr. + 1; residue nr. doesn't have to be continuous, there can be gaps */
  char *sys_in_g1;              /* 0 if atom not in group1, 1 if atom in group1, length of syslen_atoms; always allocated */
  char *sys_in_g2;              /* 0 if atom not in group2, 1 if atom in group2, length of syslen_atoms; always allocated */
  t_pf_atoms *atoms;            /* atoms, always allocated but structure only initialized if AtomBased is non-zero */
  t_pf_atoms *residues;         /* residues, always allocated but structure only initialized if ResidueBased is non-zero */
  atom_id *atom2residue;        /* stores the residue nr. for each atom; array of length syslen; only initialized if ResidueBased is non-zero */
  int ResiduesRenumber;         /* detect/force residue renumbering */
  int type;                     /* interaction types that are interesting, set based on input file; functions are supposed to test against this before calculating/storing data */
  char *ofn_atoms;              /* output file name for atoms if AtomBased is non-zero */
  FILE *of_atoms;               /* output file for atoms if AtomBased is non-zero */
  char *ofn_residues;           /* output file name for residues if ResidueBased is non-zero */
  FILE *of_residues;            /* outpuf file for residues if ResidueBased is non-zero */
  char *groupname;              /* name of group for output in compatibility mode */
  /* the following 2 should be gmx_large_int_t, but the file format is defined with int... */
  int nsteps_atoms;             /* nr. of steps for output in compatibility mode; incremented for each step written during run, also used to write the total nr. of steps at the end */
  int nsteps_residues;          /* nr. of steps for output in compatibility mode; incremented for each step written during run, also used to write the total nr. of steps at the end */
  t_pf_time_averages *time_averages;    /* only initialized when time averages are calculated */
  t_pf_per_atom_real *per_atom_real;            /* only initialized if required by user */
  t_pf_per_atom_real_int *per_atom_real_int;    /* only initialized if required by user */
  t_pf_per_atom_real *per_residue_real;          /* only initialized if required by user */
  t_pf_per_atom_real_int *per_residue_real_int; /* only initialized if required by user */
  gmx_bool no_end_zeros;	/* if True, trim the line such that the zeros at the end are not written; if False (default), all per atom/residue data is written */
} t_pf_global;

/* values for OnePair */
enum {
  PF_ONEPAIR_DETAILED,
  PF_ONEPAIR_SUMMED,
  PF_ONEPAIR_NR
};

/* values for AtomBased and ResidueBased; some have ATOM in their name, but are applied to residues as well */
enum {
  PF_FILE_OUT_NONE,		/* this should be always kept first (=zero), as the code does if (pf_global->AtomBased) */
  PF_FILE_OUT_SCALAR,
  PF_FILE_OUT_VECTOR,
  PF_FILE_OUT_COMPAT_BIN,	/* compat is always scalar */
  PF_FILE_OUT_COMPAT_ASCII,	/* compat is always scalar */
  PF_FILE_OUT_PER_ATOM_SUM,	/* scalar */
  PF_FILE_OUT_PER_ATOM_AVERAGE,	/* scalar */
  PF_FILE_OUT_PER_ATOM_MIN,	/* scalar */
  PF_FILE_OUT_PER_ATOM_MAX,	/* scalar */
  PF_FILE_OUT_NR
};

/* values for Vector2Scalar */
enum {
  PF_VECTOR2SCALAR_NORM,
  PF_VECTOR2SCALAR_PROJECTION,
  PF_VECTOR2SCALAR_NR
};

/* values for ResidueRenumber */
enum {
  PF_RESIDUE_RENUMBER_AUTO,
  PF_RESIDUE_RENUMBER_DO,
  PF_RESIDUE_RENUMBER_DONT,
  PF_RESIDUE_RENUMBER_NR
};

#define pf_in_compatibility_mode(x) ((x == PF_FILE_OUT_COMPAT_BIN) || (x == PF_FILE_OUT_COMPAT_ASCII))
#define pf_file_out_per_atom(x) ((x == PF_FILE_OUT_PER_ATOM_SUM) || (x == PF_FILE_OUT_PER_ATOM_AVERAGE) || (x == PF_FILE_OUT_PER_ATOM_MIN) || (x == PF_FILE_OUT_PER_ATOM_MAX))

/* the initial number of items allocated for an array */
#define PF_ARRAY_INITSIZE	8

#endif  /* pf_types_array_h */
