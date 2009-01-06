/**
 * \file modifications.h
 * \brief Datatypes and methods for amino acid modifications
 *
 * Two data structures define modifications.  The AA_MOD_T is the most
 * basic type.  It is the information provided by the user: mass
 * change caused by this mod, amino acids which may be modified in
 * this way, and the maximum number of this type of modification which
 * may occur on one peptide.  AA_MODs are defined here.
 * A collection of AA_MODs that may occur
 * on some peptide are represented as a PEPTIDE_MOD_T.  This stores
 * a list of AA_MODS and the net mass change experienced by the
 * peptide.  PEPTIDE_MODs are defined in peptide_modifications.h
 * AA_MODs are instantiated once after parsing the parameter file.  All
 * possible PEPTIDE_MODs are calcualted once and reused for each
 * spectrum search.  One PEPTIDE_MOD corresponds to one mass window
 * that must be searched.
 * 
 * $Revision: 1.2 $
 */
#ifndef MODIFICATION_FILE_H
#define MODIFICATION_FILE_H

#include <assert.h>
#include "utils.h"
#include "linked_list.h"
#include "objects.h"
#include "parameter.h"

/* Public constants */
enum {MAX_AA_MODS = 11};
enum {MAX_PROTEIN_SEQ_LENGTH = 40000};
enum {AA_LIST_LENGTH = 26}; // A-Z
#define GET_AA_MASK  0x001F   // 0000 0000 0001 1111
#define GET_MOD_MASK 0xFFE0   // 1111 1111 1110 0000


// this was moved to object.h b/c methods in peptide.h weren't compiling
//typedef unsigned short MODIFIED_AA_T; ///< letters in the expanded peptide
                                      ///alphabet, bits for mod1 mod2...aa
#define MOD_SEQ_NULL (MODIFIED_AA_T)('Z' - 'A' + 1) 
//enum {MOD_SEQ_NULL = 'Z' - 'A' + 1}; 
///< null terminating character of array

/*
   e.g. There are three aa_mods specified in this run and they are
   given the masks mod1  1000_0000_0000_0000
                   mod2  0100_0000_0000_0000
                   mod3  0010_0000_0000_0000
   The amino acid Y has the value
                      Y  0000_0000_0001_0011

   Suppose you set the variable aa to Y.  Do stuff.  To ask "is aa
  modified by mod2", do this. answer = mod2 AND aa; answer == 0
   To change aa to be modified by modcalle2, aa = mod2 NOR aa
 (i.e. !(mod2 || aa) )  now aa == 0100_0000_0001_0011
   If we ask again, answer == 0100_0000_0000_0000 
 */

/**
 * \brief Allocate an AA_MOD, including space for the aa_list and
 * initialize all fields to default values.  Symbol and unique
 * identifier are set according to index.  
 * \returns A heap allocated AA_MOD with default values.
 */
AA_MOD_T* new_aa_mod(int mod_idx);

/**
 * \brief Free the memory for an AA_MOD including the aa_list.
 */
void free_aa_mod(AA_MOD_T*);

/**
 * \brief Converts a MODIFIED_AA into a char, effectively unmodifying it.
 * \returns The unmodified char representation of an aa.
 */
char modified_aa_to_char(MODIFIED_AA_T aa);

/**
 * \brief Converts a char representation of an aa to an unmodified aa
 * of type MODIFIED_AA_T.  Requires 'A' <= char <= 'Z'. 
 * \returns The MODIFIED_AA_T represnetation of an aa.
 */
MODIFIED_AA_T char_aa_to_modified(char aa);

/**
 * \brief Converts a MODIFIED_AA_T* to it's textual representation,
 * i.e. a letter followed by between 0 and 11 symbols for the
 * modifications made to the amino acid.
 * \returns A newly allocated char* with amino acid and modifciation
 * symbols. 
 */
char* modified_aa_to_string(MODIFIED_AA_T aa);

/**
 * \brief Take an array of MODIFIED_AA_T's and return an array of
 * char's that includes the letter of each aa and the symbol for all
 * applied modifications.
 *
 * Assumes that the array is terminated with MOD_SEQ_NULL.
 * \returns A newly allocated array of characters, a text
 * representation of the modified sequence.
 */
char* modified_aa_string_to_string(MODIFIED_AA_T* aa_string, int length);

/**
 * \brief Allocates an array of MODIFIED_AA_T's the same length as
 * sequence and populates it with the MODIFIED_AA_T value that
 * corresponds to each sequence char value.  No modifications are
 * applied to the new array.
 *
 * \returns A newly allocated copy of the sequnce converted to type
 * MODIFIED_AA_T. 
 */
MODIFIED_AA_T* convert_to_mod_aa_seq(char* sequence);

/**
 * \brief Allocate a new MODIFIED_AA_T array and copy values into it.
 */
MODIFIED_AA_T* copy_mod_aa_seq(MODIFIED_AA_T* source, int length);

/**
 * \brief Frees memory for an array of MODIFIED_AA_Ts.  Assumes is
 * terminated with the MOD_SEQ_NULL value
 */
void free_mod_aa_seq(MODIFIED_AA_T* seq);

/**
 * \brief Gives the size of the aa_mod struct.  For serialization
 */
int get_aa_mod_sizeof();

/**
 * The new definition of a PEPTIDE_T object.
 * 
 */
/*
struct peptide{
  unsigned char length; ///< The length of the peptide
  float peptide_mass;   ///< The peptide's mass with any modifications
  PEPTIDE_SRC_T* peptide_src; ///< a linklist of peptide_src
  BOOLEAN_T is_modified;   ///< if true sequence != NULL
  MODIFIED_AA_T* sequence; ///< sequence with modifications
};
*/

/**
 * \brief checks to see if an amino acid is modified by a given mod
 * \returns TRUE if aa is modified by mod
 */
BOOLEAN_T is_aa_modified(MODIFIED_AA_T aa, AA_MOD_T* mod);

/**
 * \brief Determine if this modified amino acid can be modified by
 * this modification.
 *
 * Checks the mod list of modifiable residues to see if aa is in the
 * list.  Also checks to see if aa has already been modified by this
 * mod.  
 * \returns TRUE if it can be modified, else FALSE
 */
BOOLEAN_T is_aa_modifiable(MODIFIED_AA_T aa, AA_MOD_T* mod);

/**
 * \brief Adds a modification to a MODIFIED_AA_T.
 *
 * Assumes that the aa is modifiable, no explicit check.  If the aa is
 * already modified for the mod, no change to aa.
 */
void modify_aa(MODIFIED_AA_T* aa, AA_MOD_T* mod);

/**
 * \brief Finds the list of modifications made to an amino acid.
 *
 * Allocates a list of length(possible_mods) and fills it with pointers
 * to the modifications made to this aa as defined by
 * is_aa_modified().  Returns 0 and sets mod_list to NULL if the amino
 * acid is unmodified
 *
 * \returns The number of modifications made to this amino acid.
 */
int get_aa_mods(MODIFIED_AA_T aa, 
                AA_MOD_T* possible_mods, 
                AA_MOD_T** mod_list);


// in parameter.c
//global
//enum {MAX_AA_MODS = 11};//instead of #define forces typechecking, obeys scope
//AA_MOD_T* list_of_mods;
//int num_mods;
//AA_MOD_T* position_mods;
//int num_position_mods;

/**
 * \brief Read the paramter file and populate the static parameter
 * list of AA_MODS, inlcuding the list of position mods.
 *
 * Also updates the array of amino_masses.  Dies with an error if the
 * number of mods in the parameter file is greater than MAX_AA_MODS.
 * \returns void
 */
void read_mods_from_file(char* param_file);

/**
 * \brief Write the given aa mod to file in binary format.  Used for
 * serializing csm file headers.
 *
 * \return TRUE if written to file without error, else false.
 */
BOOLEAN_T serialize_aa_mod(AA_MOD_T* a_mod,
                           FILE* file
);

/**
 * \brief Read an aa mod from file in binary format as written by
 * serialize_aa_mod().  Overwrites any data in the passed aa_mod.
 * Used for reading in csm files. If FALSE is returned, the passed
 * a_mod will contain unpredictable values.
 *
 * \returns TRUE if aa_mod successfully read, else FALSE.  
 */
BOOLEAN_T parse_aa_mod(AA_MOD_T* a_mod,
                       FILE* file
);


/**
 * \brief Get the pointer to the list of AA_MODs requested by the
 * user.
 * \returns The number of mods pointed to by mods
 */
//int get_aa_mod_list(AA_MOD_T*** mods);

/**
 * \brief Check that the list of peptide_modifications from the file of
 * serialized PSMs matches those in the paramter file.
 *
 * If there was no parameter file or if it did not contain any mods,
 * return FALSE.  If the given mod list does not exactly match the
 * mods read from the parameter file (including the order in which
 * they are listed) return FALSE.  If returning false, print a warning
 * with the lines that should be included in the parameter file.
 *
 * \returns TRUE if the given mods are the same as those from the
 * parameter file.
 */
BOOLEAN_T compare_mods(AA_MOD_T** psm_file_mod_list, int num_mods);

/**
 * \brief Compare two mods to see if they are the same, i.e. same mass
 * change, unique identifier, position
 */
BOOLEAN_T compare_two_mods(AA_MOD_T* mod1, AA_MOD_T* mod2);

// in mass.c;
/**
 * \brief Extends the amino_masses table to include all possible
 * modifications.  
 *
 * Gets list of mods from parameter.c.  Should this fill in values for
 * both average and monoisotopic masses? 
 */
void extend_amino_masses(void);

/**
 * print all fields in mod.  For debugging
 */
void print_a_mod(AA_MOD_T* mod);

/* Setters and Getters */

/**
 * \brief Set the mass change caused by this modification.
 * \returns void
 */
void aa_mod_set_mass_change(AA_MOD_T* mod, double mass_change);
/**
 * \brief Get the mass change caused by this modification.
 * \returns The mass change caused by this modification.
 */
double aa_mod_get_mass_change(AA_MOD_T* mod);

/**
 * \brief Access to the aa_list of the AA_MOD_T struct.  This pointer
 * can be used to get or set the list of residues on which this mod
 * can be placed.
 * \returns A pointer to the list of amino acids on which this mod can
 * be placed.
 */
BOOLEAN_T* aa_mod_get_aa_list(AA_MOD_T* mod);

/**
 * \brief Set the maximum number of times this modification can be
 * placed on one peptide.
 * \returns void
 */
void aa_mod_set_max_per_peptide(AA_MOD_T* mod, int max);
/**
 * \brief Get the maximum number of times this modification can be
 * placed on one peptide.  
 * \returns The max times per peptide this mod can be placed.
 */
int aa_mod_get_max_per_peptide(AA_MOD_T* mod);

/**
 * \brief Set the maximum distance from the protein terminus that the
 * modification can be placed.  Which terminus (C or N) is determined
 * by the position type.  To indicate no position restriction, set to
 * MAX_PROTEIN_SEQ_LENGTH. 
 * \returns void
 */
void aa_mod_set_max_distance(AA_MOD_T* mod, int distance);
/**
 * \brief Get the maximum distance from the protein end that the
 * modification can be placed.  Will be MAX_PROTEIN_SEQ_LENGTH if
 * position type is ANY_POSITION.
 * \returns Maximum distance from protein terminus at which mod can be
 * placed. 
 */
int aa_mod_get_max_distance(AA_MOD_T* mod);

/**
 * \brief Set the position type of an aa_mod.
 * \returns void
 */
void aa_mod_set_position(AA_MOD_T* mod, MOD_POSITION_T position);

/**
 * \brief Where in the peptide can the modification be placed.
 * \returns ANY_POSITION for standard mods; C_TERM or N_TERM for those
 * that can only be placed at the ends of the peptide.
 */
MOD_POSITION_T aa_mod_get_position(AA_MOD_T* mod);

/**
 * \brief The character used to uniquely identify the mod in the sqt file.
 * \returns The character identifier.
 */
char aa_mod_get_symbol(AA_MOD_T* mod);

/**
 * \brief The bitmask used to uniquely identify the mod.
 * \returns The short int bitmask used to identify the mod.
 */
int aa_mod_get_identifier(AA_MOD_T* mod);

/**
 * \brief Generates a string representation of an aa_mod and returns a
 * pointer to that newly allocated string.
 */
char* aa_mod_to_string(AA_MOD_T* mod);

/**
 * \brief Create a string containing all of the amino acids that can
 * be modified by this aa_mod.  E.g. if S, T, and Y can be modified,
 * returns "STY".
 * \returns A newly allocated string.
 */
char* aa_mod_get_aa_list_string(AA_MOD_T* mod);

#endif //MODIFICATION_FILE_H
