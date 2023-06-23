// Benjamin Diament

#include <iostream>
#include <limits>
#include <gflags/gflags.h>
#include "mass_constants.h"
#include "max_mz.h"
#include "fifo_alloc.h"
#include "fixed_cap_array.h"
#include "theoretical_peak_set.h"
#include "peptide_lite.h"
#include "compiler.h"
#include "util/StringUtils.h"

PeptideLite::PeptideLite(const pb::Peptide& peptide,
        const vector<const pb::Protein*>& proteins,
        vector<const pb::AuxLocation*>* locations)
  : len_(peptide.length()), 
  mass_(peptide.mass()), 
  id_(peptide.id()),
  mods_(NULL), 
  num_mods_(0), 
  proteins_(&proteins),
  first_loc_protein_id_(peptide.first_location().protein_id()),
  first_loc_pos_(peptide.first_location().pos()), 
  protein_length_(proteins[first_loc_protein_id_]->residues().length()),
  decoyIdx_(peptide.has_decoy_index() ? peptide.decoy_index() : -1) {
    
  // Here we make sure that tide-search is compatible with old and new tide-index protocol buffers.
  // Set residues_ by pointing to the first occurrence in proteins.
  if (peptide.has_decoy_sequence() == true){  //new tide-index format
    decoy_seq_ = peptide.decoy_sequence();  // Make a copy of the string, because pb::Peptide will be reused.
    residues_ = decoy_seq_.data();
    target_residues_ = proteins[first_loc_protein_id_]->residues().data() 
                      + first_loc_pos_;
  } else {  //old tide-index format
    residues_ = proteins[first_loc_protein_id_]->residues().data() 
                    + first_loc_pos_;
    if (IsDecoy()) {
      target_residues_ = proteins[first_loc_protein_id_]->residues().data() 
                        + first_loc_pos_+len_+1;
    } else {
      target_residues_ = residues_;
    }
  }
                    
  num_mods_ = peptide.modifications_size();
  for (int i = 0; i < num_mods_; ++i)
    mods_.push_back(ModCoder::Mod(peptide.modifications(i)));

  // Add auxiliary locations;
  // This handles the old version of tide index in which the aux locations are separately stored in a vector.    
  if (peptide.aux_locations_index() && locations) {  
    const pb::AuxLocation* aux_loc = locations->at(peptide.aux_locations_index());
    for (int i = 0; i < aux_loc->location_size(); ++i) {
      aux_locations.push_back(aux_loc->location(i));
    }
  // this part handles the aux location in the current tide-index.     
  // Here the aux location are merged to the peptide pb.
  } else if (peptide.has_aux_loc() == true) {   
    const pb::AuxLocation& aux_loc = peptide.aux_loc();
    for (int i = 0; i < aux_loc.location_size(); ++i) {
      aux_locations.push_back(aux_loc.location(i));
    }
  }
  b_ions_only_ = false;
  string scoring_method = Params::GetString("score-function"); // Handle this properly with Get
  if ( scoring_method == "combined-p-values") 
    b_ions_only_ = true; 
  protein_id_str_ = string("");
  flankingAAs_ = string("");
  seq_with_mods_ = string("");
  mod_string_ = string("");
}

template<class W>
void PeptideLite::AddIons(W* workspace, bool dia_mode) {
  // Use workspace to assemble all B and Y ions. workspace will determine
  // which, if any, associated ions will be represented.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0)
    max_possible_peak = MaxBin::Global().CacheBinEnd();

  vector<double> aa_masses = getAAMasses();

  // added by Yang
  if (dia_mode) {
    // ion_mzs_.clear(); ion_mzbins_.clear(); b_ion_mzbins_.clear(); y_ion_mzbins_.clear();
    ion_mzbins_.clear(); b_ion_mzbins_.clear(); y_ion_mzbins_.clear();
  }

  // Add all charge 1 B ions.
  double total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 1);
    if (dia_mode) {
      b_ion_mzbins_.push_back(MassConstants::mass2bin(PeptideLite::MassToMz(total + MassConstants::B, 1)));
      ion_mzbins_.push_back(MassConstants::mass2bin(PeptideLite::MassToMz(total + MassConstants::B, 1)));
      // ion_mzs_.push_back(PeptideLite::MassToMz(total + MassConstants::B, 1));
    }
    total += aa_masses[i];
  }
  if (b_ions_only_ == false) {
    // Add all charge 1 Y ions.
    total = aa_masses[Len() - 1];
    for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
      workspace->AddYIon(total, 1);
      if (dia_mode) {
        y_ion_mzbins_.push_back(MassConstants::mass2bin(PeptideLite::MassToMz(total + MassConstants::Y, 1)));
        ion_mzbins_.push_back(MassConstants::mass2bin(PeptideLite::MassToMz(total + MassConstants::Y, 1)));
        // ion_mzs_.push_back(PeptideLite::MassToMz(total + MassConstants::Y, 1));
      }
      total += aa_masses[i];
    }

    // Add all charge 2 B ions.
    max_possible_peak = max_possible_peak*2 + 2;  //adjust for larger charge
    total = aa_masses[0];
    for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
      workspace->AddBIon(total, 2);
      total += aa_masses[i];
    }

    // Add all charge 2 Y ions.
    total = aa_masses[Len() - 1];
    for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
      workspace->AddYIon(total, 2);
      total += aa_masses[i];
    }
  }

  // added by Yang
  if (dia_mode) {
    sort(b_ion_mzbins_.begin(), b_ion_mzbins_.end());
    sort(y_ion_mzbins_.begin(), y_ion_mzbins_.end());
    sort(ion_mzbins_.begin(), ion_mzbins_.end());
    // sort(ion_mzs_.begin(), ion_mzs_.end());
  }
}

void PeptideLite::Compile(const TheoreticalPeakArr* peaks) {
  // Store the theoretical peak indeces for the peptide in a vector. 
  int i;
  peaks_0.clear();
  peaks_1.clear();
  peaks_0.reserve(peaks[0].size());
  peaks_1.reserve(peaks[1].size());
  
  for (i = 0; i < peaks[0].size(); ++i) {
    peaks_0.push_back(peaks[0][i]);
  }
  for (i = 0; i < peaks[1].size(); ++i) {
    peaks_1.push_back(peaks[1][i]);
  }
}

void PeptideLite::ComputeTheoreticalPeaks(TheoreticalPeakSetBYSparse* workspace, bool dia_mode) {
  AddIons<TheoreticalPeakSetBYSparse>(workspace, dia_mode);   // Generic workspace
  Compile(workspace->GetPeaks());
}

// return the amino acid masses in the current peptide
vector<double> PeptideLite::getAAMasses() const {
  vector<double> masses_charge(Len());
  const char* residue = residues_;
  for (int i = 0; i < Len(); ++i, ++residue) {
    if (i == 0) { // nterm static pep
      if(first_loc_pos_ == 0)
        masses_charge[i] = MassConstants::nprotterm_mono_table[*residue];
      else
        masses_charge[i] = MassConstants::nterm_mono_table[*residue];
    } else if (i == Len() - 1) { // cterm static pep
      if(first_loc_pos_ + len_ == protein_length_ - 1)
        masses_charge[i] = MassConstants::cprotterm_mono_table[*residue];
      else
        masses_charge[i] = MassConstants::cterm_mono_table[*residue];
    } else { // all other mods
      masses_charge[i] = MassConstants::mono_table[*residue];
    }
  }

  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    masses_charge[index] += delta;
  }
  return masses_charge;
}

string PeptideLite::SeqWithMods(int mod_precision) {
  // If the peptide is reported more than once then reuse the strings from previous calculations    
  if (seq_with_mods_.empty() == false)
    return seq_with_mods_;
  
  seq_with_mods_ = string(residues_, Len());  // Get the plain peptide sequence
  
  int mod_pos_offset = 0;
  string mod_str;
  sort(mods_.begin(), mods_.end());
  
  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    mod_str = '[' + StringUtils::ToString(delta, mod_precision) + ']';
    seq_with_mods_.insert(index + 1 + mod_pos_offset, mod_str);
    mod_pos_offset += mod_str.length();
  }
  return seq_with_mods_;  
}

/**
 * Gets the protein name with the peptide position appended. For reporting results
 */
 string  PeptideLite::GetLocationStr(const string& decoy_prefix) {
  // If the peptide is reported more than once then reuse the strings from previous calculations  

  if (protein_id_str_.empty() == false) 
    return protein_id_str_;

  string locations;
  locations = (IsDecoy()?decoy_prefix:"") + 
    proteins_->at(FirstLocProteinId())->name() + 
    "(" + std::to_string(FirstLocPos()+1) + ")";
  
  for (vector<pb::Location>::const_iterator 
    loc = aux_locations.begin(); 
    loc != aux_locations.end();
    ++loc
  ) {
    const pb::Protein* protein = proteins_->at((*loc).protein_id());
    int pos = (*loc).pos()+1;
    locations += "," + (IsDecoy()?decoy_prefix:"") + protein->name() + 
      "(" + std::to_string(pos) + ")";    
  }
  protein_id_str_ = locations;
  return locations;
}

/**
 * Gets the flanking AAs for a Tide peptide sequence for reporting results
 */
string  PeptideLite::GetFlankingAAs() {
  // If the peptide is reported more than once then reuse the strings from previous calculations  

  if (flankingAAs_.empty() == false) 
    return flankingAAs_;

  string flankingAAs;
  flankingAAs.clear();
  int prot_pos = FirstLocPos();
  const string& seq = proteins_->at(FirstLocProteinId())->residues();

  flankingAAs = ((prot_pos > 0) ? proteins_->at(FirstLocProteinId())->residues().substr(prot_pos-1, 1) : "-") +
    ((prot_pos+Len() <  proteins_->at(FirstLocProteinId())->residues().length()) ? 
    proteins_->at(FirstLocProteinId())->residues().substr(prot_pos+Len(),1) : "-");
    
  for (vector<pb::Location>::const_iterator
    loc = aux_locations.begin();
    loc != aux_locations.end();
    ++loc
  ) {
    const pb::Protein* protein = proteins_->at((*loc).protein_id());
    prot_pos = (*loc).pos();
    flankingAAs += "," + ((prot_pos > 0) ? protein->residues().substr(prot_pos-1, 1) : "-") + 
      ((prot_pos+Len() <  protein->residues().length()) ? 
      protein->residues().substr(prot_pos+Len(),1) : "-");
  }
  flankingAAs_ = flankingAAs;
  return flankingAAs;

}

// const pb::ModTable* MassConstants::mod_table_ = NULL; 
// const pb::ModTable* MassConstants::n_mod_table_ = NULL; 
// const pb::ModTable* MassConstants::c_mod_table_ = NULL; 
// const pb::ModTable* MassConstants::nprot_mod_table_ = NULL; 
// const pb::ModTable* MassConstants::cprot_mod_table_ = NULL;  

string PeptideLite::getModifications(int mod_precision) {
  // If the peptide is reported more than once then reuse the strings from previous calculations  
  if (mod_string_.empty() == false)  
    return mod_string_;

  if (seq_with_mods_.empty() == true)
    sort(mods_.begin(), mods_.end());

  vector<string> mods_list;
  string sep("_");
  int prot_pos = FirstLocPos();  
  double mod_mass;
  int var_mod_idx = 0;
  for (int i = 0; i < len_; ++i) {
    char AA = residues_[i];
    // Find static modifications 

    if (i == 0) {// peptide N-term
      if (find_static_mod(MassConstants::n_mod_table_, AA, mod_mass)) {
        //Get the string
        string mod = std::to_string(i+1) + sep + string("S") +  sep + StringUtils::ToString(mod_mass, mod_precision) + "_n";   
        mods_list.push_back(mod);
      }
    }
    if (i == 0 && prot_pos == 0) {// protein N-term
      if (find_static_mod(MassConstants::nprot_mod_table_, AA, mod_mass)) {
        //Get the string
        string mod = std::to_string(i+1) + sep + string("S") +  sep + StringUtils::ToString(mod_mass, mod_precision) + "_N";   
        mods_list.push_back(mod);
      }
    }

    if (i == len_-1) {  // peptide C-term
      if (find_static_mod(MassConstants::c_mod_table_, AA, mod_mass)) {
        //Get the string
        string mod = std::to_string(i+1) + sep + string("S") +  sep + StringUtils::ToString(mod_mass, mod_precision) + "_c";   
        mods_list.push_back(mod);
      }
    }
    if (i == len_-1 && prot_pos +i + 1 ==  proteins_->at(FirstLocProteinId())->residues().length() ) {  // protein C-term
      if (find_static_mod(MassConstants::cprot_mod_table_, AA, mod_mass)) {
        //Get the string
        string mod = std::to_string(i+1) + sep + string("S") +  sep + StringUtils::ToString(mod_mass, mod_precision) + "_C";   
        mods_list.push_back(mod);
      }
    }
    if (find_static_mod(MassConstants::mod_table_, AA, mod_mass)) { // Static mod
      //Get the string
      string mod = std::to_string(i+1) + sep + string("S") +  sep + StringUtils::ToString(mod_mass, mod_precision);
      mods_list.push_back(mod);
    }

    // Add variable mods
    if (var_mod_idx < mods_.size()) {
      int index;
      double delta;
      MassConstants::DecodeMod(mods_[var_mod_idx], &index, &delta);
      if (index == i) {
        string mod = std::to_string(i+1) + sep + string("V") +  sep + StringUtils::ToString(delta, mod_precision);
        mods_list.push_back(mod);
        var_mod_idx++;
      }
    }
  }
  mod_string_ = StringUtils::Join(mods_list, ',');
  return mod_string_;
}

bool PeptideLite::find_static_mod(const pb::ModTable* mod_table, char AA, double& mod_mass) {

  for (int i = 0; i < mod_table->static_mod_size(); i++) {
    
    const pb::Modification& mod = mod_table->static_mod(i);

    if (mod.has_delta() && mod.has_amino_acids()) {

      string AAs = mod.amino_acids();
      int AA_len = AAs.length();
      for (int j = 0; j < AA_len; ++j) {
        if (AAs[j] == AA || AAs[j] == 'X') { // Found a static mod for Amino acid AA;
          mod_mass =  mod.delta();
          return true;
        }
      }
    }
  }
  return false;  //No modification
}
