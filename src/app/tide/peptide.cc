// Benjamin Diament

#include <iostream>
#include <limits>
#include <gflags/gflags.h>
#include "mass_constants.h"
#include "max_mz.h"
#include "fifo_alloc.h"
#include "fixed_cap_array.h"
#include "theoretical_peak_set.h"
#include "peptide.h"
#include "compiler.h"
#include "util/StringUtils.h"

#ifdef DEBUG
DEFINE_int32(debug_peptide_id, -1, "Peptide id to debug.");
#endif

#if 0
DEFINE_bool(flanks, true, "Include flanking peaks.");
DEFINE_bool(dups_ok, false, "Don't remove duplicate peaks");
#endif

string Peptide::SeqWithMods() const {
  vector<char> buf(Len() + num_mods_ * 30 + 1);
  int residue_pos = 0;
  char* buf_pos = &(buf[0]);
  for (int i = 0; i < num_mods_; ++i) {
    int index;
    double delta;
    MassConstants::DecodeMod(mods_[i], &index, &delta);
    while (residue_pos <= index)
      *buf_pos++ = residues_[residue_pos++];
    buf_pos += sprintf(buf_pos, "[%s%.1f]", delta >= 0 ? "+" : "", delta);
  }
  while (residue_pos < Len())
    *buf_pos++ = residues_[residue_pos++];
  *buf_pos = '\0';
  return &(buf[0]);
}

void Peptide::Show() {
#ifdef DEBUG
  if (Id() == FLAGS_debug_peptide_id) {
    cout << Seq() << endl;
    cout << "Charge 1 Pos" << endl;
    for (int i = 0; i < peaks_charge_1_.size(); ++i)
      cout << "Theoretical Peak[" << peaks_charge_1_[i].Bin() << "] = "
           << peaks_charge_1_[i].Type() << endl;
    cout << "Charge 1 Neg" << endl;
    for (int i = 0; i < negs_charge_1_.size(); ++i)
      cout << "Theoretical Peak[" << negs_charge_1_[i].Bin() << "] = "
           << negs_charge_1_[i].Type() << endl;
    cout << "Charge 2 Pos" << endl;
    for (int i = 0; i < peaks_charge_2_.size(); ++i)
      cout << "Theoretical Peak[" << peaks_charge_2_[i].Bin() << "] = "
             << peaks_charge_2_[i].Type() << endl;
    cout << "Charge 2 Neg" << endl;
    for (int i = 0; i < negs_charge_2_.size(); ++i)
      cout << "Theoretical Peak[" << negs_charge_2_[i].Bin() << "] = "
           << negs_charge_2_[i].Type() << endl;
    cout << endl;
  }
#endif
}

template<class W>
void Peptide::AddIons(W* workspace) {
  // Use workspace to assemble all B and Y ions. workspace will determine
  // which, if any, associated ions will be represented.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0)
    max_possible_peak = MaxBin::Global().CacheBinEnd();

  vector<double> aa_masses = getAAMasses();
  // carp(CARP_DETAILED_DEBUG, "**********aa_masses:%s", StringUtils::JoinDoubleVec(aa_masses, ',').c_str() );


  // added by Yang
  ion_mzbins_.clear();
  ion_mzs_.clear();

  // Add all charge 1 B ions.
  double total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 1);
    ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::B, 1)));
    ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::B, 1));
    total += aa_masses[i];
  }

  // Add all charge 1 Y ions.
  total = aa_masses[Len() - 1];
  for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
    workspace->AddYIon(total, 1);
    ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::Y, 1)));
    ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::Y, 1));
    total += aa_masses[i];
  }

  // Add all charge 2 B ions.
  max_possible_peak = max_possible_peak*2 + 2;  //adjust for larger charge
  total = aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace->AddBIon(total, 2);
    ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::B, 2)));
    ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::B, 2));
    total += aa_masses[i];
  }

  // Add all charge 2 Y ions.
  total = aa_masses[Len() - 1];
  for (int i = Len()-2; i >= 0 && total <= max_possible_peak; --i) {
    workspace->AddYIon(total, 2);
    ion_mzbins_.push_back(MassConstants::mass2bin(Peptide::MassToMz(total + MassConstants::Y, 2)));
    ion_mzs_.push_back(Peptide::MassToMz(total + MassConstants::Y, 2));
    total += aa_masses[i];
  }

  // added by Yang
  // sort(ion_mzbins_.begin(), ion_mzbins_.end());
  // sort(ion_mzs_.begin(), ion_mzs_.end());

  /*carp(CARP_DETAILED_DEBUG, "peptide:%s", Seq().c_str() );
  for (int ion_idx=0; ion_idx<ion_mzbins_.size(); ++ion_idx) {
	  carp(CARP_DETAILED_DEBUG, "ion mzbin:%d", ion_mzbins_[ion_idx] );
  }*/

}

template<class W>
void Peptide::AddBIonsOnly(W* workspace) const {
  // Use workspace to assemble b ions only.
  // Intended primarily to support XCorr p-value calculations.
  double max_possible_peak = numeric_limits<double>::infinity();
  if (MaxBin::Global().MaxBinEnd() > 0) {
    max_possible_peak = MaxBin::Global().CacheBinEnd();
  }
  
  vector<double> aa_masses = getAAMasses();

  // Add all charge 1 B ions.
  double total = MASS_PROTON + aa_masses[0];
  for (int i = 1; i < Len() && total <= max_possible_peak; ++i) {
    workspace -> AddBIon(total);
    total += aa_masses[i];
  }
}

#ifdef DEBUG
void DisAsm(const void* prog) {
  unsigned char* pos = (unsigned char*) prog;
  while (true) {
    switch(*pos) {
    case 3: assert(pos[1] == 130); cout << "add " << *((int*) (void*) (pos+2)) << "(%edx), %eax\n"; pos += 6; break;
    case 43: assert(pos[1] == 130); cout << "sub " << *((int*) (void*) (pos+2)) << "(%edx), %eax\n"; pos += 6; break;
    case 195: cout << "ret" << endl; goto out;
    default: assert(false);
    }
  }
 out:
  return;
}
#endif

void Peptide::Compile(const TheoreticalPeakArr* peaks,
                      const pb::Peptide& pb_peptide,
                      TheoreticalPeakCompiler* compiler_prog1,
                      TheoreticalPeakCompiler* compiler_prog2) {
  int pos_size = peaks[0].size();
  prog1_ = compiler_prog1->Init(pos_size, 0);
  compiler_prog1->AddPositive(peaks[0]);
//  compiler_prog1->AddPositive(pb_peptide.peak1());
//  compiler_prog1->AddNegative(pb_peptide.neg_peak1());
  compiler_prog1->Done();

  pos_size = peaks[0].size() + peaks[1].size();
  prog2_ = compiler_prog2->Init(pos_size, 0);
  compiler_prog2->AddPositive(peaks[0]);
  compiler_prog2->AddPositive(peaks[1]);
//  compiler_prog2->AddPositive(pb_peptide.peak2());
//  compiler_prog2->AddNegative(pb_peptide.neg_peak2());
  compiler_prog2->Done();
/*    cout << Seq() << endl;
    for (int i = 0; i < peaks[0].size(); ++i)
      cout << "Theoretical Peak[" << peaks[0][i].Bin() << "] = "
           << peaks[0][i].Type() << endl;
    for (int i = 0; i < peaks[1].size(); ++i)
      cout << "Theoretical Peak[" << peaks[1][i].Bin() << "] = "
           << peaks[1][i].Type() << endl;
*/
//	exit(1);  
}

void Peptide::ComputeTheoreticalPeaks(TheoreticalPeakSet* workspace) {
  AddIons<TheoreticalPeakSet>(workspace);   // Generic workspace
#ifdef DEBUG
  Show();
#endif
}

void Peptide::ComputeBTheoreticalPeaks(TheoreticalPeakSetBIons* workspace) const {
  AddBIonsOnly<TheoreticalPeakSetBIons>(workspace);   // workspace for b ion only peak set
#ifdef DEBUG
  Show();
#endif
}

void Peptide::ComputeTheoreticalPeaks(ST_TheoreticalPeakSet* workspace,
                                      const pb::Peptide& pb_peptide,
                                      TheoreticalPeakCompiler* compiler_prog1,
                                      TheoreticalPeakCompiler* compiler_prog2) {
  // Search-time fast workspace
  AddIons<ST_TheoreticalPeakSet>(workspace);

#if 0
  TheoreticalPeakArr peaks[2];
  peaks[0].Init(2000);
  peaks[1].Init(2000);
  workspace->GetPeaks(&peaks[0], NULL, &peaks[1], NULL, NULL);
  Compile(peaks, pb_peptide, compiler_prog1, compiler_prog2);
#endif

  Compile(workspace->GetPeaks(), pb_peptide, compiler_prog1, compiler_prog2);
#ifdef DEBUG
  if (Id() == FLAGS_debug_peptide_id) {
    cout << "Prog1:" << endl;
    DisAsm(prog1_);
    cout << "Prog2:" << endl;
    DisAsm(prog2_);
  }
#endif
}

// return the amino acid masses in the current peptide
vector<double> Peptide::getAAMasses() const {
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

// Probably defunct, uses old calling format.
int NoInlineDotProd(Peptide* peptide, const int* cache, int charge) {
  const void* prog = peptide->Prog(charge);
  int result;
#ifdef _MSC_VER
  // FIXME CEG add Windows compatible inline assembly
#else
  __asm__ __volatile__("call *%[prog]\n"
                       : "=a" (result)
                       : "d" (cache), [prog] "abcSD" (prog));
#endif
  return result;
}
