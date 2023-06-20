#include <deque>
#include "peptides.pb.h"
#include "peptide_lite.h"
#include "theoretical_peak_set.h"
#include "fifo_alloc.h"
#include "spectrum_collection.h"
#include "io/OutputFiles.h"

#ifndef ACTIVE_LITE_PEPTIDE_QUEUE_H
#define ACTIVE_LITE_PEPTIDE_QUEUE_H

class TheoreticalPeakCompiler;

class ActivePeptideQueueLite {
 public:
  ActivePeptideQueueLite(RecordReader* reader,
        const vector<const pb::Protein*>& proteins,
        vector<const pb::AuxLocation*>* locations=NULL, 
        bool dia_mode = false);

  ~ActivePeptideQueueLite();

  int SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, 
        double min_range, double max_range); 

  PeptideLite* GetPeptide(int index) {
    return *(begin_ + index);
  }

  int nCandPeptides_;
  deque<PeptideLite*> queue_;
  deque<PeptideLite*>::const_iterator begin_, end_;  
  vector<bool> candidatePeptideStatus_;
  int min_candidates_;
  bool dia_mode_;

 private:
  bool isWithinIsotope(vector<double>* min_mass,
        vector<double>* max_mass, 
        double mass,
        int* isotope_idx);   

  void ComputeTheoreticalPeaksBack();    

  RecordReader* reader_;
  const vector<const pb::Protein*>& proteins_; 
  vector<const pb::AuxLocation*>* locations_;

  
  TheoreticalPeakSetBYSparse theoretical_peak_set_;
  pb::Peptide current_pb_peptide_;
};

#endif
