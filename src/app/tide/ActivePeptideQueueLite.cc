// original author: Benjamin Diament
// subsequently modified by Attila Kertesz-Farkas, Jeff Howbert
#include <deque>
#include <gflags/gflags.h>
#include "records.h"
#include "peptides.pb.h"
#include "peptide_lite.h"
#include "ActivePeptideQueueLite.h"
#include "records_to_vector-inl.h"
#include "theoretical_peak_set.h"
#include "compiler.h"
#include "app/TideLiteMatchSet.h"
#include <map> 
#define CHECK(x) GOOGLE_CHECK((x))

ActivePeptideQueueLite::ActivePeptideQueueLite(RecordReader* reader,
                                       const vector<const pb::Protein*>& proteins, 
                                       vector<const pb::AuxLocation*>* locations, 
                                       bool dia_mode)
  : reader_(reader),
    proteins_(proteins),
    theoretical_peak_set_(1000),   // probably overkill, but no harm
    locations_(locations),
    dia_mode_(dia_mode) {
  CHECK(reader_->OK());
  min_candidates_ = 30;
}

ActivePeptideQueueLite::~ActivePeptideQueueLite() {
}

// Compute the theoretical peaks of the peptide in the "back" of the queue
// (i.e. the one most recently read from disk -- the heaviest).
void ActivePeptideQueueLite::ComputeTheoreticalPeaksBack() {
  theoretical_peak_set_.Clear();
  PeptideLite* peptide = queue_.back();
  peptide->ComputeTheoreticalPeaks(&theoretical_peak_set_);
}

bool ActivePeptideQueueLite::isWithinIsotope(vector<double>* min_mass, vector<double>* max_mass, double mass, int* isotope_idx) {
  for (int i = *isotope_idx; i < min_mass->size(); ++i) {
    if (mass >= (*min_mass)[i] && mass <= (*max_mass)[i]) {
      if (i > *isotope_idx) {
        *isotope_idx = i;
      }
      return true;
    }
  }
  return false;
}

int ActivePeptideQueueLite::SetActiveRange(vector<double>* min_mass, vector<double>* max_mass, double min_range, double max_range) {
  //min_range and max_range have been introduced to fix a bug
  //introduced by m/z selection. see #222 in sourceforge
  //this has to be true:
  // min_range <= min_mass <= max_mass <= max_range

  // queue front() is lightest; back() is heaviest

  // delete anything already loaded that falls below min_range
  while (!queue_.empty() && queue_.front()->Mass() < min_range) {
    PeptideLite* peptide = queue_.front();
    queue_.pop_front();
    delete peptide;
  }

  // Enqueue all peptides that are not yet queued but are lighter than
  // max_range. For each new enqueued peptide compute the corresponding
  // theoretical peaks. Data associated with each peptide is allocated by
  // fifo_alloc_peptides_.
  bool done = false;
  //Modified for tailor score calibration method by AKF
  if (queue_.empty() || queue_.back()->Mass() <= max_range || queue_.size() < min_candidates_) {
    if (!queue_.empty()) {
      ComputeTheoreticalPeaksBack();
    }
    while (!(done = reader_->Done())) {
      // read all peptides lighter than max_range
      reader_->Read(&current_pb_peptide_);
      if (current_pb_peptide_.mass() < min_range) {
        // we would delete current_pb_peptide_;
        continue; // skip peptides that fall below min_range
      }
      PeptideLite* peptide = new PeptideLite(current_pb_peptide_, proteins_, locations_);
      assert(peptide != NULL);
      queue_.push_back(peptide);
      //Modified for tailor score calibration method by AKF
      if (peptide->Mass() > max_range && queue_.size() > min_candidates_) {
        break;
      }
      ComputeTheoreticalPeaksBack();
    }
  }
  // by now, if not EOF, then the last (and only the last) enqueued
  // peptide is too heavy
  assert(!queue_.empty() || done);

  // Set up iterator for use with HasNext(),
  // GetPeptide(), and NextPeptide(). Return the number of enqueued peptides.
  if (queue_.empty()) {
    return 0;
  }

  begin_ = queue_.begin();
  candidatePeptideStatus_.empty();
  while (begin_ != queue_.end() && (*begin_)->Mass() < min_mass->front()) {
    ++begin_;
    candidatePeptideStatus_.push_back(false);  
  }
  end_ = begin_;
  begin_ = queue_.begin();
  int* isotope_idx = new int(0);
  nCandPeptides_ = 0;
  while (end_ != queue_.end() && (*end_)->Mass() < max_mass->back() ) {
    if (isWithinIsotope(min_mass, max_mass, (*end_)->Mass(), isotope_idx)) {
      ++nCandPeptides_;
      candidatePeptideStatus_.push_back(true);
    } else {
      candidatePeptideStatus_.push_back(false);
    }
    ++end_;
  }
  delete isotope_idx;
  if (nCandPeptides_ == 0) {
    return 0;
  }
  for (; end_ != queue_.end(); ++end_) {
    if ((*end_)->peaks_0.size() == 0 || candidatePeptideStatus_.size() >= min_candidates_-1) {
      break;
    }
    candidatePeptideStatus_.push_back(false);
  }
  return nCandPeptides_;
}
