#ifndef TIDE_LITE_MATCH_SET_H
#define TIDE_LITE_MATCH_SET_H

#define  BOOST_DATE_TIME_NO_LIB
#include <boost/thread.hpp>
#include <vector>
#include "raw_proteins.pb.h"
#include "tide/records.h"
#include "tide/active_peptide_queue.h"  // no include guard
#include "tide/fixed_cap_array.h"
#include "tide/peptide.h"
#include "tide/sp_scorer.h"
#include "tide/spectrum_collection.h"
#include "tide/ActivePeptideQueueLite.h"

#include "model/Modification.h"
#include "model/PostProcessProtein.h"

using namespace std;

class TideLiteMatchSet {    
 public:
  class Scores {
   public:
    int ordinal_;
    double xcorr_score_;
    double exact_pval_;
    double refactored_xcorr_;
    int resEv_score_;
    double resEv_pval_;
    double combinedPval_;
    double tailor_; 
    int by_ion_matched_; 
    int by_ion_total_;     
    double sp_score_;
    double hyper_score_;
    double hyper_score_la_; //maybe to be removed.
    double delta_cn_;
    double delta_lcn_;
    Scores():ordinal_(0), xcorr_score_(0), exact_pval_(0), refactored_xcorr_(0), 
      resEv_score_(0), resEv_pval_(0), combinedPval_(0), tailor_(0), by_ion_matched_(0), by_ion_total_(0), 
      sp_score_(0), hyper_score_(0), hyper_score_la_(0), delta_cn_(0), delta_lcn_(0){}
  };
  typedef FixedCapacityArray<Scores> PSMScores;
  PSMScores psm_scores_;   // This one is used to gather psms during scoring.

  int n_concat_or_target_matches_;  // concat or target
  int n_decoy_matches_;
  bool psm_scores_processed_;

  static bool cmpXcorrScore(const Scores& x, const Scores& y) {  // compare PSMs by xcorr scores. Larger comes first
    return x.xcorr_score_ < y.xcorr_score_;
  }  
  static bool cmpCombinedPvalue(const Scores& x, const Scores& y) {  // compare PSMs by P-values. smaller comes first
    return x.xcorr_score_ > y.xcorr_score_;
  }  
  static bool cmpHyperScore(const Scores& x, const Scores& y) {  // compare PSMs by P-values. larger comes first
    return x.hyper_score_ < y.hyper_score_;
  }  

  TideLiteMatchSet();
  ~TideLiteMatchSet();

  string getHeader(TSV_OUTPUT_FORMATS_T format); // pass filetype;
  string getReport(TSV_OUTPUT_FORMATS_T format, 
                   SpectrumCollection::SpecCharge* sc,
                   bool target
                   );
  void calculateAdditionalScoresAndGatherTargetsDecoys(PSMScores& concat_or_target_psm_scores, PSMScores& decoy_psm_scores);  // Additional scores are:  delta_cn, delta_lcn, tailor

  /* Constants required for the tailor scoring */
  const double TAILOR_QUANTILE_TH = 0.01;
  const double TAILOR_OFFSET = 5.0 ;
  double quantile_score_;

  // Global static parameters
  static SCORE_FUNCTION_T curScoreFunction_;
  static ActivePeptideQueueLite* active_peptide_queue_;  
  static int top_matches_;
  static int decoy_num_;
  static int mass_precision_;
  static int score_precision_;
  static bool concat_;
  static string decoy_prefix_;



};

#endif
