#include <iomanip>

#include "TideLiteMatchSet.h"
#include "util/Params.h"
#include "util/StringUtils.h"
#include "tide/peptide_lite.h"

SCORE_FUNCTION_T TideLiteMatchSet::curScoreFunction_ = INVALID_SCORE_FUNCTION;
ActivePeptideQueueLite* TideLiteMatchSet::active_peptide_queue_ = 0;
int TideLiteMatchSet::top_matches_ = 0;
int TideLiteMatchSet::decoy_num_ = 0;
int TideLiteMatchSet::mass_precision_ = 0;
int TideLiteMatchSet::score_precision_ = 0;
bool TideLiteMatchSet::concat_ = false;
string TideLiteMatchSet::decoy_prefix_ = "";

TideLiteMatchSet::TideLiteMatchSet() {
  psm_scores_processed_ = false;
};

TideLiteMatchSet::~TideLiteMatchSet() {
};

string TideLiteMatchSet::getHeader(TSV_OUTPUT_FORMATS_T format) { 

  string header;  // Return variable

  // columns are defined in ./src/io/MatchColumns.h and .cpp
  int* header_cols;
  int XCorr_cols[] = {
      FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
      PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
      BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, // TODO: add a variable for the matching theoretical peak series.
      DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
      PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
      DECOY_INDEX_COL
    };    
  int Pvalues_cols[] = {  //TODO: update the colums.
      FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
      PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
      DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
      PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
      DECOY_INDEX_COL
    };    
  int Diameter_cols[] = { //TODO: update the colums.
      FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
      PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
      BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, // TODO: add a variable for the matching theoretical peak series.
      DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
      PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
      DECOY_INDEX_COL
    };    
  size_t numHeaders;
  size_t i;  
  switch (format) {
    case TIDE_SEARCH_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          header_cols = XCorr_cols;
          numHeaders = sizeof(XCorr_cols) / sizeof(int);
          break;
        case PVALUES:
          header_cols = Pvalues_cols;
          numHeaders = sizeof(Pvalues_cols) / sizeof(int);
          break;
        case DIAMETER:
          header_cols = Diameter_cols;
          numHeaders = sizeof(Diameter_cols) / sizeof(int);
          break;
      }
      break;
    case MZTAB_TSV:
      break;
    case PIN_TSV:
      break;
  }
  header += get_column_header(header_cols[0]);
  for (i = 1; i < numHeaders; ++i) {
    header += '\t';
    header += get_column_header(header_cols[i]);
  }
  return header;
}

string TideLiteMatchSet::getReport(TSV_OUTPUT_FORMATS_T format, SpectrumCollection::SpecCharge* sc, bool target) { 
  PSMScores concat_or_target_psm_scores;
  PSMScores decoy_psm_scores;
  calculateAdditionalScoresAndGatherTargetsDecoys(concat_or_target_psm_scores, decoy_psm_scores);

  string spectrum_filename = "filename";

  PSMScores* scores = &concat_or_target_psm_scores;
  if (target == false) 
    scores = &decoy_psm_scores;

  string report;
  string separator("\t");
  for (PSMScores::iterator it = scores->begin(); it != scores->end(); ++it) {
    PeptideLite* peptide = active_peptide_queue_->GetPeptide((*it).ordinal_);
    string proteinNames;
    string flankingAAs;
    proteinNames = peptide->GetLocationStr(decoy_prefix_);
    flankingAAs = peptide->GetFlankingAAs();
    string peptide_with_mods = peptide->SeqWithMods(mass_precision_);
    string modifications = peptide->getModifications(); // todo: needs to be implelented
    string peptide_with_nomods = peptide->Seq();    

    switch (format) {
      case TIDE_SEARCH_TSV:
        switch (curScoreFunction_) {
          case XCORR_SCORE:
            report += spectrum_filename + separator;                                    // File name
            report += StringUtils::ToString(sc->spectrum->SpectrumNumber(), 0) + separator;                     //scan Id
            report += StringUtils::ToString(sc->charge, 0) + separator;    //Charge
            report += StringUtils::ToString(sc->spectrum->PrecursorMZ(), mass_precision_) + separator;  //spectrum precursor mz
            report += StringUtils::ToString((sc->spectrum->PrecursorMZ() - MASS_PROTON)*sc->charge, mass_precision_) + separator;  //spectrum neutral mass
            report += StringUtils::ToString(peptide->Mass(), mass_precision_) + separator; // calculated peptide mass.
            report += StringUtils::ToString((*it).delta_cn_, score_precision_) + separator;  // delta cn
            report += StringUtils::ToString((*it).delta_lcn_, score_precision_) + separator;  //delat lcn
            report += StringUtils::ToString((*it).xcorr_score_, score_precision_) + separator;  // xcorr score
            report += StringUtils::ToString((*it).tailor_, score_precision_) + separator;      // tailor score 
            report += StringUtils::ToString((*it).by_ion_matched_, score_precision_) + separator; //by ions matched
            report += StringUtils::ToString((*it).by_ion_total_, score_precision_) + separator;   // by ions total
            report += StringUtils::ToString((*it).by_ion_matched_/(*it).by_ion_total_, score_precision_) + separator;   // by by ions matched per total (fraction)
            report += StringUtils::ToString(target ? n_concat_or_target_matches_:n_decoy_matches_, 0) + separator;
            report += peptide_with_mods + separator;
            report += modifications + separator;  
            report += peptide_with_nomods + separator;
            report += proteinNames + separator;
            report += flankingAAs + separator;
            report += string(peptide->IsDecoy()? "decoy":"target") + separator;
            report += peptide->TargetSeq() + separator;
            report += StringUtils::ToString(peptide->DecoyIdx(), 0);
            break;
          case PVALUES:
            break;
          case DIAMETER:
            break;
        }
        break;
      case MZTAB_TSV:
        switch (curScoreFunction_) {
          case XCORR_SCORE:
            break;
          case PVALUES:
            break;
          case DIAMETER:
            break;
        }
        break;
      case PIN_TSV:
        switch (curScoreFunction_) {
          case XCORR_SCORE:
            break;
          case PVALUES:
            break;
          case DIAMETER:
            break;
        }
        break;
    }
  }      
}

void TideLiteMatchSet::calculateAdditionalScoresAndGatherTargetsDecoys(PSMScores& concat_or_target_psm_scores, PSMScores& decoy_psm_scores){ 
  if (psm_scores_processed_ == true)
    return;

  bool (*comp)(const Scores& x, const Scores& y);

  switch (curScoreFunction_) {
  case XCORR_SCORE:
    comp = &cmpXcorrScore;
    break;
  case PVALUES:
    comp = &cmpCombinedPvalue;
    break;
  case DIAMETER:
    break;
  }

  quantile_score_ = 1.0;

  // Calculate Tailor scores. Get the 99th quantile:
  if (curScoreFunction_ == XCORR_SCORE) {

    int quantile_pos = (int)(TAILOR_QUANTILE_TH*(double)psm_scores_.size()+0.5);
    if (quantile_pos < 3) 
      quantile_pos = 3;
    if (quantile_pos >= psm_scores_.size()) 
      quantile_pos = psm_scores_.size()-1;

    make_heap(psm_scores_.begin(), psm_scores_.end(), cmpXcorrScore);
    for (int i = 0; i < quantile_pos; ++i){
      pop_heap(psm_scores_.begin()+i, psm_scores_.end(), cmpXcorrScore);
    }
    quantile_score_ = psm_scores_[quantile_pos].xcorr_score_ +TAILOR_OFFSET; // Make sure scores positive
  }

  // Gather targets and decoys
  concat_or_target_psm_scores.clear();
  decoy_psm_scores.clear();

  map<int, int> decoyWriteCount;
  int gatherSize = top_matches_ + 1;  
  
  make_heap(psm_scores_.begin(), psm_scores_.end(), comp);

  for (PSMScores::iterator i = psm_scores_.end(); i != psm_scores_.begin(); ) {
    pop_heap(psm_scores_.begin(), i--, comp);

    PeptideLite* peptide = (active_peptide_queue_->GetPeptide((*i).ordinal_));
    if (concat_ || !peptide->IsDecoy()) {
      if (concat_or_target_psm_scores.size() < gatherSize) {
        concat_or_target_psm_scores.push_back(*i);
      }
    } else {
      int idx = peptide->DecoyIdx();
      map<int, int>::iterator j = decoyWriteCount.find(idx);
      if (j == decoyWriteCount.end()) {
        j = decoyWriteCount.insert(make_pair(idx, 0)).first;
      }
      if (j->second < gatherSize) {
        j->second++;
        decoy_psm_scores.push_back(*i);
      }
    }
    if ((concat_ && concat_or_target_psm_scores.size() >= gatherSize) || (!concat_ && concat_or_target_psm_scores.size() >= gatherSize && decoy_psm_scores.size() >= gatherSize*decoy_num_)){
      break;
    }
  }  

  //Calculate delta_cn and delta_lcn
  Scores last_psm = *std::min_element(psm_scores_.begin(), psm_scores_.end(), comp);
  double last_value = last_psm.xcorr_score_;

  sort(concat_or_target_psm_scores.begin(), concat_or_target_psm_scores.end(), comp);
  sort(decoy_psm_scores.begin(), decoy_psm_scores.end(), comp);

  for (PSMScores::iterator it = concat_or_target_psm_scores.begin(); it != concat_or_target_psm_scores.end(); ++it){
    (*it).tailor_ = (*it).xcorr_score_ / quantile_score_;
    (*it).delta_cn_ = (*it).xcorr_score_ / (*(it+1).xcorr_score_ ;
    (*it).delta_lcn_ = (*it).xcorr_score_ / last_value;
  }

  for (PSMScores::iterator it = decoy_psm_scores.begin(); it != decoy_psm_scores.end(); ++it){
    (*it).tailor_ = (*it).xcorr_score_ / quantile_score_;
    (*it).delta_cn_ = (*it).xcorr_score_ / (*(it+1).xcorr_score_ ;
    (*it).delta_lcn_ = (*it).xcorr_score_ / last_value;
  }

  psm_scores_processed_ = true;

}

