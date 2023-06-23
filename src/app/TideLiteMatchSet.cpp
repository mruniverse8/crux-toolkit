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
int TideLiteMatchSet::mod_precision_ = 0;
bool TideLiteMatchSet::concat_ = false;
string TideLiteMatchSet::decoy_prefix_ = "";


// columns are defined in ./src/io/MatchColumns.h and .cpp
int TideLiteMatchSet::XCorr_cols[] = {
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, BY_IONS_FRACTION_COL, // TODO: add a variable for the matching theoretical peak series.
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
 int TideLiteMatchSet::Pvalues_cols[] = {  //TODO: update the colums.
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, REFACTORED_SCORE_COL, EXACT_PVALUE_COL, 
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    
int TideLiteMatchSet::Diameter_cols[] = { //TODO: update the colums.
    FILE_COL, SCAN_COL, CHARGE_COL, SPECTRUM_PRECURSOR_MZ_COL, SPECTRUM_NEUTRAL_MASS_COL,
    PEPTIDE_MASS_COL, DELTA_CN_COL, DELTA_LCN_COL, XCORR_SCORE_COL, TAILOR_COL, 
    BY_IONS_MATCHED_COL, BY_IONS_TOTAL_COL, // TODO: add a variable for the matching theoretical peak series.
    DISTINCT_MATCHES_SPECTRUM_COL, SEQUENCE_COL, MODIFICATIONS_COL, UNMOD_SEQUENCE_COL,
    PROTEIN_ID_COL, FLANKING_AA_COL, TARGET_DECOY_COL, ORIGINAL_TARGET_SEQUENCE_COL,
    DECOY_INDEX_COL
  };    


TideLiteMatchSet::TideLiteMatchSet() {
  psm_scores_processed_ = false;
};

TideLiteMatchSet::~TideLiteMatchSet() {
};

int* TideLiteMatchSet::getColumns(TSV_OUTPUT_FORMATS_T format, size_t& numHeaders){
  switch (format) {
    case TIDE_SEARCH_TSV:
      switch (curScoreFunction_) {
        case XCORR_SCORE:
          numHeaders = sizeof(XCorr_cols) / sizeof(int);
          return XCorr_cols;
        case PVALUES:
          numHeaders = sizeof(Pvalues_cols) / sizeof(int);
          return Pvalues_cols;
        case DIAMETER:
          numHeaders = sizeof(Diameter_cols) / sizeof(int);
          return Diameter_cols;
      }
      break;
    case MZTAB_TSV:
      break;
    case PIN_TSV:
      break;
  }
}

string TideLiteMatchSet::getHeader(TSV_OUTPUT_FORMATS_T format) { 

  string header;  // Return variable

  size_t numHeaders;
  int* header_cols = getColumns(format, numHeaders);

  header += get_column_header(header_cols[0]);
  for (size_t i = 1; i < numHeaders; ++i) {
    header += '\t';
    header += get_column_header(header_cols[i]);
  }
  return header;
}

void TideLiteMatchSet::getReport(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, SpectrumCollection::SpecCharge* sc, string &concat_or_target_report, string& decoy_report) { 
  
  PSMScores concat_or_target_psm_scores;
  PSMScores decoy_psm_scores;

  // Get the top_n matches target and decoy PSMs
  gatherTargetsDecoys(concat_or_target_psm_scores, decoy_psm_scores); // output variables. decoy_psm_scores is empty in case of concat=True or PIN or MZTAB formats
  
  //calculate tailor, delta_cn, and delta_lcn for the top n matches
  calculateAdditionalScores(concat_or_target_psm_scores);  
  calculateAdditionalScores(decoy_psm_scores);  // decoy_psm_scores is empty in case of concat or in case of PIN or MZTAB formats

  return;
  // Prepare the results in a string
  printResults(format, spectrum_filename, sc, concat_or_target_psm_scores, concat_or_target_report);
  printResults(format, spectrum_filename, sc, decoy_psm_scores, decoy_report); // decoy_report is empty string if decoy_psm_scores is empty 

}

void TideLiteMatchSet::gatherTargetsDecoys(PSMScores& concat_or_target_psm_scores, PSMScores& decoy_psm_scores){ 
  if (psm_scores_processed_ == true)
    return;

  // reset output variables
  concat_or_target_psm_scores.clear();
  decoy_psm_scores.clear();

  if (psm_scores_.size() == 0) {
    psm_scores_processed_ = true;
    return;
  }
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

    int quantile_pos = (int)(TAILOR_QUANTILE_TH*(double)psm_scores_.size()+0.5)-1; // zero indexed
    printf("quantile:%d\n", quantile_pos);
    if (quantile_pos < 2) 
      quantile_pos = 2;  // the third element
    printf("quantile:%d\n", quantile_pos);
    if (quantile_pos >= psm_scores_.size()) 
      quantile_pos = psm_scores_.size()-1; // the last element

    printf("quantile:%d\n", quantile_pos);
    printf("Size:%d\n", psm_scores_.size());
    make_heap(psm_scores_.begin(), psm_scores_.end(), cmpXcorrScore);
    for (int i = 0; i <= quantile_pos; ++i){
      pop_heap(psm_scores_.begin(), psm_scores_.end()-i, cmpXcorrScore);
      Scores back = psm_scores_[psm_scores_.size()-1-i];
      printf("%lf\n", back.xcorr_score_);
    }
    quantile_score_ = psm_scores_[psm_scores_.size()-1-quantile_pos].xcorr_score_ +TAILOR_OFFSET; // Make sure scores positive
    printf("quantile score:%lf\n", quantile_score_);
  }

  // get the value of the last score for the delta_lcn scores
  last_psm_ = std::min_element(psm_scores_.begin(), psm_scores_.end(), comp);
  printf("last value:%lf\n", (*last_psm_).xcorr_score_);
  
  // Gather target and decoy PSMs
  int gatherSize = top_matches_ + 1; // Get one more psms than the top-matches, so the delta_cn can be calculated correctly for the last element.
  map<int, int> decoyWriteCount;
  
  make_heap(psm_scores_.begin(), psm_scores_.end(), comp);

  for (PSMScores::iterator i = psm_scores_.end(); i != psm_scores_.begin(); ) {
    pop_heap(psm_scores_.begin(), i--, comp);
    printf("pop: %lf\n", (*i).xcorr_score_);

// TODO: uncommented for debugging, put it back
    PeptideLite* peptide = (active_peptide_queue_->GetPeptide((*i).ordinal_));
    if (concat_ || !peptide->IsDecoy()) {
//    if (concat_ ) {
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
}

void TideLiteMatchSet::calculateAdditionalScores(PSMScores& psm_scores) {  // Additional scores are:  delta_cn, delta_lcn, tailor;
  // The  gatherTargetsDecoys must be run before calling this function.

  switch (curScoreFunction_) {
  case XCORR_SCORE:
  case DIAMETER:
    for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it){
      (*it).tailor_ = ((*it).xcorr_score_  + TAILOR_OFFSET )/ quantile_score_;
      (*it).delta_lcn_ = ((*it).xcorr_score_ - (*last_psm_).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = ((*it).xcorr_score_ - (*(it+1)).xcorr_score_)/max((*it).xcorr_score_, 1.0);
      else 
        (*it).delta_cn_ = 0.0;
      printf("xcorr:%lf, tailor:%lf, delta_cn: %lf, delta_lcn: %lf\n", (*it).xcorr_score_, (*it).tailor_, (*it).delta_cn_, (*it).delta_lcn_);
    }
    break;
  case PVALUES:
  case PVALUES_HR:
  case PVALUES_LR:
    for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it){
      (*it).delta_lcn_ = -log10((*it).combinedPval_) + log10((*last_psm_).combinedPval_);
      if (it != psm_scores.end()-1)
        (*it).delta_cn_ = -log10((*it).combinedPval_) + log10((*(it+1)).combinedPval_);
      else 
        (*it).delta_cn_ = 0.0;
      printf("xcorr:%lf, tailor:%lf, delta_cn: %lf, delta_lcn: %lf\n", (*it).xcorr_score_, (*it).tailor_, (*it).delta_cn_, (*it).delta_lcn_);
    }
    break;
  case HYPERSCORE:
  case HYPERSCORE_LA:

    break;
  }

  psm_scores_processed_ = true;

}

void TideLiteMatchSet::printResults(TSV_OUTPUT_FORMATS_T format, string spectrum_filename, SpectrumCollection::SpecCharge* sc, PSMScores& psm_scores, string& report) { 
  // The order of the fields of the results is solely based on the column order
  size_t numHeaders;
  int* header_cols = getColumns(format, numHeaders);
  int cnt = 0;
  for (PSMScores::iterator it = psm_scores.begin(); it != psm_scores.end(); ++it, ++cnt) {
    if (cnt > top_matches_)
      break;
    PeptideLite* peptide = active_peptide_queue_->GetPeptide((*it).ordinal_);

    string proteinNames = peptide->GetLocationStr(decoy_prefix_);
    string flankingAAs = peptide->GetFlankingAAs();
    string peptide_with_mods = peptide->SeqWithMods(mod_precision_);
    string modifications = peptide->getModifications(mod_precision_); // todo: needs to be tested
    string peptide_without_mods = peptide->Seq();   
    // The order of the fields depends on the columns defined at the beginning of this file.
    // The fields below can be in arbitrary order.
    for (size_t i = 0; i < numHeaders; ++i) {   
      switch (header_cols[i]){
      case FILE_COL:
        report += spectrum_filename;
        break;
      case SCAN_COL:
        report += StringUtils::ToString(sc->spectrum->SpectrumNumber(), 0); // Scan Id
        break;
      case CHARGE_COL:
        report += StringUtils::ToString(sc->charge, 0);                     // Charge
        break;
      case SPECTRUM_PRECURSOR_MZ_COL:
        report += StringUtils::ToString(sc->spectrum->PrecursorMZ(), mass_precision_);                             //spectrum precursor mz
        break;
      case SPECTRUM_NEUTRAL_MASS_COL:
        report += StringUtils::ToString((sc->spectrum->PrecursorMZ() - MASS_PROTON)*sc->charge, mass_precision_);  // spectrum neutral mass
        break;
      case PEPTIDE_MASS_COL:
        report += StringUtils::ToString(peptide->Mass(), mass_precision_);                                         // spectrum neutral mass
        break;
      case DELTA_CN_COL:
        report += StringUtils::ToString((*it).delta_cn_, score_precision_);         // delta_cn
        break;
      case DELTA_LCN_COL:
        report += StringUtils::ToString((*it).delta_lcn_, score_precision_);        // delta_lcn
        break;
      case XCORR_SCORE_COL:
        report += StringUtils::ToString((*it).xcorr_score_, score_precision_);      // xcorr score
        break;
      case TAILOR_COL:
        report += StringUtils::ToString((*it).tailor_, score_precision_);           // tailor score
        break;
      case BY_IONS_MATCHED_COL:
        report += StringUtils::ToString((*it).by_ion_matched_, score_precision_);   // by ions matched
        break;
      case BY_IONS_TOTAL_COL:
        report += StringUtils::ToString((*it).by_ion_total_, score_precision_);     // total no of by ions
        break;
      case BY_IONS_FRACTION_COL:
        report += StringUtils::ToString((*it).by_ion_matched_/(*it).by_ion_total_, score_precision_);  // fraction of the matched per total by-ions
        break;
      case DISTINCT_MATCHES_SPECTRUM_COL:
        report += "distinct matches";//StringUtils::ToString(target ? n_concat_or_target_matches_:n_decoy_matches_, 0);  // no candidate peptides
        break;
      case SEQUENCE_COL:
        report += peptide_with_mods;        // peptide sequence with modifications
        break;
      case MODIFICATIONS_COL:
        report += modifications;            // the list of modifications in the peptide
        break;
      case UNMOD_SEQUENCE_COL:
        report += peptide_without_mods;     // plain peptide sequence stripped with mods
        break;
      case PROTEIN_ID_COL:
        report += proteinNames;             // protein IDs
        break;
      case FLANKING_AA_COL:
        report += flankingAAs;              // flanking Amino Acids 
        break;
      case TARGET_DECOY_COL:
        report += string(peptide->IsDecoy()? "decoy":"target");  // target or decoy
        break;
      case ORIGINAL_TARGET_SEQUENCE_COL:
        report += peptide->TargetSeq();     // original target sequences
        break;
      case DECOY_INDEX_COL:
        report += StringUtils::ToString(peptide->DecoyIdx(), 0);  // decoy index
        break;
      }
      if (i < numHeaders-1)  // If not the last column, add a column separator
        report += '\t';
    }
  }  
}