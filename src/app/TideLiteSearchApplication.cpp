#include <cstddef>
#include <cstdio>
#include "app/tide/abspath.h"
#include "app/tide/records_to_vector-inl.h"

#include "io/carp.h"
#include "parameter.h"
#include "io/SpectrumRecordWriter.h"
//#include "TideIndexApplication.h"
#include "TideLiteSearchApplication.h"
#include "ParamMedicApplication.h"
#include "PSMConvertApplication.h"
#include "tide/mass_constants.h"
#include "TideLiteMatchSet.h"
#include "tide/spectrum_collection.h"
#include "util/Params.h"
#include "util/FileUtils.h"
#include "util/StringUtils.h"
#include <math.h> 
#include <map>
#include "tide/ActivePeptideQueueLite.h"

//here :O
#include <iostream>
#include <set>

/* This constant is the product of the original "magic number" (10000,
 * on line 4622 of search28.c) that was used to rescale the XCorr
 * score, and the integerization constant used by Benjamin Diament in
 * Tide.  In the Tide publication, that constant was reported as 10^7.
 * --WSN, 10 March 2015 */
const double TideLiteSearchApplication::XCORR_SCALING = 100000000.0;

/* This constant is used to put the refactored XCorr back into the
 * same range as the original XCorr score.  It is the XCorr "magic
 * number" (10000) divided by the EVIDENCE_SCALE_INT (defined in
 * tide/spectrum_preprocess2.cc). */
const double TideLiteSearchApplication::RESCALE_FACTOR = 20.0;



TideLiteSearchApplication::TideLiteSearchApplication() {
  remove_index_ = "";
  spectrum_flag_ = NULL;
  decoy_num_ = 0;
  num_range_skipped_ = 0;
  num_precursors_skipped_ = 0;
  num_isotopes_skipped_ = 0;
  num_retained_ = 0;
  total_candidate_peptides_ = 0;
  precursor_window_ = 0;
  spectrum_min_mz_ = 0; 
  spectrum_max_mz_ = 0;
  min_scan_ = 0;
  max_scan_ = 0;
  min_peaks_ = 0;
  min_precursor_charge_ = 0;
  max_precursor_charge_ = 0;
  tsv_output_ = false;
  mztab_output_ = false;
  pin_output_ = false;
  target_file_ = NULL;
  decoy_file_ = NULL;

  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    locks_array_.push_back(new boost::mutex());
  }
}

TideLiteSearchApplication::~TideLiteSearchApplication() {
  for (int i = 0; i < NUMBER_LOCK_TYPES; i++) {
    delete locks_array_[i];
  }
}

int TideLiteSearchApplication::main(int argc, char** argv) {
  return main(Params::GetStrings("tide spectra file"));
}

int TideLiteSearchApplication::main(const vector<string>& input_files) {
  return main(input_files, Params::GetString("tide database"));
}

int TideLiteSearchApplication::main(const vector<string>& input_files, const string input_index) {
  carp(CARP_INFO, "Running tide-lite-search...");

  print_interval_ = Params::GetInt("print-search-progress");  
  
  bin_width_  = Params::GetDouble("mz-bin-width");
  bin_offset_ = Params::GetDouble("mz-bin-offset");

  use_neutral_loss_peaks_ = Params::GetBool("use-neutral-loss-peaks");
  use_flanking_peaks_ = Params::GetBool("use-flanking-peaks");

  negative_isotope_errors_ = getNegativeIsotopeErrors();  

  window_type_= string_to_window_type(Params::GetString("precursor-window-type"));
  precursor_window_ = Params::GetDouble("precursor-window");

  spectrum_min_mz_ =  Params::GetDouble("spectrum-min-mz") ;
  spectrum_max_mz_ = Params::GetDouble("spectrum-max-mz") ;
  min_peaks_ = Params::GetInt("min-peaks");
  min_precursor_charge_ = Params::GetInt("min-precursor-charge");
  max_precursor_charge_ = Params::GetInt("max-precursor-charge");  

  top_matches_ = Params::GetInt("top-match");  

  // Check scan-number parameter
  string scan_range = Params::GetString("scan-number");
  if (scan_range.empty()) {
    min_scan_ = 0;
    max_scan_ = BILLION;
    carp(CARP_DEBUG, "Searching all scans");
  } else if (scan_range.find('-') == string::npos) {
    // Single scan
    min_scan_ = max_scan_ = atoi(scan_range.c_str());
    carp(CARP_INFO, "Searching single scan %d", min_scan_);
  } else {
    if (!get_range_from_string(scan_range.c_str(), min_scan_, max_scan_)) {
      carp(CARP_FATAL, "The scan number range '%s' is invalid. "
           "Must be of the form <first>-<last>.", scan_range.c_str());
    }
    if (min_scan_ > max_scan_) {
      carp(CARP_FATAL, "Invalid scan range: %d to %d.", min_scan_, max_scan_);
    }
    carp(CARP_INFO, "Searching scan range %d to %d.", min_scan_, max_scan_);
  }

  // Determine which score function to use for scoring PSMs and store in SCORE_FUNCTION enum. TODO: update this in ./src/util/crux-utils.cpp: 68
  curScoreFunction_ = string_to_score_function_type(Params::GetString("score-function"));  
 
  // Get a peptide reader to the peptide index datasets along with proteins, auxlocs. 
  vector<const pb::Protein*> proteins;
  vector<const pb::AuxLocation*> locations;
  pb::Header peptides_header;
  string peptides_file = FileUtils::Join(input_index, "pepix");  
  HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
  getPeptideIndexData(input_index, proteins, locations, peptides_header);
  // Create active peptide queue
  active_peptide_queue_ = new ActivePeptideQueueLite(peptide_reader.Reader(), proteins, &locations);  

  // Create the output files
  createOutputFiles(&target_file_, &decoy_file_);

  TideLiteMatchSet::curScoreFunction_ = curScoreFunction_;
  TideLiteMatchSet::active_peptide_queue_ = active_peptide_queue_;
  TideLiteMatchSet::top_matches_ = top_matches_;
  TideLiteMatchSet::decoy_num_ = decoy_num_;
  TideLiteMatchSet::mass_precision_ =  Params::GetInt("mass-precision");
  TideLiteMatchSet::score_precision_ = Params::GetInt("precision");
  TideLiteMatchSet::mod_precision_ = Params::GetInt("mod-precision");
  TideLiteMatchSet::concat_ = Params::GetBool("concat");
  TideLiteMatchSet::decoy_prefix_ = top_matches_;

  TSV_OUTPUT_FORMATS_T format = TIDE_SEARCH_TSV;
  string header = TideLiteMatchSet::getHeader(format);
  printf("%s\n", header.c_str());

  TideLiteMatchSet psm_scores;
  for (int i = 0; i < 1500; ++i){
    TideLiteMatchSet::Scores currScore;
    currScore.xcorr_score_ = i+15.0;
    currScore.ordinal_ = i;
    psm_scores.psm_scores_.push_back(currScore);
  }
  SpectrumCollection::SpecCharge* sc = NULL; // get this data from the queue
  string concat_or_target_report;
  string decoy_report;
//  psm_scores.getReport(format, string("spectrum_filename"), sc, concat_or_target_report, decoy_report);  // get targets
  printf("%s\n", concat_or_target_report.c_str());




  // Convert the input spectrum data files to spectrumRecords if needed
/*  vector<InputFile> sr = getInputFiles(input_files);
  // Loop through spectrum files
  // Here :O
  vector<vector<Spectrum*>*> vect_spectra;
  vector<int> idxs_spectra(sr.size(),0);
  struct pair_mass_spectra{
    pair_mass_spectra(vector<Spectrum*> *ptr2, int idxx): ptr(ptr2), idx(idxx) {
    }
    vector<Spectrum*>* ptr;
    int idx;
    bool operator<(const pair_mass_spectra& other) const{
      return ptr->at(idx)->PrecursorMZ() < other.ptr->at(idx)->PrecursorMZ();
    }
  };
  set<pair_mass_spectra> queue_spectra;
  for (vector<InputFile>::const_iterator f = sr.begin(); f != sr.end(); f++) {
    string spectra_file = f->SpectrumRecords;
    carp(CARP_INFO, "Reading spectrum file %s.", spectra_file.c_str());
    // here changes
    SpectrumCollection* spectra = NULL;
    spectra = loadSpectra(spectra_file); // heap memory!

    vect_spectra.push_back(spectra->Spectra());
    if (!vect_spectra.back()->empty()) {
      queue_spectra.insert(pair_mass_spectra(vect_spectra.back(), 0));
    }
  }
  std::cout<<"Printing in order!!!\n";
  while(!queue_spectra.empty()) {
    auto it = queue_spectra.begin();
    auto idx = it->idx;
    auto ptr = it->ptr;
    std::cout<<ptr->at(idx)->PrecursorMZ()<<" ";
    queue_spectra.erase(it);
    idx++;
    if (ptr->size() == idx) {
      continue;
    }
    queue_spectra.insert(pair_mass_spectra(ptr, idx));
  }
  //merge done!!!
    // Delete temporary spectrumrecords file
  for (vector<InputFile>::const_iterator f = sr.begin(); f != sr.end(); f++) {
    if (!f->Keep) {
      string spectra_file = f->SpectrumRecords;      
      carp(CARP_DEBUG, "Deleting %s", spectra_file.c_str());
      remove(spectra_file.c_str());
    }
  }
*/
  delete active_peptide_queue_;
  return 0;
}


void TideLiteSearchApplication::search() {
  
  // This is the main search loop.

                           
  SpectrumCollection::SpecCharge* sc; // get this data from the queue.

  Spectrum* spectrum = sc->spectrum;
  double precursor_mz = spectrum->PrecursorMZ();
  int charge = sc->charge;
  int scan_num = spectrum->SpectrumNumber();

  if (precursor_mz < spectrum_min_mz_|| 
      precursor_mz > spectrum_max_mz_ || 
      scan_num < min_scan_ || 
      scan_num > max_scan_ ||
      spectrum->Size() < min_peaks_  ||
      charge < min_precursor_charge_ || 
      charge >max_precursor_charge_ ) {
    //continue; ///switch to continue in a loop
    return;
  }

  if (spectrum_flag_ != NULL) {  // TODO: Do something, possibly for cascade search
  }

  vector<double>* min_mass = new vector<double>();
  vector<double>* max_mass = new vector<double>();
  double min_range, max_range;
  
  computeWindow(*sc, min_mass, max_mass, &min_range, &max_range);
  active_peptide_queue_->SetActiveRange(min_mass, max_mass, min_range, max_range);

  if (active_peptide_queue_->nCandPeptides_ == 0) { // No peptides to score.
    return; ///switch to continue in a loop
  }
  locks_array_[LOCK_CANDIDATES]->lock();
  total_candidate_peptides_ += active_peptide_queue_->nCandPeptides_;
  locks_array_[LOCK_CANDIDATES]->unlock();  

  TideLiteMatchSet psm_scores;
  switch (curScoreFunction_) {
    case XCORR_SCORE:
      XCorrScoring(sc, psm_scores);
      break;
    case PVALUES:
      PValueScoring();
      break;
  }
  TSV_OUTPUT_FORMATS_T format = TIDE_SEARCH_TSV;
  string concat_or_target_report;
  string decoy_report;
//  psm_scores.getReport(format, string("spectrum_filename"), sc, concat_or_target_report, decoy_report);  // get targets
//  string report = psm_scores.getReport(format, sc, false);  // get decoys
  locks_array_[LOCK_RESULTS]->lock();
  //write restuls
  // target_file << results
  locks_array_[LOCK_RESULTS]->unlock();
}

void TideLiteSearchApplication::XCorrScoring(SpectrumCollection::SpecCharge* sc, TideLiteMatchSet& psm_scores){
  // scpectrum preprocessing for xcorr scoring

  ObservedPeakSet observed(bin_width_, bin_offset_,
                           use_neutral_loss_peaks_,
                           use_flanking_peaks_);

  long num_range_skipped = 0;
  long num_precursors_skipped = 0;
  long num_isotopes_skipped = 0;
  long num_retained = 0;

  int charge = sc->charge;  

  observed.PreprocessSpectrum(*(sc->spectrum), charge, &num_range_skipped,
                              &num_precursors_skipped,
                              &num_isotopes_skipped, &num_retained);

  // locks_array_[LOCK_REPORTING]->lock();
  num_range_skipped_ += num_range_skipped;
  num_precursors_skipped_ += num_precursors_skipped;
  num_isotopes_skipped_ += num_isotopes_skipped;
  num_retained_ += num_retained;
  // locks_array_[LOCK_REPORTING]->unlock();

  // Score the inactive peptides in the peptide queue if the number of nCadPeptides 
  // is less than the minimum for Tailor scoring

  bool score_inactive_peptides = true;
  if (active_peptide_queue_->min_candidates_ < active_peptide_queue_->nCandPeptides_)
    score_inactive_peptides = false;
  
  //Actual Xcorr Scoring 
  const int* cache = observed.GetCache();        
  int cnt = 0;
  for (deque<PeptideLite*>::const_iterator iter = active_peptide_queue_->begin_; 
    iter != active_peptide_queue_->end_;
    ++iter, ++cnt) {
    if (active_peptide_queue_->candidatePeptideStatus_[cnt] == false || score_inactive_peptides) 
      continue;
    TideLiteMatchSet::Scores currScore;
    int xcorr = 0;
    int match_cnt = 0;

    // Score with single charged theoretical peaks
    for (vector<unsigned int>::const_iterator iter_uint = (*iter)->peaks_0.begin();
      iter_uint != (*iter)->peaks_0.end();
      iter_uint++) {
      xcorr += cache[*iter_uint];
      if (cache[*iter_uint] > 0) {
        ++match_cnt;
      }
    }
    // Score with double charged theoretical peaks
    if (charge > 2){
      for (vector<unsigned int>::const_iterator iter_uint = (*iter)->peaks_1.begin();
        iter_uint != (*iter)->peaks_1.end();
        iter_uint++) {
        xcorr += cache[*iter_uint];
        if (cache[*iter_uint] > 0) {
          ++match_cnt;
        }
      }
    }
    currScore.ordinal_ = cnt;
    currScore.xcorr_score_ = (double)xcorr/XCORR_SCALING;// (double)(it->first / XCORR_SCALING)
    currScore.by_ion_matched_ = match_cnt;
    currScore.by_ion_total_ = (*iter)->peaks_0.size();
    if (charge > 2){
      currScore.by_ion_total_ += (*iter)->peaks_1.size();
    }
    psm_scores.psm_scores_.push_back(currScore);
  } 
}

void TideLiteSearchApplication::PValueScoring(){
//  double precursorMass = sc->neutral_mass;  //Added by Andy Lin (needed for residue evidence); TODO: rename neutral_mass

}

void TideLiteSearchApplication::computeWindow(
  const SpectrumCollection::SpecCharge& sc,
  vector<double>* out_min,
  vector<double>* out_max,
  double* min_range,
  double* max_range
) {
  double unit_dalton = BIN_WIDTH;
  switch (window_type_) {
  case WINDOW_MASS:
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back(sc.neutral_mass + (*ie * unit_dalton) - precursor_window_);
      out_max->push_back(sc.neutral_mass + (*ie * unit_dalton) + precursor_window_);
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors_.front() * unit_dalton)) - precursor_window_;
    *max_range = (sc.neutral_mass + (negative_isotope_errors_.back() * unit_dalton)) + precursor_window_;
    break;
  case WINDOW_MZ: {
    double mz_minus_proton = sc.spectrum->PrecursorMZ() - MASS_PROTON;
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back((mz_minus_proton - precursor_window_) * sc.charge + (*ie * unit_dalton));
      out_max->push_back((mz_minus_proton + precursor_window_) * sc.charge + (*ie * unit_dalton));
    }
    // N.B. Next two lines used to have max_charge instead of sc.charge. Not sure why. --WSN 9 Jun 2021
    *min_range = (mz_minus_proton*sc.charge + (negative_isotope_errors_.front() * unit_dalton)) - precursor_window_*sc.charge;
    *max_range = (mz_minus_proton*sc.charge + (negative_isotope_errors_.back() * unit_dalton)) + precursor_window_*sc.charge;
    break;
  }
  case WINDOW_PPM: {
    double tiny_precursor = precursor_window_ * 1e-6;
    for (vector<int>::const_iterator ie = negative_isotope_errors_.begin(); ie != negative_isotope_errors_.end(); ++ie) {
      out_min->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 - tiny_precursor));
      out_max->push_back((sc.neutral_mass + (*ie * unit_dalton)) * (1.0 + tiny_precursor));
    }
    *min_range = (sc.neutral_mass + (negative_isotope_errors_.front() * unit_dalton)) * (1.0 - tiny_precursor);
    *max_range = (sc.neutral_mass + (negative_isotope_errors_.back() * unit_dalton)) * (1.0 + tiny_precursor);
    break;
  }
  default:
    carp(CARP_FATAL, "Invalid window type");
  }
  carp(CARP_DETAILED_DEBUG, "Scan=%d Charge=%d Mass window=[%f, %f]",
       sc.spectrum->SpectrumNumber(), sc.charge, (*out_min)[0], (*out_max)[0]);
}

vector<int> TideLiteSearchApplication::getNegativeIsotopeErrors() {
  string isotope_errors_string = Params::GetString("isotope-error");
  if (isotope_errors_string[0] == ',') {
    carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s)",
         isotope_errors_string.c_str());
  }
  for (string::const_iterator i = isotope_errors_string.begin(); i != isotope_errors_string.end(); i++) {
    if (*i == ',' && (i+1 == isotope_errors_string.end() || *(i+1) == ',')) {
      carp(CARP_FATAL, "Error in isotope_error parameter formatting: (%s) ", isotope_errors_string.c_str());
    }
  }

  vector<int> negative_isotope_errors(1, 0);
  if (!isotope_errors_string.empty()) {
    vector<int> isotope_errors = StringUtils::Split<int>(isotope_errors_string, ',');
    for (vector<int>::iterator i = isotope_errors.begin(); i != isotope_errors.end(); i++) {
      if (*i < 0) {
        carp(CARP_FATAL, "Found a negative isotope error: %d.", *i);
      } else if (find(negative_isotope_errors.begin(), negative_isotope_errors.end(),
                      -(*i)) != negative_isotope_errors.end()) {
        carp(CARP_FATAL, "Found duplicate when parsing isotope_error parameter: %d", *i);
      }
      negative_isotope_errors.push_back(-(*i));
    }
  }
  sort(negative_isotope_errors.begin(), negative_isotope_errors.end());
  return negative_isotope_errors;
}


void TideLiteSearchApplication::getPeptideIndexData(const string input_index, ProteinVec& proteins, vector<const pb::AuxLocation*>& locations, pb::Header& peptides_header){

  string peptides_file = FileUtils::Join(input_index, "pepix");  
  string proteins_file = FileUtils::Join(input_index, "protix");
  string auxlocs_file = FileUtils::Join(input_index, "auxlocs");  

  // Read protein index file
  carp(CARP_INFO, "Reading index %s", input_index.c_str());

  pb::Header protein_header;

  if (!ReadRecordsToVector<pb::Protein, const pb::Protein>(&proteins, proteins_file, &protein_header)) {
    carp(CARP_FATAL, "Error reading index (%s)", proteins_file.c_str());
  }
  // There shouldn't be more than one header in the protein pb.
  pb::Header_Source headerSource = protein_header.source(0);  
  string decoy_prefix = "";
  if (headerSource.has_decoy_prefix()){
    decoy_prefix = headerSource.decoy_prefix();
  } else {
    carp(CARP_WARNING, "You are using an index generated by an old version of tide-index."
                       "This will not affect your results, but this index may need to be "
                       "re-created to work with future versions of tide-index. ");
  }
  
  // Read auxlocs index file
  ReadRecordsToVector<pb::AuxLocation>(&locations, auxlocs_file);

  // Read peptides index file

  if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
      !peptides_header.has_peptides_header()) {
    carp(CARP_FATAL, "Error reading index (%s)", peptides_file.c_str());
  }

  const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();
  decoy_num_ = pepHeader.has_decoys_per_target() ? pepHeader.decoys_per_target() : 0;

  // Initizalize the Mass Constants class
  MassConstants::Init(&pepHeader.mods(), &pepHeader.nterm_mods(), &pepHeader.cterm_mods(),
      &pepHeader.nprotterm_mods(), &pepHeader.cprotterm_mods(), bin_width_, bin_offset_);

  // Calculate the Amino Acid Frequencies for the P-value calculation
  // pb::Peptide current_pb_peptide_;
  // HeadedRecordReader peptide_reader = HeadedRecordReader(peptides_file, &peptides_header);
  // while (!(peptide_reader.Done())) { //read all peptides form index
  //   peptide_reader.Read(&current_pb_peptide_);
  //   //Peptide* peptide = new(&fifo_alloc_peptides_) Peptide(current_pb_peptide_, proteins_, &fifo_alloc_peptides_);
  //   PeptideLite* peptide = new PeptideLite(current_pb_peptide_, proteins, NULL);
  //   string proteinNames = peptide->GetLocationStr("decoy_");
  //   string flankingAAs = peptide->GetFlankingAAs();
  //   string peptide_with_mods = peptide->SeqWithMods(4);
  //   string modifications = peptide->getModifications(4); // todo: needs to be tested
  //   printf("%s\t", peptide_with_mods.c_str());
  //   printf("%s\t", modifications.c_str());
  //   printf("%s\t", flankingAAs.c_str());
  //   printf("%s\n", proteinNames.c_str());
  // }
  // ModificationDefinition::ClearAll();

  // const pb::ModTable *modtable = &pepHeader.mods();
  // int stat_mod_size = modtable->static_mod_size();
  // printf("static mod size:%d", stat_mod_size);
  // for (int i = 0; i < modtable->static_mod_size(); ++i) {
  //   char aa = modtable->static_mod(i).amino_acids()[0];
  //   double delta = modtable->static_mod(i).delta();
  //   printf("aa:%c, mass:%lf\n", aa, delta);
  // }

  // modtable = &pepHeader.nterm_mods();
  // stat_mod_size = modtable->static_mod_size();
  // printf("static mod size:%d", stat_mod_size);
  // for (int i = 0; i < modtable->static_mod_size(); ++i) {
  //   char aa = modtable->static_mod(i).amino_acids()[0];
  //   double delta = modtable->static_mod(i).delta();
  //   printf("aa:%c, mass:%lf\n", aa, delta);
  // }


  // Initialize the TideLiteMatchSet for reporting search results
/*  TideMatchSet::decoy_prefix_ = decoy_prefix;  

  stringstream ss;
  ss << Params::GetString("enzyme") << '-' << Params::GetString("digestion");
  TideMatchSet::CleavageType = ss.str();

  TideMatchSet::initModMap(pepHeader.mods(), ANY);
  TideMatchSet::initModMap(pepHeader.nterm_mods(), PEPTIDE_N);
  TideMatchSet::initModMap(pepHeader.cterm_mods(), PEPTIDE_C);
  TideMatchSet::initModMap(pepHeader.nprotterm_mods(), PROTEIN_N);
  TideMatchSet::initModMap(pepHeader.cprotterm_mods(), PROTEIN_C);
*/
  //TODO: do the amino acid frequency calculation here. Shouldn't it be done in tide-index?

}

void TideLiteSearchApplication::createOutputFiles(ofstream **target_file, ofstream **decoy_file) {
  
  // Create output search results files
  bool compute_sp = false;   //TODO: reconsider later    bool compute_sp = Params::GetBool("compute-sp");
  
  string concat_file_name = make_file_path("tide-search.txt");
  string target_file_name = make_file_path("tide-search.target.txt");
  string decoy_file_name  = make_file_path("tide-search.decoy.txt");
  bool overwrite = Params::GetBool("overwrite");  
  if (overwrite) {
    remove(concat_file_name.c_str());  
    remove(target_file_name.c_str());  
    remove(decoy_file_name.c_str());  
  }
  if (Params::GetBool("concat")) {

    *target_file = create_stream_in_path(concat_file_name.c_str(), NULL, overwrite);
    output_file_name_ = concat_file_name;
  
  } else {
  
    *target_file = create_stream_in_path(target_file_name.c_str(), NULL, overwrite);
    output_file_name_ = target_file_name;
    if (decoy_num_ > 0) {
      *decoy_file = create_stream_in_path(decoy_file_name.c_str(), NULL, overwrite);
    }
  }  
  if (*target_file) {
    TideMatchSet::writeHeaders(*target_file, false, decoy_num_ > 1, compute_sp);
    TideMatchSet::writeHeaders(*decoy_file, true, decoy_num_ > 1, compute_sp);
  }

}

string TideLiteSearchApplication::getOutputFileName() {
  return output_file_name_;
}

// In order to add more options, you need to add them to ./src/util/Params.cpp
vector<string> TideLiteSearchApplication::getOptions() const {
  string arr[] = {
    "auto-mz-bin-width",
    "auto-precursor-window",
//    "compute-sp",   // To be removed
    "concat",
    "deisotope",
    "elution-window-size",
    "evidence-granularity", // To be removed
    "file-column",
    "fileroot",
    "fragment-tolerance",
    "isotope-error",
    "mass-precision",
    "max-precursor-charge",
    "min-precursor-charge",
    "min-peaks",
    "mod-precision",
    "mz-bin-offset",
    "mz-bin-width",
    "mzid-output",
    "num-threads",
    "output-dir",
    "override-charges",
    "overwrite",
    "parameter-file",
//    "peptide-centric-search", // To be removed
    "pepxml-output",
    "pin-output",
    "pm-charges",
    "pm-max-frag-mz",
    "pm-max-precursor-delta-ppm",
    "pm-max-precursor-mz",
    "pm-max-scan-separation",
    "pm-min-common-frag-peaks",
    "pm-min-frag-mz",
    "pm-min-peak-pairs",
    "pm-min-precursor-mz",
    "pm-min-scan-frag-peaks",
    "pm-pair-top-n-frag-peaks",
    "pm-top-n-frag-peaks",
    "precision",
    "precursor-window",
    "precursor-window-type",
    "print-search-progress",
    "remove-precursor-peak",
    "remove-precursor-tolerance",
    "scan-number",
    "score-function",
    "skip-preprocessing",
    "spectrum-max-mz",
    "spectrum-min-mz",
    "spectrum-parser",
    "sqt-output",
    "store-index",
    "store-spectra",
    "top-match",
    "txt-output",
    "use-flanking-peaks",
    "use-neutral-loss-peaks",
    "use-z-line",
    "verbosity"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}
/*
Parameters to be removed from ./Params.cpp"
brief-output
peptide-centric-search
use-tailor-calibration
exact-p-value
charge-state
*/

string TideLiteSearchApplication::getName() const {
  return "tide-lite-search";
}

string TideLiteSearchApplication::getDescription() const {
  return
    "[[nohtml:Search a collection of spectra against a sequence database, "
    "returning a collection of peptide-spectrum matches (PSMs). This is a "
    "fast search engine but requires that you first build an index with "
    "tide-index.]]"
    "[[html:<p>Tide is a tool for identifying peptides from tandem mass "
    "spectra. It is an independent reimplementation of the SEQUEST<sup>&reg;"
    "</sup> algorithm, which assigns peptides to spectra by comparing the "
    "observed spectra to a catalog of theoretical spectra derived from a "
    "database of known proteins. Tide's primary advantage is its speed. Our "
    "published paper provides more detail on how Tide works. If you use Tide "
    "in your research, please cite:</p><blockquote>Benjamin J. Diament and "
    "William Stafford Noble. <a href=\"http://dx.doi.org/10.1021/pr101196n\">"
    "&quot;Faster SEQUEST Searching for Peptide Identification from Tandem "
    "Mass Spectra&quot;</a>. <em>Journal of Proteome Research</em>. "
    "10(9):3871-9, 2011.</blockquote> "
    "<p>When <code>tide-search</code> runs, it performs "
    "several intermediate steps, as follows:</p><ol>"
    "<li>If a FASTA file was provided, convert it to an index using "
    "<code>tide-index</code>.</li>"
    "<li>Convert the given "
    "fragmentation spectra to a binary format.</li><li>Search the spectra "
    "against the database and store the results in binary format.</li><li>"
    "Convert the results to one or more requested output formats.</li></ol><p>"
    "By default, the intermediate binary files are stored in the output "
    "directory and deleted when Tide finishes execution. If you plan to search "
    "against given database more than once or search a given set of spectra "
    "more than once, then you can direct Tide to save the binary spectrum "
    "files using the <code>--store-index</code> and "
    "<code>--store-spectra</code> options. "
    "Subsequent runs of the program will go faster "
    "if provided with inputs in binary format.</p>]]";
}

vector<string> TideLiteSearchApplication::getArgs() const {
  string arr[] = {
    "tide spectra file+",
    "tide database"
  };
  return vector<string>(arr, arr + sizeof(arr) / sizeof(string));
}


vector< pair<string, string> > TideLiteSearchApplication::getOutputs() const {
  vector< pair<string, string> > outputs;
  outputs.push_back(make_pair("tide-lite-search.target.txt",
    "a tab-delimited text file containing the target PSMs. See <a href=\""
    "../file-formats/txt-format.html\">txt file format</a> for a list of the fields."));
  outputs.push_back(make_pair("tide-lite-search.decoy.txt",
    "a tab-delimited text file containing the decoy PSMs. This file will only "
    "be created if the index was created with decoys."));
  outputs.push_back(make_pair("tide-lite-search.params.txt",
    "a file containing the name and value of all parameters/options for the "
    "current operation. Not all parameters in the file may have been used in "
    "the operation. The resulting file can be used with the --parameter-file "
    "option for other Crux programs."));
  outputs.push_back(make_pair("tide-lite-search.log.txt",
    "a log file containing a copy of all messages that were printed to the "
    "screen during execution."));
  return outputs;
}
bool TideLiteSearchApplication::needsOutputDirectory() const {
  return true;
}

COMMAND_T TideLiteSearchApplication::getCommand() const {
  return TIDE_LITE_SEARCH_COMMAND;
}

void TideLiteSearchApplication::processParams() {
  const string index = Params::GetString("tide database");

  if (!FileUtils::Exists(index)) {
    carp(CARP_FATAL, "'%s' does not exist", index.c_str());
  
  } else if (FileUtils::IsRegularFile(index)) {
    // Index is FASTA file

    carp(CARP_INFO, "Creating index from '%s'", index.c_str());
    string targetIndexName = Params::GetString("store-index");
    if (targetIndexName.empty()) {
      targetIndexName = FileUtils::Join(Params::GetString("output-dir"),
                                        "tide-search.tempindex");
      remove_index_ = targetIndexName;
    }
    TideIndexApplication indexApp;
    indexApp.processParams();
    if (indexApp.main(index, targetIndexName) != 0) {
      carp(CARP_FATAL, "tide-index failed.");
    }
    Params::Set("tide database", targetIndexName);
  
  } else {
    // Index is Tide index directory

    pb::Header peptides_header;
    string peptides_file = FileUtils::Join(index, "pepix");
    HeadedRecordReader peptide_reader(peptides_file, &peptides_header);
    if ((peptides_header.file_type() != pb::Header::PEPTIDES) ||
        !peptides_header.has_peptides_header()) {
      carp(CARP_FATAL, "Error reading index (%s).", peptides_file.c_str());
    }

    const pb::Header::PeptidesHeader& pepHeader = peptides_header.peptides_header();

    Params::Set("enzyme", pepHeader.enzyme());
    const char* digestString =
      digest_type_to_string(pepHeader.full_digestion() ? FULL_DIGEST : PARTIAL_DIGEST);
    Params::Set("digestion", digestString);
    Params::Set("isotopic-mass", pepHeader.monoisotopic_precursor() ? "mono" : "average");
  }
  // run param-medic?
/*  const string autoPrecursor = Params::GetString("auto-precursor-window");
  const string autoFragment = Params::GetString("auto-mz-bin-width");
  if (autoPrecursor != "false" || autoFragment != "false") {
    if (autoPrecursor != "false" && Params::GetString("precursor-window-type") != "ppm") {
      carp(CARP_FATAL, "Automatic peptide mass tolerance detection is only supported with ppm "
                       "units. Please re-run with auto-precursor-window set to 'false' or "
                       "precursor-window-type set to 'ppm'.");
    }
    ParamMedic::RunAttributeResult errorCalcResult;
    ParamMedicApplication::processFiles(Params::GetStrings("tide spectra file"),
      true, false, &errorCalcResult, NULL);

    if (autoPrecursor != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_PRECURSOR_PREDICTION));
        carp(CARP_INFO, "precursor ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "precursor error estimate (ppm): %.2f", prediction);
        Params::Set("precursor-window", prediction);
      } else {
        carp(autoPrecursor == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate precursor error: %s", fail.c_str());
      }
    }
    if (autoFragment != "false") {
      string fail = errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_FAILURE);
      if (fail.empty()) {
        double sigma = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_SIGMA));
        double prediction = StringUtils::FromString<double>(
          errorCalcResult.getValue(ParamMedic::ErrorCalc::KEY_FRAGMENT_PREDICTION));
        carp(CARP_INFO, "fragment ppm standard deviation: %f", sigma);
        carp(CARP_INFO, "Fragment bin size estimate (Th): %.4f", prediction);
        Params::Set("mz-bin-width", prediction);
      } else {
        carp(autoFragment == "fail" ? CARP_FATAL : CARP_ERROR,
             "failed to calculate fragment error: %s", fail.c_str());
      }
    }
  }  
  */
}

vector<TideLiteSearchApplication::InputFile> TideLiteSearchApplication::getInputFiles(const vector<string>& filepaths) const {
  // Try to read all spectrum files as spectrumrecords, convert those that fail
  vector<TideLiteSearchApplication::InputFile> input_sr;
  for (vector<string>::const_iterator f = filepaths.begin(); f != filepaths.end(); f++) {
    SpectrumCollection spectra;
    pb::Header spectrum_header;
    string spectrumrecords = *f;
    bool keepSpectrumrecords = true;
    if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
      // Failed, try converting to spectrumrecords file
      carp(CARP_INFO, "Converting %s to spectrumrecords format", f->c_str());
      carp(CARP_INFO, "Elapsed time starting conversion: %.3g s", wall_clock() / 1e6);
      spectrumrecords = Params::GetString("store-spectra");
      keepSpectrumrecords = !spectrumrecords.empty();
      if (!keepSpectrumrecords) {
        spectrumrecords = make_file_path(FileUtils::BaseName(*f) + ".spectrumrecords.tmp");
      } else if (filepaths.size() > 1) {
        carp(CARP_FATAL, "Cannot use store-spectra option with multiple input "
                         "spectrum files");
      }
      carp(CARP_DEBUG, "New spectrumrecords filename: %s", spectrumrecords.c_str());
      if (!SpectrumRecordWriter::convert(*f, spectrumrecords)) {
        carp(CARP_FATAL, "Error converting %s to spectrumrecords format", f->c_str());
      }
      carp(CARP_DEBUG, "Reading converted spectrum file %s", spectrumrecords.c_str());
      // Re-read converted file as spectrumrecords file
      if (!spectra.ReadSpectrumRecords(spectrumrecords, &spectrum_header)) {
        carp(CARP_DEBUG, "Deleting %s", spectrumrecords.c_str());
        FileUtils::Remove(spectrumrecords);
        carp(CARP_FATAL, "Error reading spectra file %s", spectrumrecords.c_str());
      }
    }
    input_sr.push_back(TideLiteSearchApplication::InputFile(*f, spectrumrecords, keepSpectrumrecords));
  }
  return input_sr;
}

SpectrumCollection* TideLiteSearchApplication::loadSpectra(const string& file) {

}
