#ifndef TIDELITESEARCHAPPLICATION_H
#define TIDELITESEARCHAPPLICATION_H

#include "CruxApplication.h"
#include "TideMatchSet.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <gflags/gflags.h>
#include "peptides.pb.h"
#include "spectrum.pb.h"
#include "tide/theoretical_peak_set.h"
#include "tide/max_mz.h"
#include "util/MathUtil.h"
#include "tide/ActivePeptideQueueLite.h"
#include "TideIndexApplication.h"
#include "TideLiteMatchSet.h"

using namespace std;

class TideLiteSearchApplication : public CruxApplication {
 private:
  struct InputFile {
    std::string OriginalName;
    std::string SpectrumRecords;
    bool Keep;
    InputFile(const std::string& name,
              const std::string& spectrumrecords,
              bool keep):
      OriginalName(name), SpectrumRecords(spectrumrecords), Keep(keep) {}
  };
 protected:
  static const double XCORR_SCALING;
  static const double RESCALE_FACTOR;
  static const double TAILOR_QUANTILE_TH;
  static const double TAILOR_OFFSET;

  map<pair<string, unsigned int>, bool>* spectrum_flag_;
  string output_file_name_;
  std::string remove_index_;  
  double bin_width_;
  double bin_offset_;
  bool use_neutral_loss_peaks_;
  bool use_flanking_peaks_;

  long int num_range_skipped_;
  long int num_precursors_skipped_;
  long int num_isotopes_skipped_;
  long int num_retained_;
  long int total_candidate_peptides_;
  int print_interval_;
  vector<int> negative_isotope_errors_;
  WINDOW_TYPE_T window_type_;
  double precursor_window_;
  double spectrum_min_mz_;
  double spectrum_max_mz_;
  double min_scan_;
  double max_scan_;
  double min_peaks_;
  double min_precursor_charge_;
  double max_precursor_charge_;
  int num_threads_;

  bool tsv_output_;  // original tide-search output format in tab-delimited text files (txt)
  bool mztab_output_;// mzTAB output format
  bool pin_output_;  // pin output format for percolator.

  ofstream* target_file_;
  ofstream* decoy_file_;

  vector<boost::mutex *> locks_array_;  

  vector<TideLiteSearchApplication::InputFile> getInputFiles(const vector<string>& filepaths) const;
  static SpectrumCollection* loadSpectra(const std::string& file);
  void getPeptideIndexData(string, ProteinVec& proteins, vector<const pb::AuxLocation*>& locations, pb::Header& peptides_header);
  void createOutputFiles(std::ofstream **targe_file, std::ofstream **decoy_file);
  vector<int> getNegativeIsotopeErrors();

  /**
   * Function that contains the search algorithm and performs the search
   */
  void search();
  void XCorrScoring(SpectrumCollection::SpecCharge* sc, TideLiteMatchSet& psm_scores);
  void PValueScoring();

  void computeWindow(
      const SpectrumCollection::SpecCharge& sc,
      vector<double>* out_min,
      vector<double>* out_max,
      double* min_range,
      double* max_range
    );
 public:

  SCORE_FUNCTION_T curScoreFunction_;
  ActivePeptideQueueLite* active_peptide_queue_; 
  int decoy_num_;  // Number of decoys per peptide;
  int top_matches_;

  vector<double> dAAFreqN_;
  vector<double> dAAFreqI_;
  vector<double> dAAFreqC_;
  vector<double> dAAMass_;
  map<double, std::string> mMass2AA_;


  /**
   * Constructor
   */
  TideLiteSearchApplication();

  /**
   * Destructor
   */
  ~TideLiteSearchApplication();

  /**
   * Main methods
   */
  virtual int main(int argc, char** argv);

  int main(const vector<string>& input_files);

  int main(const vector<string>& input_files, const string input_index);

  /**
   * Returns the command name
   */
  virtual string getName() const;

  /**
   * Returns the command description
   */
  virtual string getDescription() const;

  /**
   * Returns the command arguments
   */
  virtual vector<string> getArgs() const;

  /**
   * Returns the command options
   */
  virtual vector<string> getOptions() const;

  /**
   * Returns the command outputs
   */
  virtual vector< pair<string, string> > getOutputs() const;

  /**
   * Returns whether the application needs the output directory or not. (default false)
   */
  virtual bool needsOutputDirectory() const;

  /**
   * Returns the command ID 
   */
  virtual COMMAND_T getCommand() const;

  /**
   * Processes the output file names
   */
  string getOutputFileName();

  /**
   * Processes the parameters
   */
  virtual void processParams();

};

#endif

/*
 * Local Variables:
 * mode: c
 * c-basic-offset: 2
 * End:
 */
