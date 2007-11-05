#!PYTHON
# FILE: compare_runtime.py
# AUTHOR: CHRIS PARK
# CREATE DATE: 7/30/2007
# SEQUEST="./sequest-one-xcorr"
SEQUEST="./sequest27"
TIME_PREFIX="real "

"""
This script compares the runtime between Crux and Sequest
"""
import os
import sys
import commands
from optparse import OptionParser

#-------------------

def plot_compare_data(crux_array, crux_no_index_array, sequest_array,
    mass_windows, number_of_spectra, fasta_file):
  """compares runtime for each scoring method """

  fasta_file = os.path.basename(fasta_file)
  print "generating figures"
  
  fh = open("crux-fast-%s.xy" % fasta_file, "w")
  for idx in range(len(mass_windows)):
    fh.write("%.8f\t%.8f\n" % (mass_windows[idx], crux_array[idx]))
  fh.close()

  fh = open("crux-no-index-%s.xy" % fasta_file, "w")
  for idx in range(len(mass_windows)):
    fh.write("%.8f\t%.8f\n" % (mass_windows[idx], crux_no_index_array[idx]))
  fh.close()

  fh = open("sequest-full-xcorr-%s.xy" % fasta_file, "w")
  for idx in range(len(mass_windows)):
    fh.write("%.8f\t%.8f\n" % (mass_windows[idx], sequest_array[idx]))
  fh.close()


#-------------------------------------------
# By, Christopher Park
# This script compares the runtime between of Sequest and CRUX
#-------------------------------------------

# Process command line options
usage = "Usage: compare_runtime <ms2 file> <fasta_file>"
option_parser = OptionParser(usage)
(options, args) = option_parser.parse_args()

print "number of args: %d\n" % len(args)

if not len(args) == 2:
  print usage
  sys.exit(1)

ms2_file = args[0]
fasta_file = args[1]

#mass windows to test runtime
mass_windows = ["0.1", "1", "3"]
# mass_windows = ["0.1"]

#runtime result arrays for each method
sequest_results = []
crux_results = []
crux_no_index_results = []

number_of_spectra = 0

charge_list = [1,2,3]

#first run Sequest with varying mass windows
for window in mass_windows:

  # 1, run Sequest
  seq_time = 0.0
  for charge in charge_list:
    command = "time -p %s -D%s -Psequest.params_%s *.%i.dta" \
        % (SEQUEST, fasta_file, window, charge)
    print >>sys.stderr, "\nRunning %s\n" % command
    (exit_code, result) = commands.getstatusoutput(command)
    #debug
    print result
    #print exit_code
  
    if exit_code == "1":
      print "%s %s" % ("failed to run Sequest on mass window:", window)
      sys.exit(1)
    else:
      #now parse the runtime from the result output
      result = result.split('\n')
      for line in result:
        #get user time
        #FIXME is it user or real?
        if line.startswith(TIME_PREFIX):
          fields = line. rstrip('\n').split()
          seq_time += float(fields[1])
  sequest_results.append(seq_time)

        
  # 2, now run Crux
  command = "time -p ./match_search \
      --output-mode sqt \
      --sqt-output-file %s.sqt \
      --parameter-file crux.params_%s \
      --number-decoy-set 0 \
      %s %s" % (window, window, ms2_file, fasta_file)
  print >>sys.stderr, "\nRunning %s\n" % command

  (exit_code, result) = \
        commands.getstatusoutput(command)
  #debug
  print result
  #print exit_code
  
  if exit_code == "1":
    print "%s %s" % ("failed to run Crux on mass window:", window)
    sys.exit(1)
  else:
    #now parse the runtime from the result output
    result = result.split('\n')
    for line in result:
      #get user time
      if line.startswith(TIME_PREFIX):
        fields = line.rstrip('\n').split()
        crux_results.append(float(fields[1]))

  # 3, now run Crux without an index
  command = "time -p ./match_search \
      --output-mode sqt \
      --sqt-output-file %s.sqt \
      --parameter-file crux_no_index.params_%s \
      --number-decoy-set 0 \
      %s %s" % (window, window, ms2_file, fasta_file)
  print >>sys.stderr, "\nRunning %s\n" % command

  (exit_code, result) = \
        commands.getstatusoutput(command)
  #debug
  print result
  #print exit_code
  
  if exit_code == "1":
    print "%s %s" % ("failed to run Crux on mass window:", window)
    sys.exit(1)
  else:
    #now parse the runtime from the result output
    result = result.split('\n')
    for line in result:
      #get user time
      if line.startswith(TIME_PREFIX):
        fields = line.rstrip('\n').split()
        crux_no_index_results.append(float(fields[1]))



#Debug
#print crux_results, sequest_results, mass_windows

# calculate number of spectra
dtas = filter(lambda x: x.endswith("dta"), os.listdir(os.getcwd()))


#now plot the results
plot_compare_data(crux_results, crux_no_index_results, sequest_results, [ float(i) for i in mass_windows], len(dtas), fasta_file)
