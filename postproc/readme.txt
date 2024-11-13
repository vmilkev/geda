
1. The execution script (that calls a required program) is here: /usr/home/qgg/vimi/gedaplus/run_gedaplus.sh
2. The script can be copied to a convenient location and launched from there.
3. How to run/launch the script (three variants):
   i. path_to_script/run_gedaplus.sh geda_resulting_file
      Example: gedaplus/run_gedaplus.sh gedaplus/aligned_init_sniffer_101.txt
      NOTE: this variant should be used only if entire data in the input file is reliable.
   ii. path_to_script/run_gedaplus.sh geda_resulting_file integer_number
       Where the integer_number is a window size in hours for baseline estimation/correction; the default is 6 hours;
       Example: ./run_gedaplus.sh aligned_init_sniffer_101.txt 4
       NOTE: this variant should be used only if entire data in the input file is reliable.
   iii. path_to_script/run_gedaplus.sh geda_resulting_file integer_number file_with_reliable_signals
        Where file_with_reliable_signals is a file with signals IDs for which phenotypes need to be calculated;
        Example: ./run_gedaplus.sh aligned_init_sniffer_4.txt 3 ids_list_sniffer_4.txt
        NOTE: this variant should be used if only the part of the data in the input file is reliable; the reliable signals are provided in the file_with_reliable_signals.
4. Where to find files examples: /usr/home/qgg/vimi/gedaplus/
5. Output is a file consisting of time stamp, animal id and 4 variants of traits; each gas results are in a separate file.
   NOTE: units of calculated traits are vapue_per_minute!