V0.32

PARSES: A Pipeline for Analysis of RNA-Sequencing Exogenous Sequences

Step by step instructions for OS X (exclude `backquotes` when typing a command):
1. Open a terminal.
2. Execute `which gcc`, if /usr/bin/gcc is returned then continue. Otherwise, you have to download and install Xcode. If you have Leopard download Xcode from here: http://adcdownload.apple.com/Developer_Tools/xcode_3.1.3_developer_tools/xcode313_2736_developerdvd.dmg. If you have Snow Leopard download Xcode from here: http://adcdownload.apple.com/ios/ios_sdk_4.1__final/xcode_3.2.4_and_ios_sdk_4.1.dmg. Simply open the dmg file when it is finished downloading and agree to everything with the default settings.
3. Execute `which rake`, if a directory is returned skip to step 9.
4. Open a browser and navigate to 'http://rubyforge.org/frs/download.php/56872/rake-0.8.7.tgz', save the file to your home directory.
5. Open a terminal in your home directory.
6. Execute `tar -xzf rake-0.8.7.tgz`
7. Execute `cd rake-0.8.7`
8. Execute `sudo ruby install.rb`
9. Extract PARSES into your ~/bin folder. (make the directory if it does not exist. You can do this through the GUI)
10. Open a terminal in ~/bin.
11. Execute `chmod +x *.pl`
12. Change directory to directory containing FASTQ sequence (you may want to make a backup of your data, my program shouldn't do anything to it but I would hate to be responsible for such a thing.). For instance if the data is contained in /usr/share execute `cd /usr/share`
13. Execute `sudo rake -f ~/bin/rakefile.rb seq=mysequence file=reads.fastq type=solexa install`. Please note your sequence name should not have any spaces. Keep all of the default settings for the MEGAN install. If any of the installations fail manually install them.
14. Execute `rake -f ~/bin/rakefile.rb seq=mysequence file=reads.fastq type=solexa`.
15. Register MEGAN if you have not yet. The program should open and you can follow the instructions to get emailed a key.
16. Update the LCA Paramaters by selecting Options->Change LCA Parameters in MEGAN.
17. Right click the nodes of each of the taxonomies individually you are interested in and select 'Extract Reads by Taxa' for each. Keep the default name. Exit MEGAN.
18. If the pipeline continues running, good, otherwise, execute `rake -f ~/bin/rakefile.rb seq=mysequence`
19. MEGAN should open a file with all of the reads you extracted de novo aligned and binned into their original taxonomies. You can change a file anywhere in the process and the pipeline will run from that point to the end (or whatever task you select) when you execute `rake -f ~/bin/rakefile.rb seq=mysequence`.

Requirements:
Root priveleges if installing software.
gcc - http://adcdownload.apple.com/ios/ios_sdk_4.1__final/xcode_3.2.4_and_ios_sdk_4.1.dmg (for Mac)
Ruby - http://www.ruby-lang.org/en/downloads/
Rake - http://rake.rubyforge.org/
abyssKmerOptimizer.pl marked as executable in same folder as rakefile
fac.pl marked as executable in same folder as rakefile
Xextractspans.pl marked as executable in same folder as rakefile
Xfilterspans.pl marked as executable in same folder as rakefile
Xnovotonm.pl marked as executable in same folder as rakefile
Working directory is the directory in which the fasta sequence file is contained.

If you want to get all applications used in the pipeline installed upfront example:
sudo rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fasta type=illumina1.3 install

First run example:
rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fasta type=illumina1.3

Subsequent runs example:
rake -f /rake/file/location seq=NameYouGiveToYourSequence

Run up to a certain point (stop after executing a specified task):
rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fasta type=illumina1.3 taskName

To get all task names (list of all tasks are below):
rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fasta type=illumina1.3 -T

rake abyssInstall              # Install latest version of ABySS with Googl...
rake alignSequence             # Novoalign - Align reads in to base genome
rake blastInstall              # Install latest version of BLAST+.
rake bowtieIndex               # Create an index for bowtie/tophat of the h...
rake bowtieInstall             # Install latest version of Bowtie.
rake denovoAssembleCluster     # ABySS - Assemble reads associated with clu...
rake hgInstall                 # Install latest version of human genome dat...
rake install                   # Install latest version of everything.
rake localAlignContigs         # BLAST - Associate contigs with organisms.
rake localAlignReads           # BLAST - Associate reads with organisms.
rake meganInstall              # Install latest version of MEGAN.
rake metaGenomeAnalyzeContigs  # MEGAN - Separate contigs into taxonomies.
rake metaGenomeAnalyzeReads    # MEGAN - Separate reads into taxonomies.
rake novoIndex                 # Create an index for novoalign of the human...
rake novoalignInstall          # Install latest version of Novoalign.
rake ntInstall                 # Install latest version of the NT database.
rake parallelIteratorInstall   # Install latest version of Parallel::Iterat...
rake removeHuman               # Xnovotonm - Harvest non-base organism reads
rake removeSpans               # Tophat - Align spanning reads in order to ...
rake reserialize               # Automatically saving any settings changes ...
rake samtoolsInstall           # Install latest version of Samtools.
rake tophatInstall             # Install latest version of Tophat.

Most parameters are set by default. This can be changed in the .sequenceName file which is produced automatically from the seq argument specificed and stored in your home directory. Editing this file will change the behavior of the program. It is highly recommended that none of the path strings are changed. PARSES does not yet have the ability to determine where execution should begin based on parameters reset. You can manually force PARSES to reprocess a file by executing the touch command, `touch filename`, on the file you wish to be reprocessed with the new settings.

You may be interested in viewing the log file produced for your sequence called sequenceName.log.

PARSES will find any human genome databases you have as long as the files have the form chr*.fa. In addition, any novoalign index will be located and used if it has the form '*hg*.ndx' as well as any bowtie index with the form '*hg*.ebwt' and any NT database if its location is specified by the BLASTDB environment database. Any database or index not found will be downloaded/built.
