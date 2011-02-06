#PARSES V0.41
For more information see [wiki](https://github.com/Lythimus/PARSES/wiki)

#Requirements
PARSES is intended to be executed with Solexa data on a *NIX-based desktop computer. It may not be used for any sort of financial gain. Licenses for both MEGAN and Novoalign are strictly for non-profit research at non-profit institutions and academic usage.

##Supported Data Types
* fasta
* sanger
* solexa
* illumina1.3
* illumina1.5

##System Requirements
PARSES's system requirements are directly dependent on the size of the data set being processed. It is recommended to be run on a Linux or OS X 64-bit machine with at least 4GBs of memory. You must also have root privileges to the machine if you are installing software.

##Software Requirements
* 
* Working directory is the directory in which the FASTQ sequence file is contained.
* [gcc](http://adcdownload.apple.com/ios/ios_sdk_4.1__final/xcode_3.2.4_and_ios_sdk_4.1.dmg)
* [Ruby](http://www.ruby-lang.org/en/downloads)
* [Rake](http://rake.rubyforge.org)
* abyssKmerOptimizer.pl marked as executable in same folder as rakefile
* fac.pl marked as executable in same folder as rakefile
* addTaxon.pl marked as executable in same folder as rakefile
* parallelBlast.sh marked as executable in same folder as rakefile
* Xextractspans.pl marked as executable in same folder as rakefile
* Xfilterspans.pl marked as executable in same folder as rakefile
* Xnovotonm.pl marked as executable in same folder as rakefile
___

#Installation
Execute all commands from the directory of the data. Place all scripts for PARSES into a single directory and mark as executable. All latest versions of programs will be installed, in the event of an error during installation a repository of the programs may be used which is not guaranteed to be up to date. The repository can be activated by including the repo=true command. Installation will automatically be performed during any execution but it can be manually performed by evoking any of the following installation commands:

	sudo rake -f /rake/file/location install #installs and indexes all resources.
	sudo rake -f /rake/file/location novoalignInstall
	sudo rake -f /rake/file/location bowtieInstall
	sudo rake -f /rake/file/location hgInstall
	sudo rake -f /rake/file/location novoIndex
	sudo rake -f /rake/file/location bowtieIndex
	sudo rake -f /rake/file/location samtoolsInstall
	sudo rake -f /rake/file/location tophatInstall
	sudo rake -f /rake/file/location abyssInstall
	sudo rake -f /rake/file/location blastInstall
	sudo rake -f /rake/file/location ntInstall
	sudo rake -f /rake/file/location meganInstall
	sudo rake -f /rake/file/location parallelIteratorInstall

___

#Examples
Example executions.

##First Execution
`rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fastq type=illumina1.3`

##Subsequent Executions
`rake -f /rake/file/location seq=NameYouGiveToYourSequence`

##Install using repository of links
`rake -f /rake/file/location repo=true install`

##Run to Specified Point
`rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fastq type=illumina1.3 localAlignContigs`

##Run truncated version of PARSES (only execute the specified task and ignore prerequisites)
`rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourSequenceFileName.fastq type=illumina1.3 truncate=true localAlignContigs`

##Run truncated version of PARSES (only execute the specified task and ignore prerequisites) and override the file naming schema
`rake -f /rake/file/location seq=NameYouGiveToYourSequence file=YourFileToProcess type=illumina1.3 truncate=true forcefile=true localAlignContigs`

##List All Tasks
`rake -f /rake/file/location -T`

##Output Files

datafile.fastq

* alignSequence *(Novoalign)*

datafile.fastq.novo

* removeHuman *(Xnovotonm)*

datafile.fastq.novo.NM.fasta

* removeSpans *(Tophat)*

datafile.fastq.novo.NM.fasta.nospans

* localAlignReads *(BLAST+)*

datafile.fastq.novo.NM.fasta.nospans.blast

* metaGenomeAnalyzeReads *(MEGAN)*

datafile.fastq.novo.NM.fasta.nospans.blast.megan.rma

* denovoAssembleCluster *(ABySS)*

datafile.fastq.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa

* localAlignContigs *(BLAST+)*

datafile.fastq.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast

* metaGenomeAnalyzeContigs *(MEGAN)*

datafile.fastq.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast.megan.rma

___

#Task List

##Basic Tasks
1. alignSequence             # Novoalign - Align reads in to base genome
2. removeHuman               # Xnovotonm - Harvest non-base organism reads
3. removeSpans               # Tophat - Align spanning reads in order to remove base organism reads
4. localAlignReads           # BLAST - Associate reads with organisms.
5. metaGenomeAnalyzeReads    # MEGAN - Separate reads into taxonomies.
6. denovoAssembleCluster     # ABySS - Assemble reads associated with clusters of taxonomies.
7. localAlignContigs         # BLAST - Associate contigs with organisms.
8. metaGenomeAnalyzeContigs  # MEGAN - Separate contigs into taxonomies.

##Special Tasks
* clean                     # Remove any temporary products.
* clobber                   # Remove any generated file.
* reserialize               # Automatically saving any settings changes which may have been made

##Installation Tasks
* abyssInstall              # Install latest version of ABySS with Google Sparsehash.
* blastInstall              # Install latest version of BLAST+.
* bowtieIndex               # Create an index for bowtie/tophat of the human genome database.
* bowtieInstall             # Install latest version of Bowtie.
* hgInstall                 # Install latest version of human genome database.
* install                   # Install latest version of everything.
* meganInstall              # Install latest version of MEGAN.
* novoIndex                 # Create an index for novoalign of the human genome database.
* novoalignInstall          # Install latest version of Novoalign.
* ntInstall                 # Install latest version of the NT database.
* parallelIteratorInstall   # Install latest version of Parallel::Iterator for perl.
* samtoolsInstall           # Install latest version of Samtools.
* tophatInstall             # Install latest version of Tophat.

___

#Analyzing Results

##MEGAN
The primary method for viewing results involves perusing the RMA file produced by MEGAN--though a PDF is also produced with standard LCA parameters. To familiarize yourself with MEGAN you can watch the [MEGAN Introduction Youtube video](http://www.youtube.com/watch?v=i7-OCW-DctY).

Your results will resemble the following:

<img src="http://dl.dropbox.com/u/8709992/minscore55minsupport30.jpg" width="900">

##Report
Reporting is currently stored in logs.

___

#Logging
A log file is produced in the data directory with the filename `sequenceName.log`. It not only logs all commands executed, in addition to errors returned, but it also contains report information on the results of the data analysis. Below is a list of all information logged at each step of PARSES.

##Non-task specific
* Sequence name.
* Read length.
* Data type.
* Timing information for each portion of PARSES.
* Arguments used for each program.

##Initialization
* Information regarding the specifications of the machine on which PARSES is executing. (log)
* Command used to invoke PARSES. (log)
* Total number of reads. (log)

##alignSequence *(Novoalign)*
* Number of reads which align to the human genome. (log)
* Percentage of total reads which align to the human genome. (log)

##removeHuman *(Xnovotonm)*
* Number of reads left after removing human reads. (log)
* Percentage of total reads left. (log)

##removeSpans *(Tophat)*
* Number of reads left after removing splicing junctions of human reads. (log)
* Percentage of reads left after removing splicing regions. (log)

##localAlignReads *(BLAST+)*

##metaGenomeAnalyzeReads *(MEGAN)*
* Potential number of exogenous reads. (log)
* Number of reads assigned to a taxonomy. (log)
* Percentage of total reads assigned to a taxonomy. (log)
* Number of reads not assigned to a taxonomy. (log)
* Percentage of total reads not assigned to a taxonomy. (log)
* Number of reads which did not significantly align. (log)
* Percentage of total reads which did not significantly align. (log)
* Image of taxonomy overview with default LCA parameters.

##denovoAssembleCluster *(ABySS)*
* Coverage threshold. (log)
* Median kmer coverage. (log)
* Number of kmers in contigs. (log)
* Number of contigs. (log)
* Optimum kmer. (log)

##localAlignContigs *(BLAST+)*

##metaGenomeAnalyzeContigs *(MEGAN)*
* Potential number of exogenous contigs. (log)
* Number of contigs assigned to a taxonomy. (log)
* Number of contigs not assigned to a taxonomy. (log)
* Number of contigs which did not significantly align. (log)

Below is an example of a log file:

	# Logfile created on Mon Nov 29 20:21:37 -0600 2010 by logger.rb/22285
	I, [2010-12-27T14:54:16.761320 #94508]  INFO -- : Begin run for seq=akata file=s_4_sequence_Akata.txt type=illumina1.3 task=-f
	I, [2010-12-27T14:54:16.849927 #94508]  INFO -- : Executing PARSES v0.30 in a 0 environment with 64GB of memory and 24 cores with a 64-bit architecture.
	I, [2010-12-27T17:56:20.865114 #94508]  INFO -- : tophat -p 24 --solexa1.3-quals --output-dir akata_tophat_out /usr/share/hgChrAll s_4_sequence_Akata.txt
	I, [2010-12-27T17:56:20.865442 #94508]  INFO -- : samtools view -h -o akata_tophat_out/accepted_hits.sam akata_tophat_out/accepted_hits.bam
	I, [2010-12-27T17:56:20.865475 #94508]  INFO -- : Xextractspans.pl akata_tophat_out/accepted_hits.sam
	I, [2010-12-27T17:56:20.865506 #94508]  INFO -- : Xfilterspans.pl s_4_sequence_Akata.txt.novo.NM.fasta akata_tophat_out/accepted_hits.sam.spans

___

#Configuration

##PARSES Configuration
A configuration file is produced for PARSES in `$HOME/.PARSES`. It contains the paths to the human genome and NT databases as well as the paths to the Bowtie and TopHat indices. It has the following form:

	--- !ruby/object:Settings
	bowtieIndex: /usr/share/hgChrAll
	humanGenomeDatabase: /usr/share
	novoIndex: /usr/share/hgChrAll.ndx
	ntDatabase: /usr/share/nt/nt

##Sequence Configuration
In addition, configuation files are produced for each sequence in `$(pwd)/.sequenceName`. Settings for each program are chosen by default but can be changed via the sequence file which has the following form:

	--- !ruby/object:Sequence
	abyssPath: s_4_sequence_Akata.txt.novo.NM.fasta.nospans.blast.megan.rma
	abyssPathGlob: reads-*.fasta
	blast1Path: s_4_sequence_Akata.txt.novo.NM.fasta.nospans
	blast2Path: s_4_sequence_Akata.txt.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa
	blastOutputFormat: 0
	blastPathGlob: reads-*.fasta.*.kmer.contigs.fa.kmerOptimized.fa
	dataType: illumina1.3
	eValue1: 1.0e-06
	eValue2: 100
	expansionNumber: 10
	filePath: s_4_sequence_Akata.txt
	imageFileType: jpg
	maxKmerLength: 38
	maxMatches: 0
	megan1Path: s_4_sequence_Akata.txt.novo.NM.fasta.nospans.blast
	megan2Path: s_4_sequence_Akata.txt.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast
	minKmerLength: 7
	minScoreByLength: 0
	minSupport: 5
	novoalignPath: s_4_sequence_Akata.txt.novo
	pipeEndPath: s_4_sequence_Akata.txt.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast.megan.rma
	readLength: 38
	removeNonMappedPath: s_4_sequence_Akata.txt.novo.NM.fasta
	topPercent: 10.0
	useCogs: "false"
	useGos: "false"
	winScore: 0.0

It is recommended the path variables not be adjusted. Everything except expansionNumber, which is the number of times to execute the expansion command when generating a picture from MEGAN is straight-forward.

##System Configuration
In addition, the amount of RAM, number of CPU cores, CPU architecture, operating system, default shell, and existence of locate database is automatically computed each execution but not stored.
