require 'logger'
require 'yaml'
require 'net/http'
require 'net/ftp'
require 'uri'
require 'rake/clean'
require 'csv'

#### DISCLAIMER
# This product may not be used for any sort of financial gain. Licenses for both MEGAN and Novoalign are strictly for non-profit research at non-profit institutions and academic usage.

PROG_NAME = 'PARSES'
VER = '0.43'
PROG_DIR = File.dirname(__FILE__)
MEGAN_EXPANSION = 'expand direction=vertical; update;'

# Execute a bash command, time it, and log it.
def safeExec(command, log, sequence, errorMessage)
	results=`{ /usr/bin/time #{command} } 2> .time`
	if $?.exitstatus != 0
		puts "#{errorMessage} (status = #{$?.exitstatus})"
		log.info("#{errorMessage} (status = #{$?.exitstatus})")
	else
		timeFile = File.open(".time")
		File.open(".time").first =~ /(\d+).\d+ real/
		timeFile.close
		sequence.executionTime = sequence.executionTime + $1.to_i
		log.info("#{command} execution time: #{Time.at($1.to_i).gmtime.strftime('%R:%S')}")
	end
	return $?.exitstatus, results
end

# Determine if software is installed
def command?(name)
	`which #{name}`
	$?.success?
end

# Determine if installing an application
installMode = false
for arg in ARGV do
	if arg =~ /install|index|-T|clobber|clean/i
		installMode = true
		break
	end
end

# Columnar representation of a Novoalign file
class NOVOFILE
	READ_HEADER=0
	READ_TYPE=1
	SEQUENCE=2
	QUALITY=3
	ALIGNMENT_TYPE=4
	ALIGNMENT_SCORE=5
	ALIGNMENT_QUALITY=6
	REVERSE_COMPLEMENT_ALIGNMENT_SCORE=7
	ALIGNMENT_HEADER=8
	ALIGNMENT_OFFSET=9
	ALIGNMENT_STRANDEDNESS=10
	PAIRED_END1=11
	PAIRED_END2=12
	PAIRED_END3=13
	MISMATCHES_INSERTIONS_DELETEIONS=14
	HUMAN_TAXONOMY_NUMBER=9606
end

# Links to each program which may be used should the tool not be automatically downloaded properly
class ProgramRepository
	attr_accessor :megan, :tophat, :abyss, :novoalign, :blast, :bowtie, :samtools, :hg
end

# Create MEGAN parsable CSV file from a Novoalign FASTA output file.
def novoToMeganCsv(novoFile)
	output = File.open( "#{novoFile}.csv", 'wb' )
	CSV.foreach(novoFile) do |row|
		rowSplit = row.to_s.split("\t")
		output.write("#{rowSplit[NOVOFILE::READ_HEADER]},#{NOVOFILE::HUMAN_TAXONOMY_NUMBER},#{rowSplit[NOVOFILE::ALIGNMENT_SCORE]}\n") unless rowSplit[0][0].chr == '#'
	end
	output.close
end

# Follow redirecting links. Used to obtain latest versions of files.
def fetch(uri_str, limit = 10)
	raise ArgumentError, 'HTTP redirect too deep' if limit == 0
	response = Net::HTTP.get_response(URI.parse(uri_str))
	case response
	when Net::HTTPSuccess     then response
	when Net::HTTPRedirection then fetch(response['location'], limit - 1)
	else response.error!
	end
end

# Program settings such as paths to various databases and indices
class Settings
	attr_accessor :humanGenomeDatabase, :ntDatabase, :giTaxIdNuclDatabase, :bowtieIndex, :novoIndex
end

# Load program configuration file
if File::exists?(File.expand_path("~/.#{PROG_NAME}")) #Program has been executed before
	progSettingsFile = File.open(File.expand_path("~/.#{PROG_NAME}"))
	progSettings = YAML.load(progSettingsFile)
	progSettingsFile.close
else
	progSettings = Settings.new
end



### GATHER SEQUENCE RUN INFORMATION FROM PREVIOUS RUNS
seqName = "#{ENV['seq']}".chomp
abort 'Must specify a sequence to analyze via the seq= option' if seqName.empty? and !installMode
abort 'Perl must be installed.' if !command? 'perl' and !installMode
abort 'gcc must be installed' if !command? 'gcc'
seqFileName = "#{ENV['file']}".chomp
dataType = "#{ENV['type']}".chomp
useRepo = "#{ENV['repo']}".chomp
truncate = "#{ENV['truncate']}".chomp
forceFile = "#{ENV['forcefile']}".chomp
indexName = "#{ENV['indexName']}".chomp
indexGlob = "#{ENV['indexGlob']}".chomp

# Specify properties about this individual data set/run
class Sequence
	attr_accessor :readLength, :dataType, :executionTime, :filePath, :numOfReads, :novoalignPath, :removeNonMappedPath, :blast1Path, :megan1Path, :abyssPath, :abyssPathGlob, :minKmerLength, :maxKmerLength, :blastPathGlob, :blast2Path, :eValue1, :eValue2, :blastOutputFormat, :megan2Path, :expansionNumber, :maxMatches, :minScoreByLength, :topPercent, :winScore, :minSupport, :useCogs, :useGos, :imageFileType, :pipeEndPath
	def initialize(filePath, dataType)
		if File.exists?(filePath)
			f = File.new(filePath)
			line = f.gets
			line = f.gets while line[0] =~ /\W/
			@readLength=line.length
			f.close
			@numOfReads="#{`egrep -c '^[ACTGN]' "#{filePath}"`}".chomp.to_i
		else
			filePath=''
		end
		@dataType=dataType
		@executionTime = 0.0
#PATHS
		@filePath=filePath
		@novoalignPath="#{filePath}.novo"
		@removeNonMappedPath="#{novoalignPath}.NM.fasta"
		@blast1Path="#{removeNonMappedPath}.nospans"
		@megan1Path="#{blast1Path}.blast"
		@abyssPath="#{megan1Path}.megan.rma"
		@abyssPathGlob="reads-*.fasta"
		@blastPathGlob="#{abyssPathGlob}.*.kmer.contigs.fa.kmerOptimized.fa"
		@blast2Path="#{abyssPath}.kmerOptimized.fa"
		@megan2Path="#{blast2Path}.blast"
		@pipeEndPath="#{megan2Path}.megan.rma"
#ABYSS
		@minKmerLength=7
		@maxKmerLength=readLength
#BLAST
		@eValue1=0.001
		@eValue2=100
		@blastOutputFormat=6
#MEGAN
		@expansionNumber=10
		@minScoreByLength=25
		@topPercent=10.0
		@winScore=0.0
		@minSupport=5
		@useCogs='false'
		@useGos='false'
		@imageFileType='PDF'
	end
end

# Novoalign supports PRB, PRBnSEQ, QSEQ Illumina from Bustard, ABI Solid color space FASTA, ABI Solid color space with qual file, and color space FASTQ.
# Tophat supports colorspace FASTA and space-delimited interger files
# Data types supported
class DataType
	attr_accessor :novoalign, :tophat, :abyss
	def initialize(dataType)
		case dataType
		when 'fasta'
			@novoalign='-F FA'
			@tophat=''
			@abyss = ''
		when 'sanger'
			@novoalign='-F STDFQ'
			@tophat=''
			@abyss='--standard-quality'
		when 'solexa'
			@novoalign='-F SLXFQ'
			@tophat='--solexa-quals'
			@abyss='--illumina-quality'
		when 'illumina1.3'
			@novoalign='-F ILMFQ'
			@tophat='--solexa1.3-quals'
			@abyss='--illumina-quality'
		when 'illumina1.5'
			@novoalign='-F ILMFQ'
			@tophat='--solexa1.3-quals'
			@abyss='--illumina-quality'
		end
	end
end
# Load information concerning this data set/run
if (!installMode)
	if File::exists?(".#{seqName}") #Sequence has been run before
		seqFile = File.open(".#{seqName}")
		sequence = YAML.load(sequenceSettingsFile)
		#Was a new file or file type entered from commandline?
		if seqFileName.empty?
			seqFileName = sequence.filePath
		else
			sequence.filePath = seqFileName
		end
		if dataType.empty?
			dataType = sequence.dataType
		else
			sequence.dataType = dataType
		end
	elsif seqFileName.empty? or dataType.empty? #Sequence has not been run and a file or data type is not specified
		abort "Must specify a file of data to analyze for the first execution of the process as well as it's data type via the file= and type= options" if !installMode
	else #Sequence has not been run, but file and file type are specified
		seqFile = File.open(".#{seqName}")
		sequence = Sequence.new(seqFileName, dataType)
	end
	seqFile.close
	setDataTypes = DataType.new(sequence.dataType)

	log = Logger.new("#{seqName}.log")
	log.info("Begin run for seq=#{seqName} file=#{seqFileName} type=#{dataType} task=#{ARGV[-1]}")
else
	sequence = Sequence.new('~', 'solexa') #garbage data just to get things installed
end

# Clean and clobber
if installMode
	seqNameEmpty = seqName.empty?
	seqName = FileList[".[a-zA-Z0-9]*"] if seqNameEmpty #CLEAN or CLOBBER all files if no sequence is specified
	seqName.each{ | sn |
		seqFile = File.open(".#{sn}")
		seqFileName = YAML.load(seqFile).filePath
		seqFileName = YAML.load(seqFile).filePath if seqNameEmpty and 
		seqFile.close
		CLEAN.include("#{sn}_tophat_out")
		CLEAN.include("#{seqFileName}.novo.counts")
		CLEAN.include("#{seqFileName}.novo.NM.fasta.nospans.details")
		CLEAN.include("#{seqFileName}.novo.NM.fasta.nospans.mergedBlast")
		CLEAN.include("#{seqFileName}.novo.NM.fasta.nospans.[a-z][a-z]")
		CLEAN.include("#{seqFileName}.novo.NM.fasta.nospans.[a-z][a-z].blast")
		CLEAN.include("reads-*.fasta.*.kmer.contigs.fa")
		CLEAN.include("reads-*.fasta.*.kmer.contigs.coverage")
		CLEAN.include(".time")
		CLOBBER.include("#{seqFileName}.*")
		CLOBBER.include("reads-*")
		CLOBBER.exclude("#{seqFileName}.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast.megan.rma")
		CLOBBER.exclude("#{seqFileName}.novo.NM.fasta.nospans.blast.megan.rma.kmerOptimized.fa.blast.pdf")
		CLOBBER.include(".#{sn}")
		CLOBBER.include("#{sn}.log")
	}
end


### GATHER SYSTEM INFORMATION

# Determine which program to use to find files and find the first hit
locate = (command? 'locate') && (!ENV['LOCATE_PATH'].to_s.empty? or File.exists? '/var/lib/mlocate/mlocate.db')
def findFile(filename, locate)
	path = ''
	path = `locate #{filename} | head -1`.chomp if locate
	path = `find / -name #{filename} | head -1`.chomp if path.to_s.empty?
	return path.to_s
end

# Find a novo index of the specificied criteria, make one if not found
def buildNovoIndex(name, pathGlob)
	novoIndex=findFile("#{name}*.ndx", locate)
	if novoIndex.to_s.empty?
		novoIndex="#{File.dirname(pathGlob)}/#{name}.ndx"
		`novoindex "#{novoIndex}" "#{pathGlob}"`
	end
	return novoIndex.to_s.chomp
end

# Find a bowtie index of the specificied criteria, make one if not found
def buildBowtieIndex(name, pathGlob)
	bowtieIndex=findFile("#{name}*.ebwt", locate)
	bowtieIndex= $1 if bowtieIndex =~ /(.*#{name}.*)\.\d+\.ebwt/i
	resourceFiles=''
	if bowtieIndex.to_s.empty?
		bowtieIndex="#{File.dirname(pathGlob)}/#{name}"
		FileList[pathGlob].each do |filename|
			resourceFiles << filename + ','
		end
		resourceFiles.chomp!(',')
		sh %{
			cd "#{File.dirname(pathGlob)}";
			bowtie-build "#{resourceFiles}" "#{name}";
		}
	end
	return bowtieIndex.to_s.chomp
end

# shell=`ps -p $$ | tail -1 | awk '{print $NF}'` #supposedly more accurate method to return shell, but is returning sh instead of bash
shell = File.basename(ENV['SHELL']).chomp

# Represents operating system of current environment
class OS
	OSX=0
	LINUX=1
	BSD=2
	SOLARIS=3
	WINDOWS=4
end

# Determine OS
osName = `uname -s`.chomp!
case osName
when 'Darwin'
	os = OS::OSX
when 'Linux'
	os = OS::LINUX
when /BSD/
	os = OS::BSD
when 'SunOS'
	os = OS::SOLARIS
when /CYGWIN/
	os = OS::WINDOWS
end

arch = `uname -m` =~ /64/ ? 64 : 32 #Determine architecture
ncpu = -1
memInGigs = -1
progRepo = ProgramRepository.new

#Download program links repositories
case os
when OS::OSX
	ncpu = `sysctl -n hw.ncpu`.chomp.to_i # determine number of cpus for MacOS
	"#{`/usr/sbin/system_profiler SPHardwareDataType | grep Memory`}" =~ /Memory:\s+(\d+)\s*GB/; # determine amount of memory in MacOS
	memInGigs = $1.chomp.to_i
	if useRepo
		Net::HTTP.start('cloud.github.com', 80) { |http| # download repository for MacOS
			progRepo = YAML.load(http.get('/downloads/Lythimus/PARSES/.programRepositoryMac.txt').body)
		}
	end
when OS::LINUX
	`cat /proc/cpuinfo | grep -G processor.*:.* | tail -n 1` =~ /processor.*:.*(\d+)/
	ncpu = $1.chomp.to_i + 1 # determine number of cpus for Linux
	memInGigs = %x[echo `cat /proc/meminfo | grep MemTotal` | sed  "s/[^0-9]//g"].to_i/2**20 # determine amount of memory in Linux
	if useRepo
		Net::HTTP.start('cloud.github.com', 80) { |http| # download repository for Linux
			progRepo = YAML.load(http.get('/downloads/Lythimus/PARSES/.programRepositoryLinux.txt').body)
		}
	end
end

## consider adding time of processing to all logs
log.info("Executing #{PROG_NAME + ' v' + VER} in a #{osName} environment with #{memInGigs}GB of memory and #{ncpu} processing cores with a #{arch}-bit architecture.") if !installMode

######### PIPELINE
desc 'Novoalign - Align reads in to base genome'
task :alignSequence => [:novoalignInstall, :hgInstall, :novoIndex]
file sequence.novoalignPath => sequence.filePath do
	puts 'Sequence Alignment'
	if indexName.empty? or indexGlob.empty?
		novoIndex = progSettings.novoIndex
	else
		novoIndex = buildNovoIndex(indexName, indexGlob)
	end
	seqFileName = sequence.filePath if forceFile != 'true'
	exitStatus = safeExec("novoalign -d \"#{novoIndex}\" #{setDataTypes.novoalign} -f \"#{seqFileName}\" > \"#{sequence.novoalignPath}\";", log, sequence,
			'Novoalign sequence alignment not performed')
	if exitStatus == 0
		readsLeft=`egrep -c '@' "#{sequence.novoalignPath}"`
		log.info("#{readsLeft} human reads discovered which is #{sequence.numOfReads/readsLeft*100}% of total.")
	end
end

desc 'Xnovotonm - Harvest non-base organism reads'
task :removeHuman => :alignSequence
file sequence.removeNonMappedPath => sequence.novoalignPath do
	puts 'RemoveHuman'
	seqFileName = sequence.novoalignPath if forceFile != 'true'
	exitStatus = safeExec("\"#{PROG_DIR}/Xnovotonm.pl\" \"#{seqFileName}\";", log, sequence,
			'Removal of Human mapped reads from novoalign results not performed')
	if exitStatus == 0
		readsLeft=`egrep -c '^[ACTGN]' "#{sequence.removeNonMappedPath}"`
		log.info("Reduced to #{readsLeft} non-mapped reads which is #{sequence.numOfReads/readsLeft*100}% of total.")
	end
end

desc 'Tophat - Align spanning reads in order to remove base organism reads'
task :removeSpans => [:bowtieIndex, :tophatInstall, :removeHuman]
file sequence.blast1Path => sequence.removeNonMappedPath do
	puts 'RemoveSpans'
	## DIDN'T IMPLEMENT forceFile BECAUSE TOO COMPLICATED
	if indexName.empty? or indexGlob.empty?
		bowtieIndex = progSettings.bowtieIndex
	else
		bowtieIndex = buildBowtieIndex(indexName, indexGlob)
	end
	safeExec("tophat -p #{ncpu} #{setDataTypes.tophat} --output-dir \"#{ENV['seq']}_tophat_out\" \"#{bowtieIndex}\" \"#{sequence.filePath}\";", log, sequence,
			'TopHat sequence alignment not performed')
	safeExec("samtools view -h -o \"#{ENV['seq']}_tophat_out/accepted_hits.sam\" \"#{ENV['seq']}_tophat_out/accepted_hits.bam\";", log, sequence,
			'Samtools conversion not performed')
	safeExec("\"#{PROG_DIR}/Xextractspans.pl\" \"#{ENV['seq']}_tophat_out/accepted_hits.sam\";", log, sequence,
			'Spanning region extraction not performed')
	exitStatus = safeExec("\"#{PROG_DIR}/Xfilterspans.pl\" \"#{sequence.removeNonMappedPath}\" \"#{ENV['seq']}_tophat_out/accepted_hits.sam.spans\";", log, sequence,
			'Filter of spanning regions not performed')
	if exitStatus == 0
		readsLeft=`egrep -cv '>' "#{sequence.blast1Path}"`
		log.info("Reduced to #{readsLeft} reads which is #{sequence.numOfReads/readsLeft*100}% of total.")
	end
end

desc 'BLAST - Associate reads with organisms.'
task :localAlignReads => [:blastInstall, :ntInstall, :removeSpans]
file sequence.megan1Path => sequence.blast1Path do
	puts 'localAlignReads'
	dust = "-dust no" if sequence.readLength < 50 # Dust filtering should be disabled for short reads
	seqFileName = sequence.blast1Path if forceFile != 'true'
	pieces = ([(memInGigs/30), ncpu].min)
	pieces = 1 if pieces == 0
	pieceSize = `wc -l "#{seqFileName}"`.chomp.to_i / pieces
	pieceSize = pieceSize + (pieceSize % 2) # This needs to be changed if paired-end read suppor tis ever implemented
	`split -l #{pieceSize} "#{seqFileName}" "#{seqFileName}."`
	fileCount = pieces
	@blastCommands = []
	FileList["#{seqFileName}.[a-zA-Z0-9][a-zA-Z0-9]"].each { | blastPiece |
		@blastCommands << "blastn -db \'#{progSettings.ntDatabase}\' -soft_masking true #{dust} -num_threads #{ncpu} -evalue #{sequence.eValue1} -outfmt #{sequence.blastOutputFormat} -query \'#{blastPiece}\' -out \'#{blastPiece}.blast\' &"
		fileCount = fileCount - 1
		if fileCount == 0
			safeExec("\"#{PROG_DIR}/parallelBlast.sh\" #{@blastCommands.join(' ')}", log, sequence, #This is not a space in the join, rather looks like a space. It is used by the BASH script as a delimiter.
					'BLAST of reads not performed')
			@blastCommands = []
			fileCount = pieces
		end
	}
	`cat #{seqFileName}.[a-z][a-z].blast > #{seqFileName}.mergedBlast`
	safeExec("\"#{PROG_DIR}/addTaxon.pl\" \"#{progSettings.giTaxIdNuclDatabase.to_s}\" \"#{seqFileName}.mergedBlast\" \"#{seqFileName}\";", log, sequence,
			'Adding taxon to end of file not performed')
end

desc 'MEGAN - Separate reads into taxonomies.'
task :metaGenomeAnalyzeReads => [ :meganInstall, :localAlignReads]
file sequence.abyssPath => sequence.megan1Path do
	puts 'metaGenomeAnalyzeReads'
	## DIDN'T IMPLEMENT forceFile BECAUSE TOO COMPLICATED
	`MEGAN +g -V -E -x 'CRASHPROGRAMPLS'` =~ /MEGAN.*version\s*(\d+\.?\d*)/
	if ($1.to_f < 4.0)
		safeExec("MEGAN +g -E -x \"import blastfile='#{sequence.megan1Path}' readfile='#{sequence.blast1Path}' meganfile='#{sequence.abyssPath}' minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format='#{sequence.imageFileType}' file='#{sequence.megan1Path + '.' + sequence.imageFileType.downcase}' REPLACE=true; quit;\";", log, sequence,
			'MEGAN processing of BLASTed reads not performed')
		exitStatus, results = safeExec("MEGAN -f \"#{sequence.abyssPath}\" -x \"#{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all;\";", log, sequence,
			'Opening MEGAN file not performed')
	else
		safeExec("MEGAN +g -E -x \"import blastfile='#{sequence.megan1Path}' readfile='#{sequence.blast1Path}' meganfile='#{sequence.abyssPath}' minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; set context=seedviewer; #{MEGAN_EXPANSION*sequence.expansionNumber} select nodes=all; uncollapse subtrees; update; exportimage format='#{sequence.imageFileType}' file='#{sequence.megan1Path + '.' + sequence.imageFileType.downcase}' REPLACE=true; quit;\";", log, sequence,
			'MEGAN processing of BLASTed reads not performed')
		exitStatus, results = safeExec("MEGAN -f \"#{sequence.abyssPath}\" -x \"set context=seedviewer; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all;\";", log, sequence,
			'Opening MEGAN file not performed')
	end
	if exitStatus == 0
		totalReads = $_.chomp.to_i if (results =~ /Total reads:\s*(\d+)/)
		assignedReads = $_.chomp.to_i if (results =~ /Assigned reads:\s*(\d+)/)
		unassignedReads = $_.chomp.to_i if (results =~ /Unassigned reads:\s*(\d+)/)
		noHits = $_.chomp.to_i if (results =~ /Reads with no hits:\s*(\d+)/)
		log.info("Potential Exogenous Reads: #{totalReads}. Assigned Reads: #{assignedReads} which is #{assignedReads/sequence.numOfReads}% of entire data set. Unassigned Reads: #{unassignedReads} which means LCA hides #{unassignedReads/totalReads}% of the data. No hits: #{noHits} which is #{noHits/sequence.numOfReads} of the entire data set.")
	end
end

desc 'ABySS - Assemble reads associated with clusters of taxonomies.'
task :denovoAssembleCluster => [:abyssInstall, :parallelIteratorInstall, :metaGenomeAnalyzeReads]
file sequence.blast2Path => FileList["#{sequence.abyssPathGlob}"] do
	puts 'denovoAssemblyCluster'
	seqFileName = sequence.abyssPathGlob if forceFile != 'true'
	exitStatus=-1
	FileList["#{seqFileName}"].each { | abyssFiles |
		exitStatus = safeExec("\"#{PROG_DIR}/abyssKmerOptimizer.pl\" #{abyssFiles} #{sequence.minKmerLength} #{sequence.maxKmerLength} #{setDataTypes.abyss};", log, sequence,
			'ABySS not performed')
	}
	`cat #{sequence.blastPathGlob} > #{sequence.blast2Path}` if forceFile != 'true'
	if exitStatus == 0
		coverageThreshold = $_ if (results =~ /(Using a coverage threshold of \d+)/)
		medianKmerCoverage = $_ if (results =~ /(The median k-mer coverage is \d+)/)
		assembly = $_ if (results =~ /(Assembled \d+ k-mer in \d+ contigs)/)
		FileList["#{sequence.blastPathGlob}"].first =~ /#{sequence.abyssPathGlob}.(\d+)\.kmer\.contigs\.fa\.kmerOptimized\.fa/ if forceFile != 'true'
		kmer = $_
		log.info("#{coverageThreshold}. #{medianKmerCoverage}. #{assembly}. The optimum kmer is #{kmer}.")
	end
end

desc 'BLAST - Associate contigs with organisms.'
task :localAlignContigs => [:blastInstall, :ntInstall, :denovoAssembleCluster]
file sequence.megan2Path => sequence.blast2Path do
	puts 'localAlignContigs'
	seqFileName = sequence.blast2Path if forceFile != 'true'
	safeExec("blastn -db \"#{progSettings.ntDatabase}\" -soft_masking true -num_threads #{ncpu} -evalue #{sequence.eValue2} -outfmt #{sequence.blastOutputFormat} -query \"#{seqFileName}\" -out \"#{sequence.megan2Path}.noTax\";", log, sequence,
			'BLAST of contigs not performed')
	safeExec("\"#{PROG_DIR}/addTaxon.pl\" \"#{progSettings.giTaxIdNuclDatabase.to_s}\" \"#{sequence.megan2Path}.noTax\" \"#{seqFileName}\";", log, sequence,
			'Adding taxon to end of BLAST contigs file not performed')
end

desc 'MEGAN - Separate contigs into taxonomies.'
task :metaGenomeAnalyzeContigs => [ :meganInstall, :localAlignContigs]
file sequence.pipeEndPath => sequence.megan2Path do
	## DIDN'T IMPLEMENT forceFile BECAUSE TOO COMPLICATED
	puts 'metaGenomeAnalyzeContigs'
	`MEGAN +g -V -E -x 'CRASHPROGRAMPLS'` =~ /MEGAN.*version\s*(\d+\.?\d*)/
	if ($1.to_f < 4.0)
		safeExec("MEGAN +g -E -x \"import blastfile='#{sequence.megan2Path}' readfile='#{sequence.blast2Path}' meganfile='#{sequence.pipeEndPath}' minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format='#{sequence.imageFileType}' file='#{sequence.megan2Path + '.' + sequence.imageFileType.downcase}' REPLACE=true; quit;\";", log, sequence,
				'MEGAN processing of BLASTed contigs not performed')
		exitStatus, results = safeExec("MEGAN -f \"#{sequence.pipeEndPath}.rma\";", log, sequence,
				'Opening MEGAN file not performed')
	else
		safeExec("MEGAN +g -E -x \"import blastfile='#{sequence.megan2Path}' readfile='#{sequence.blast2Path}' meganfile='#{sequence.pipeEndPath}' minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; update; set context=seedviewer; #{MEGAN_EXPANSION*sequence.expansionNumber} select nodes=all; uncollapse subtrees; update; exportimage format='#{sequence.imageFileType}' file='#{sequence.megan2Path + '.' + sequence.imageFileType.downcase}' REPLACE=true; quit;\";", log, sequence,
				'MEGAN processing of BLASTed contigs not performed')
		exitStatus, results = safeExec("MEGAN -f \"#{sequence.pipeEndPath}.rma\";", log, sequence,
				'Opening MEGAN file not performed')
	end
	if exitStatus == 0
		totalReads = $_.chomp.to_i if (results =~ /Total reads:\s*(\d+)/)
		assignedReads = $_.chomp.to_i if (results =~ /Assigned reads:\s*(\d+)/)
		unassignedReads = $_.chomp.to_i if (results =~ /Unassigned reads:\s*(\d+)/)
		noHits = $_.chomp.to_i if (results =~ /Reads with no hits:\s*(\d+)/)
		log.info("Potential Exogenous Reads: #{totalReads}. Assigned Reads: #{assignedReads}. Unassigned Reads: #{unassignedReads}. No hits: #{noHits}.")
	end
end



######### INSTALLATIONS

desc 'Install latest version of human genome database.'
task :hgInstall do
	if progSettings.humanGenomeDatabase.to_s.chomp.empty?
		progSettings.humanGenomeDatabase=File.dirname(findFile('chr*.fa', locate))
		if progSettings.humanGenomeDatabase.to_s.chomp == '.'
			hg = ''
			ftp = Net::FTP.new('hgdownload.cse.ucsc.edu')
			ftp.login()
			ftp.passive=true
			ftp.chdir('apache/htdocs/goldenPath')
			ftp.list("hg[0-9][0-9]").last =~ /.*(hg\d+).*/
			hg = $1
			ftp.chdir(hg + '/bigZips')
			ftp.getbinaryfile('chromFa.tar.gz')
			ftp.close
			sh %{
				mkdir /usr/share/#{hg};
				tar -xzf chromFa.tar.gz -C /usr/share/#{hg};
				rm chromFa.tar.gz;
			}
			progSettings.humanGenomeDatabase = '/usr/share/' + hg
		end
	end
end

desc 'Create an index for novoalign of the human genome database.'
task :novoIndex => [:hgInstall, :novoalignInstall] do
	if progSettings.novoIndex.to_s.empty?
		progSettings.novoIndex=buildNovoIndex('hgChrAll', "#{progSettings.humanGenomeDatabase.to_s}/chr[0-9XY]*.fa")
	end
end

desc 'Create an index for bowtie/tophat of the human genome database.'
task :bowtieIndex => [:hgInstall, :bowtieInstall] do
	if progSettings.bowtieIndex.to_s.empty?
		progSettings.bowtieIndex=buildBowtieIndex('hgChrAll', "#{progSettings.humanGenomeDatabase.to_s}/chr[0-9XY]*.fa")
	end
end

desc 'Install latest version of Novoalign.'
task :novoalignInstall do
	if !command? 'novoalign'
		novoalign = ''
		Net::HTTP.start('www.novocraft.com', 80) { |http|
			http.get('/main/releases.php', 'Referer' => 'http://www.novocraft.com/').body =~ /(V\d+\.\d+\.\d+)/
			novoalign = $1
			File.open('novocraft.tar.gz', 'w'){ |file|
				if useRepo == true
					link = progRepo.novoalign
				else
					case os
					when OS::OSX then
						link = '/downloads/' + novoalign + '/novocraft' + novoalign + '.MacOSX.tar.gz'
					when OS::LINUX then
						link = '/downloads/' + novoalign + '/novocraft' + novoalign + '.gcc.tar.gz'
					end
				end
				file.write(http.get(link, 'Referer' => 'http://www.novocraft.com/').body)
			}
		}
		sh %{
			tar -xzf novocraft.tar.gz -C /usr/bin;
			ln -s /usr/bin/novocraft/isnovoindex /usr/bin;
			ln -s /usr/bin/novocraft/novo2maq /usr/bin;
			ln -s /usr/bin/novocraft/novo2paf /usr/bin;
			ln -s /usr/bin/novocraft/novoalign /usr/bin;
			ln -s /usr/bin/novocraft/novobarcode /usr/bin;
			ln -s /usr/bin/novocraft/novoindex /usr/bin;
			ln -s /usr/bin/novocraft/novoutil /usr/bin;
		}
	end
end

desc 'Install latest version of Bowtie.'
task :bowtieInstall do
	if !command? 'bowtie'
		bowtie = ''
		Net::HTTP.start('bowtie-bio.sourceforge.net', 80) { |http|
			http.get('/index.shtml').body =~ /https:\/\/sourceforge\.net\/projects\/bowtie-bio\/files\/bowtie\/(\d*\.\d*\.\d*)/
		}
		Net::HTTP.start('softlayer.dl.sourceforge.net', 80) { |http|
			bowtie = $1
			File.open('bowtie-' + bowtie + '-src.zip', 'w'){ |file|
				if useRepo == true
					link = progRepo.bowtie
				else
					link = '/project/bowtie-bio/bowtie/' + bowtie + '/bowtie-' + bowtie + '-src.zip'
				end
				file.write(http.get(link).body)
			}
		}
		sh %{
			unzip -d /usr/bin bowtie-#{bowtie}-src.zip;
			cd /usr/bin/bowtie-#{bowtie};
			make;
			ln -s /usr/bin/bowtie-#{bowtie}/bowtie /usr/bin;
			ln -s /usr/bin/bowtie-#{bowtie}/bowtie-build /usr/bin;
			ln -s /usr/bin/bowtie-#{bowtie}/bowtie-inspect /usr/bin;
		}
	end
end

desc 'Install latest version of Samtools.'
task :samtoolsInstall do
	if !command? 'samtools'
		samtools = ''
		File.open('samtools.tar.bz2', 'w'){ |file|
			if useRepo == true
				link = progRepo.samtools
			else
				link = 'http://sourceforge.net/projects/samtools/files/latest'
			end
			file.write(fetch(link).body)
		}
		`tar -jxf samtools.tar.bz2 -C /usr/bin`
		`rm samtools.tar.bz2`
		samtools = `ls /usr/bin | grep samtools | head -1`.chomp!
		sh %{
			cd /usr/bin/#{samtools};
			make;
			cp libbam.a /usr/lib;
			mkdir /usr/include/bam;
			cp *.h /usr/include/bam;
			ln -s /usr/bin/#{samtools}/samtools /usr/bin;
		}
	end
end

desc 'Install latest version of Tophat.'
task :tophatInstall => [:samtoolsInstall, :bowtieInstall] do
	if !command? 'tophat'
		tophat = ''
		Net::HTTP.start('tophat.cbcb.umd.edu', 80) { |http|
			http.get('/index.html').body =~ /\.\/downloads\/(tophat-\d+\.\d+\.\d+\.tar\.gz)/
			tophat = $1
			File.open(tophat, 'w'){ |file|
				if useRepo == true
					link = progRepo.tophat
				else
					link = '/downloads/' + tophat
				end
				file.write(http.get(link).body)
			}
		}
		sh %{
			tar -xzf #{tophat} -C /usr/bin;
			rm #{tophat};
			cd /usr/bin/#{tophat.chomp!('.tar.gz')};
			./configure;
			make;
			make install;
		}
	end
end

desc 'Install latest version of ABySS with Google Sparsehash.'
task :abyssInstall do
	if !command? 'ABYSS'
		Net::HTTP.start('code.google.com', 80) { |http|
			http.get('/p/google-sparsehash/').body =~ /(http:\/\/google-sparsehash\.googlecode\.com)(\/files\/sparsehash-\d+\.\d+\.tar\.gz)/
		}
		base = File.basename($1)
		gsh = $2
		Net::HTTP.start(base, 80) { |http|
			File.open(File.basename(gsh), 'w'){ |file|
				if useRepo == true
					link = progRepo.abyss
				else
					link = gsh
				end
				file.write(http.get(link).body)
			}
		}
		sh %{
			tar -xzf #{File.basename(gsh)} -C /usr/bin;
			rm #{gsh};
			cd /usr/bin/#{File.basename(gsh, '.tar.gz')};
			./configure;
			make;
			make install;
		}

		abyss = ''
		Net::HTTP.start('www.bcgsc.ca', 80) { |http|
			http.get('/downloads/abyss/?C=N;O=D').body =~ /(abyss-\d+\.\d+\.\d+\.tar\.gz)/
			abyss = $1
			File.open(abyss, 'w'){ |file|
				file.write(http.get('/downloads/abyss/' + abyss).body)
			}
		}
		sh %{
			tar -xzf #{abyss} -C /usr/bin;
			rm #{abyss};
			cd /usr/bin/#{abyss.chomp('.tar.gz')};
			./configure CPPFLAGS=-I/usr/local/include --prefix=/usr/bin --enable-maxk=96 && make && make install;
		}
	end
end

desc 'Install latest version of BLAST+.'
task :blastInstall do
	blast = ''
	if !command? 'blastn'
		ftp = Net::FTP::new('ftp.ncbi.nlm.nih.gov')
		ftp.login()
		ftp.passive=true
		ftp.chdir('blast/executables/blast+/LATEST')
		ftp.list("ncbi-blast*+-src.tar.gz").last =~ /.*(ncbi-blast.*\+\-src\.tar\.gz).*/
		blast = $1
		ftp.getbinaryfile(blast)
		ftp.close
		sh %{
			tar -xzf #{blast} -C /usr/bin;
			cd /usr/bin/#{blast.chomp!('.tar.gz')}/c++;
			./configure --without-debug --with-mt --with-build-root=ReleaseMT;
			cd ReleaseMT/build;
			make all_r;
			
			#sudo rm /usr/bin/blastdb_aliastool /usr/bin/dustmasker /usr/bin/rpstblastn /usr/bin/blastdbcheck /usr/bin/gene_info_reader /usr/bin/segmasker /usr/bin/blastdbcmd /usr/bin/gumbelparams /usr/bin/seqdb_demo /usr/bin/blast_formatter /usr/bin/srsearch /usr/bin/blastn /usr/bin/makeblastdb /usr/bin/tblastn /usr/bin/blastp /usr/bin/makembindex /usr/bin/tblastx /usr/bin/blastx /usr/bin/project_tree_builder /usr/bin/convert2blastmask /usr/bin/psiblast /usr/bin/windowmasker /usr/bin/datatool /usr/bin/rpsblast
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastdb_aliastool /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/dustmasker /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/rpstblastn /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastdbcheck /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/gene_info_reader /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/segmasker /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastdbcmd /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/gumbelparams /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/seqdb_demo /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blast_formatter /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/srsearch /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastn /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/makeblastdb /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/tblastn /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastp /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/makembindex /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/tblastx /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/blastx /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/project_tree_builder /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/convert2blastmask /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/psiblast /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/windowmasker /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/datatool /usr/bin;
			ln -s /usr/bin/#{blast}/c++/ReleaseMT/bin/rpsblast /usr/bin;
		}
	end
end

desc 'Install latest version of the NT database.'
task :ntInstall => :blastInstall do
	if ENV['BLASTDB'].to_s.empty?
		if progSettings.ntDatabase.to_s.empty?
			progSettings.ntDatabase = findFile('nt.nal', locate).to_s.chomp('.nal')
			if progSettings.ntDatabase.to_s.empty?
				progSettings.ntDatabase='/usr/share/nt/nt'
				sh %{
					mkdir /usr/share/nt;
					cd /usr/share/nt;
					/usr/bin/ncbi-blast*/c++/src/app/blast/update_blastdb.pl nt;
					for i in nt*.tar.gz; do tar -xzf $i; rm $i; done;
					echo "export BLASTDB=#{progSettings.ntDatabase.to_s}" >> ~/.#{shell}rc;
				}
			end
		else
			`echo "export BLASTDB=#{progSettings.ntDatabase.to_s}" >> ~/.#{shell}rc`
			`export BLASTDB=#{progSettings.ntDatabase.to_s}`
		end
	else
		progSettings.ntDatabase = ENV['BLASTDB'].to_s.chomp('.nal') if progSettings.ntDatabase.to_s.empty?
	end
end

desc 'Install latest version of MEGAN.'
task :meganInstall do
	if !command? 'MEGAN'
		megan = ''
		Net::HTTP.start('www-ab.informatik.uni-tuebingen.de', 80) { |http|
			http.get('/data/software/megan/download/welcome.html').body =~ /(V\d+_\d+\/MEGAN_unix_\d+_\d+\.sh)/
			megan = $1
			File.open("#{File.basename(megan)}", 'w'){ |file|
			if useRepo == true
				link = progRepo.megan
			else
				link = '/data/software/megan/download/' + megan
			end
			file.write(http.get(link).body)
			}
		}
		megan = File.basename(megan)
		sh %{
			chmod +x #{megan};
			./#{megan};
			rm #{megan};
		}
		megan = `which MEGAN`.chomp
		`ln -s "#{findFile(MEGAN)}" /usr/bin` if megan.empty?
		megan = `which MEGAN`.chomp
		if (arch == 64) # If CPU architecture is 64-bit, allow for more than 2GB of RAM and force 64-bit Java.
			text = File.new(megan).read.gsub(/"\$prg_dir\/\$progname" "-server" "-Xms\d+." "-Xmx\d+."/, "\"$prg_dir/$progname\" \"-server\" \"-d64\" \"-Xms#{memInGigs}G\" \"-Xmx#{memInGigs}G\"")
			File.open(megan, 'w+'){ |file| file.write(text) }
		end
	end
end

desc 'Install GI to Taxonomy ID database.'
task :giTaxIdNuclInstall do
	progSettings.giTaxIdNuclDatabase=findFile('gi_taxid_nucl.dmp', locate)
	if progSettings.giTaxIdNuclDatabase.to_s.chomp == '.'
		ftp = Net::FTP::new('ftp://ftp.ncbi.nih.gov')
		ftp.login()
		ftp.passive=true
		ftp.chdir('pub/taxonomy')
		ftp.getbinaryfile('gi_taxid_nucl.zip')
		ftp.close
		`unzip -d /usr/share gi_taxid_nucl.zip`
		progSettings.giTaxIdNuclDatabase='/usr/share/gi_taxid_nucl.dmp'
	end
end

desc 'Install latest version of Parallel::Iterator for perl.'
task :parallelIteratorInstall do
	`perl -MParallel::Iterator -e 1`
	if $?.exitstatus != 0
		`perl -MCPAN -e 'install Parallel::Iterator'`
	end
end

desc 'Index the genome of a specified organism with Novo and Bowtie.'
task :otherIndex do
	puts 'OtherIndex'
	if !indexName.empty? and !indexGlob.empty?
		buildNovoIndex(indexName, indexGlob)
		buildBowtieIndex(indexName, indexGlob)
	end
end


task :default do
	Rake::Task[:metaGenomeAnalyzeContigs].invoke
end

desc 'Install latest version of everything.'
task :install do
	#$ replace with regular expression on for each task at some point
	Rake::Task[:novoalignInstall].invoke
	Rake::Task[:bowtieInstall].invoke
	Rake::Task[:hgInstall].invoke
	Rake::Task[:novoIndex].invoke
	Rake::Task[:bowtieIndex].invoke
	Rake::Task[:samtoolsInstall].invoke
	Rake::Task[:tophatInstall].invoke
	Rake::Task[:abyssInstall].invoke
	Rake::Task[:blastInstall].invoke
	Rake::Task[:ntInstall].invoke
	Rake::Task[:meganInstall].invoke
	Rake::Task[:giTaxIdNuclInstall].invoke
	Rake::Task[:parallelIteratorInstall].invoke
end

desc 'Automatically saving any settings changes which may have been made'
task :reserialize do
	# Reserialize object in case any changes have been made
	seqFile = File.open(".#{seqName}", 'w')
	progSettingsFile = File.open(File.expand_path("~/.#{PROG_NAME}"), 'w')
	YAML.dump(sequence, seqFile) if !installMode
	YAML.dump(progSettings, progSettingsFile)
	seqFile.close
	progSettingsFile.close
end

#Allow for pipeline override, chopping off the front of the pipeline
if truncate != 'true'
	task :alignSequence => sequence.novoalignPath if !File.exists?(sequence.novoalignPath)
	task :removeHuman => sequence.removeNonMappedPath if !File.exists?(sequence.removeNonMappedPath)
	task :removeSpans => sequence.blast1Path if !File.exists?(sequence.blast1Path)
	task :localAlignReads => sequence.megan1Path if !File.exists?(sequence.megan1Path)
	task :metaGenomeAnalyzeReads => sequence.abyssPath if !File.exists?(sequence.abyssPath)
	task :denovoAssembleCluster => sequence.blast2Path if !File.exists?(sequence.blast2Path)
	task :localAlignContigs => sequence.megan2Path if !File.exists?(sequence.megan2Path)
	task :metaGenomeAnalyzeContigs => sequence.pipeEndPath if !File.exists?(sequence.pipeEndPath)
end

# invoke reserialize after all tasks have been executed
current_tasks =  Rake.application.top_level_tasks
current_tasks << :reserialize
Rake.application.instance_variable_set(:@top_level_tasks, current_tasks)
