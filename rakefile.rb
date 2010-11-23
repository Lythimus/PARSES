require 'logger'
require 'yaml'
require 'net/http'
require 'net/ftp'
require 'uri'
require 'rake/clean'

#### DISCLAIMER
# This product may not be used for any sort of financial gain. Licenses for both MEGAN and Novoalign are strictly for non-profit research at non-profit institutions and academic usage.

PROG_NAME = "PARSES"
VER = "0.27"
PROG_DIR = File.dirname(__FILE__)
MEGAN_EXPANSION = 'expand direction=vertical; update;'

def command?(name)
  `which #{name}`
  $?.success?
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

class Settings
	attr_accessor :humanGenomeDatabase, :ntDatabase, :bowtieIndex, :novoIndex
end

if File::exists?(File.expand_path("~/.#{PROG_NAME}")) #Program has been executed before
	progSettings = YAML.load_file(File.expand_path("~/.#{PROG_NAME}"))
else
	progSettings = Settings.new
end

### GATHER SEQUENCE RUN INFORMATION FROM PREVIOUS RUNS
seqName = "#{ENV['seq']}".chomp
abort "Must specify a sequence to analyze via the seq= option" if seqName.empty?
abort "Perl must be installed." if !command? "perl"
abort "gcc must be installed" if !command? "gcc"
seqFileName = "#{ENV['file']}".chomp
dataType = "#{ENV['type']}".chomp

class Sequence
	attr_accessor :readLength, :dataType, :filePath, :novoalignPath, :removeNonMappedPath, :blast1Path, :megan1Path, :abyssPath, :abyssPathGlob, :minKmerLength, :maxKmerLength, :blastPathGlob, :eValue1, :eValue2, :blastOutputFormat, :megan2PathGlob, :expansionNumber, :maxMatches, :minScoreByLength, :topPercent, :winScore, :minSupport, :useCogs, :useGos, :imageFileType, :pipeEndGlob
	def initialize(filePath, dataType)
		f = File.new(filePath)
		line = f.gets
		line = f.gets while line[0] =~ /\W/
		@readLength = line
		f.close
		@dataType = dataType
#PATHS
		@filePath = filePath
		@novoalignPath="#{filePath}.novo"
		@removeNonMappedPath="#{novoalignPath}.NM.fasta"
		@blast1Path="#{removeNonMappedPath}.nospans"
		@megan1Path="#{blast1Path}.blast"
		@abyssPath="#{megan1Path}.megan.rma"
		@abyssPathGlob="reads-*.fasta"
		@blastPathGlob="#{abyssPathGlob}.*.kmer.contigs.fa.kmerOptimized.fa"
		@megan2PathGlob="#{blastPathGlob}.blast"
		@pipeEndGlob="#{megan2PathGlob}.megan.rma"
#ABYSS
		@minKmerLength=7
		@maxKmerLength=readLength
#BLAST
		@eValue1=0.001
		@eValue2=100
		@blastOutputFormat=5
#MEGAN
		@expansionNumber=10
		@maxMatches=0
		@minScoreByLength=0
		@topPercent=10.0
		@winScore=0.0
		@minSupport=5
		@useCogs="false"
		@useGos="false"
		@imageFileType="jpg"
	end
end

# Novoalign supports PRB, PRBnSEQ, QSEQ Illumina from Bustard, ABI Solid color space FASTA, ABI Solid color space with qual file, and color space FASTQ.
# Tophat supports colorspace FASTA and space-delimited interger files
class DataType
	attr_accessor :novoalign, :tophat, :abyss
	def initialize(dataType)
		case dataType
		when 'fasta'
			@novoalign='-F FA'
			@tophat=''
			@abyss = ''
		when 'sanger'
			@novoalign="-F STDFQ"
			@tophat=''
			@abyss="--standard-quality"
		when 'solexa'
			@novoalign='-F SLXFQ'
			@tophat='--solexa-quals'
			@abyss="--illumina-quality"
		when 'illumina1.3'
			@novoalign='-F ILMFQ'
			@tophat='--solexa1.3-quals'
			@abyss="--illumina-quality"
		when 'illumina1.5'
			@novoalign='-F ILMFQ'
			@tophat='--solexa1.3-quals'
			@abyss="--illumina-quality"
		end
	end
end


if File::exists?(".#{seqName}") #Sequence has been run before
	sequence = YAML.load_file(".#{seqName}")
	seqFileName = sequence.filePath
	dataType = sequence.dataType
elsif seqFileName.empty? or dataType.empty? #Sequence has not been run and a file or data type is not specified
	abort "Must specify a file of data to analyze for the first execution of the process as well as it's data type via the file= and type= options"
else #Sequence has not been run, but file and file type are specified
	seqFile = File.open(".#{seqName}", "a")
	sequence = Sequence.new(seqFileName, dataType)
end
setDataTypes = DataType.new(sequence.dataType)

log = Logger.new("#{seqName}.log")
log.info("Begin run for seq=#{seqName} file=#{seqFileName} type=#{dataType} task=#{ARGV[0]}")

# Clean and clobber are not functioning at the moment
#CLEAN.include(sequence.filePath + '\S*')
#CLOBBER.include(sequence.filePath + '\S*')

### GATHER SYSTEM INFORMATION

if command? "locate" and (!ENV['LOCATE_PATH'].to_s.empty? or File.exists? '/var/lib/mlocate/mlocate.db')
	find="locate"
else
	find="find / -name"
end
shell=`ps -p $$ | tail -1 | awk '{print $NF}'`

class OS
	OSX=0
	LINUX=1
	BSD=2
	SOLARIS=3
	WINDOWS=4
end

os = `uname -s`.chomp!
case os
when 'Darwin'
	os = OS::OSX
when 'Linux'
	os = OS::LINUX
when /.*BSD.*/
	os = OS::BSD
when 'SunOS'
	os = OS::SOLARIS
when /.*CYGWIN.*/
	os = OS::WINDOWS
end

arch = `uname -m` =~ /.*64.*/ ? 64 : 32
ncpu = -1
memInGigs = -1

case os
when OS::OSX
	ncpu = `sysctl -n hw.ncpu`.chomp.to_i
	"#{`/usr/sbin/system_profiler SPHardwareDataType | grep Memory`}" =~ /Memory:\w+(\d+).*/;
	memInGigs = $1
when OS::LINUX
	`cat /proc/cpuinfo | grep -G processor.*:.* | tail -n 1` =~ /processor.*:.*(\d+)/
	ncpu = $1.chomp.to_i + 1
	memInGigs = %x[echo `cat /proc/meminfo | grep MemTotal` | sed  "s/[^0-9]//g"].to_i/2**20
end

## consider adding time of processing to all logs
log.info("Executing #{PROG_NAME + ' v' + VER} in a #{os} environment with #{memInGigs}GB of memory and #{ncpu} cores with a #{arch}-bit architecture.")

######### PIPELINE
desc "Novoalign - Align reads in to base genome"
task :alignSequence => [:novoalignInstall, :hgInstall, :novoIndex, sequence.novoalignPath]
file sequence.novoalignPath => sequence.filePath do
	puts "Sequence Alignment"
	`novoalign -d #{progSettings.novoIndex} #{setDataTypes.novoalign} -f #{sequence.filePath} > #{sequence.novoalignPath}`
	if $?.exitstatus != 0
		puts "Novoalign sequence alignment not performed (status = #{$?.exitstatus})"
		log.info("Novoalign sequence alignment not performed (status = #{$?.exitstatus})")
	else
		log.info("novoalign -d #{progSettings.novoIndex} #{setDataTypes.novoalign} -f #{sequence.filePath} > #{sequence.novoalignPath}")
	end
end

desc "Xnovotonm - Harvest non-base organism reads"
task :removeHuman => [:alignSequence, sequence.removeNonMappedPath]
file sequence.removeNonMappedPath => sequence.novoalignPath do
	puts "RemoveHuman"
	`#{PROG_DIR}/Xnovotonm.pl #{sequence.novoalignPath}`
	if $?.exitstatus != 0
		puts "Removal of Human mapped reads from novoalign results not performed (status = #{$?.exitstatus})"
		log.info("Removal of Human mapped reads from novoaign results not performed (status = #{$?.exitstatus})")
	else
		log.info(PROG_DIR + "Xnovotonm.pl #{sequence.removeNonMappedPath} - #{`egrep -cv '>' #{sequence.novoalignPath}`} total sequences reduced to #{`egrep -cv '>' #{sequence.removeNonMappedPath}`} after human reads removed.")
	end
end

desc "Tophat - Align spanning reads in order to remove base organism reads"
task :removeSpans => [:bowtieIndex, :tophatInstall, :removeHuman, sequence.blast1Path]
file sequence.blast1Path => sequence.removeNonMappedPath do
	puts "RemoveSpans"
	sh %{
		tophat -p #{ncpu} #{setDataTypes.tophat} --output-dir #{ENV['seq']}_tophat_out #{progSettings.bowtieIndex} #{sequence.filePath};
		samtools view -h -o #{ENV['seq']}_tophat_out/accepted_hits.sam #{ENV['seq']}_tophat_out/accepted_hits.bam;
		{PROG_DIR}/Xextractspans.pl "#{ENV['seq']}_tophat_out/accepted_hits.sam";
		{PROG_DIR}/Xfilterspans.pl "#{sequence.removeNonMappedPath}" "#{ENV['seq']}_tophat_out/accepted_hits.sam.spans";
	}
	if $?.exitstatus != 0
		puts "Removal of spans from source with Tophat not performed (status = #{$?.exitstatus})"
		log.info("Removal of spans from source with Tophat not performed (status = #{$?.exitstatus})")
	else
		log.info("tophat -p #{ncpu} #{setDataTypes.tophat} --output-dir #{ENV['seq']}_tophat_out #{progSettings.bowtieIndex} #{sequence.filePath}")
		log.info("samtools view -h -o #{ENV['seq']}_tophat_out/accepted_hits.sam #{ENV['seq']}_tophat_out/accepted_hits.bam")
		log.info("Xextractspans.pl #{ENV['seq']}_tophat_out/accepted_hits.sam")
		log.info("Xfilterspans.pl #{sequence.removeNonMappedPath} #{ENV['seq']}_tophat_out/accepted_hits.sam.spans")
	end
end

desc "BLAST - Associate reads with organisms."
task :localAlignReads => [:blastInstall, :ntInstall, :removeSpans, sequence.megan1Path]
file sequence.megan1Path => sequence.blast1Path do
	puts "localAlignReads"
	`blastn -db #{progSettings.ntDatabase} -soft_masking true -num_threads #{ncpu} -outfmt #{sequence.blastOutputFormat} -evalue #{sequence.eValue1} -query #{sequence.blast1Path} -out #{sequence.blast1Path}.blast`
	if $?.exitstatus != 0
		puts "BLAST of reads not performed (status = #{$?.exitstatus})"
		log.info("BLAST of reads not performed (status = #{$?.exitstatus})" )
	else
		#count of blast results
		log.info("blastn -db #{progSettings.ntDatabase} -soft_masking true -num_threads #{ncpu} -outfmt #{sequence.blastOutputFormat} -evalue #{sequence.eValue1} -query #{sequence.blast1Path} -out #{sequence.blast1Path}.blast")
	end
end

desc "MEGAN - Separate reads into taxonomies."
task :metaGenomeAnalyzeReads => [ :meganInstall, :localAlignReads, sequence.abyssPath]
file sequence.abyssPath => sequence.megan1Path do
	puts "metaGenomeAnalyzeReads"
	`MEGAN +g -x "import blastfile=#{sequence.megan1Path} readfile=#{sequence.blast1Path} meganfile=#{sequence.abyssPath} maxmatches=#{sequence.maxMatches} minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format=#{sequence.imageFileType} file=#{sequence.megan1Path + '.' + sequence.imageFileType.downcase}; quit;"`
	`MEGAN -f "#{sequence.abyssPath}" -x "#{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all;"`
	if $?.exitstatus != 0
		puts "MEGAN of BLASTed reads not performed (status = #{$?.exitstatus})"
		log.info("MEGAN of BLASTed reads not performed (status = #{$?.exitstatus})")
	else
		#LCA parameters used, number of statistically discernible sub groups, total number of matches over a threshold, number of not assigned reads, number of no hit reads.
		log.info("MEGAN +g -x import blastfile=#{sequence.megan1Path} readfile=#{sequence.blast1Path} meganfile=#{sequence.abyssPath} maxmatches=#{sequence.maxMatches} minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format=#{sequence.imageFileType} file=#{sequence.megan1Path + '.' + sequence.imageFileType.downcase}; quit;")
	end
end

desc "ABySS - Assemble reads associated with clusters of taxonomies."
task :denovoAssembleCluster => [:abyssInstall, :parallelIteratorInstall, :metaGenomeAnalyzeReads, FileList["#{sequence.blastPathGlob}"]]
file FileList["#{sequence.blastPathGlob}"] => FileList["#{sequence.abyssPathGlob}"] do
	puts "denovoAssemblyCluster"
	FileList["#{sequence.abyssPathGlob}"].each { | abyssFiles |
		`#{PROG_DIR}/abyssKmerOptimizer.pl #{abyssFiles} #{sequence.minKmerLength} #{sequence.maxKmerLength} #{setDataTypes.abyss}`
	}
	if $?.exitstatus != 0
		puts "DeNovo Assembly not performed (status = #{$?.exitstatus})"
		log.info("DeNovo Assembly not performed (status = #{$?.exitstatus})")
	else
		# log contig number, n50 score, contig average length, and winning kmer.
		log.info("#{PROG_DIR}/abyssKmerOptimizer.pl #{sequence.abyssPathGlob} #{sequence.minKmerLength} #{sequence.maxKmerLength} #{setDataTypes.abyss}")
	end
end

desc "BLAST - Associate contigs with organisms."
task :localAlignContigs => [:blastInstall, :ntInstall, :denovoAssembleCluster, FileList["#{sequence.megan2PathGlob}"]]
file FileList["#{sequence.megan2PathGlob}"] => FileList["#{sequence.blastPathGlob}"] do
	puts "localAlignContigs"
	FileList["#{sequence.blastPathGlob}"].each { | blastFiles |
		`blastn -db #{progSettings.ntDatabase} -soft_masking true -num_threads #{ncpu} -outfmt #{sequence.blastOutputFormat} -evalue #{sequence.eValue2} -query #{blastFiles} -out #{blastFiles}.blast`
	}
	if $?.exitstatus != 0
		puts "BLAST of contigs not performed (status = #{$?.exitstatus})" 
		log.info("BLAST of contigs not performed (status = #{$?.exitstatus})" )
	else
		#count of blast results
		log.info("blastn -db #{progSettings.ntDatabase} -soft_masking true -num_threads #{ncpu} -outfmt #{sequence.blastOutputFormat} -evalue #{sequence.eValue2} -query #{sequence.blastPathGlob} -out #{sequence.blastPathGlob}.blast")
	end
end

desc "MEGAN - Separate contigs into taxonomies."
task :metaGenomeAnalyzeContigs => [ :meganInstall, :localAlignContigs, FileList["#{sequence.pipeEndGlob}"]]
file FileList["#{sequence.pipeEndGlob}"] => FileList["#{sequence.megan2PathGlob}"] do
	puts "metaGenomeAnalyzeContigs"
	comparisonString = 'compare mode=absolute'
	FileList["#{sequence.blastPathGlob}"].each { | blastFiles |
		comparisonString += " meganfile=#{sequence.megan2PathGlob}"
		`MEGAN +g -x "import blastfile=#{sequence.megan2PathGlob} readfile=#{sequence.blastPathGlob} meganfile=#{pipeEndGlob} maxmatches=#{sequence.maxMatches} minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format=#{sequence.imageFileType} file=#{megan2PathGlob + '.' + sequence.imageFileType.downcase}; quit;"`
	}
	`MEGAN +g -x "#{comparisonString}"; save meganfile=#{seq}.rma`
	`MEGAN -f "#{seq}.rma"`
	if $?.exitstatus != 0
		puts "MEGAN of BLASTed contigs not performed (status = #{$?.exitstatus})"
		log.info("MEGAN of BLASTed contigs not performed (status = #{$?.exitstatus})")
	else
		#LCA parameters used, number of statistically discernible sub groups, total number of matches over a threshold, number of not assigned reads, number of no hit reads.
		log.info("MEGAN +g -x import blastfile=#{blastFiles + '.blast'} readfile=#{blastFiles} meganfile=#{blastFiles + '.blast.megan.rma'} maxmatches=#{sequence.maxMatches} minscore=#{sequence.minScoreByLength} toppercent=#{sequence.topPercent} winscore=#{sequence.winScore} minsupport=#{sequence.minSupport} summaryonly=false usecompression=true usecogs=#{sequence.useCogs} usegos=#{sequence.useGos} useseed=false; #{MEGAN_EXPANSION*sequence.expansionNumber} uncollapse all; update; exportgraphics format=#{sequence.imageFileType} file=#{blastFiles + '.blast.' + sequence.imageFileType.downcase}; quit;")
		log.info("MEGAN +g -x #{comparisonString}; save meganfile=#{seq}.rma")
	end
end



######### INSTALLATIONS

desc "Install latest version of human genome database."
task :hgInstall do
	if progSettings.humanGenomeDatabase.to_s.empty?
		progSettings.humanGenomeDatabase=File.dirname(`#{find} "chr*.fa" | head -1`)
		if progSettings.humanGenomeDatabase.to_s.empty?
			hg = ''
			ftp = Net::FTP.new("hgdownload.cse.ucsc.edu")
			ftp.login()
			ftp.passive=true
			ftp.chdir("apache/htdocs/goldenPath")
			ftp.list("hg[0-9][0-9]").last =~ /.*(hg\d+).*/
			hg = $1
			ftp.chdir(hg + "/bigZips")
			ftp.getbinaryfile("chromFa.tar.gz")
			ftp.close
			sh %{
				mkdir /usr/share/#{hg};
				tar -xzf chromFa.tar.gz -C /usr/share/#{hg};
				rm chromFa.tar.gz;
			}
			progSettings.humanGenomeDatabase = "/usr/share/" + hg
			puts progSettings.humanGenomeDatabase
		end
	end
end

desc "Create an index for novoalign of the human genome database."
task :novoIndex => [:hgInstall, :novoalignInstall] do
	if progSettings.novoIndex.to_s.empty?
		progSettings.novoIndex=`#{find} "*hg*.ndx" | head -1`
		if progSettings.novoIndex.empty?
			progSettings.novoIndex="#{progSettings.humanGenomeDatabase}/hgChrAll.ndx"
			`novoindex #{progSettings.novoIndex} #{progSettings.humanGenomeDatabase}/chr[0-9XY].fa #{progSettings.humanGenomeDatabase}/chr[0-9][0-9].fa`
		end
	end
end

desc "Create an index for bowtie/tophat of the human genome database."
task :bowtieIndex => [:hgInstall, :bowtieInstall] do
	if progSettings.bowtieIndex.to_s.empty?
		bowtieIndex=`#{find} "*hg*\.ebwt" | head -1`
		progSettings.bowtieIndex= $1 if bowtieIndex =~ /(.*hg.*)\.\d+\.ebwt/i
		resourceFiles=""
		if progSettings.bowtieIndex.to_s.empty?
			progSettings.bowtieIndex=progSettings.humanGenomeDatabase.to_s + '/hgChrAll'
			FileList[progSettings.humanGenomeDatabase.to_s + '/*.fa'].each do |filename|
				resourceFiles << $1 + "," if filename =~ /.*(chr[0-9XY]+.fa)/
			end
			resourceFiles.chomp!(',')
			sh %{
				cd #{progSettings.humanGenomeDatabase};
				bowtie-build #{resourceFiles} #{progSettings.bowtieIndex};
				`echo "export BOWTIEINDEX=#{progSettings.bowtieIndex}" >> ~/.#{shell}rc;`;
			}
		end
	end
end

desc "Install latest version of Novoalign."
task :novoalignInstall do
	if !command? "novoalign"
		novoalign = ""
		Net::HTTP.start('www.novocraft.com', 80) { |http|
			http.get('/main/releases.php', 'Referer' => 'http://www.novocraft.com/').body =~ /(V\d+\.\d+\.\d+)/
			novoalign = $1
			File.open('novocraft.tar.gz', 'w'){ |file|
				case os
				when OS::OSX then
					file.write(http.get('/downloads/' + novoalign + '/novocraft' + novoalign + '.MacOSX.tar.gz', 'Referer' => 'http://www.novocraft.com/').body)
				when OS::LINUX then
					file.write(http.get('/downloads/' + novoalign + '/novocraft' + novoalign + '.gcc.tar.gz', 'Referer' => 'http://www.novocraft.com/').body)
				end
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

desc "Install latest version of Bowtie."
task :bowtieInstall do
	if !command? "bowtie"
		bowtie = ""
		Net::HTTP.start('bowtie-bio.sourceforge.net', 80) { |http|
			http.get('/index.shtml').body =~ /https:\/\/sourceforge\.net\/projects\/bowtie-bio\/files\/bowtie\/(\d*\.\d*\.\d*)/
		}
		Net::HTTP.start('softlayer.dl.sourceforge.net', 80) { |http|
			bowtie = $1
			File.open('bowtie-' + bowtie + '-src.zip', 'w'){ |file|
				file.write(http.get('/project/bowtie-bio/bowtie/' + bowtie + '/bowtie-' + bowtie + '-src.zip').body)
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

desc "Install latest version of Samtools."
task :samtoolsInstall do
	if !command? "samtools"
		samtools = ""
		File.open('samtools.tar.bz2', 'w'){ |file|
			file.write(fetch('http://sourceforge.net/projects/samtools/files/latest').body)
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

desc "Install latest version of Tophat."
task :tophatInstall => [:samtoolsInstall, :bowtieInstall] do
	if !command? "tophat"
		tophat = ""
		Net::HTTP.start('tophat.cbcb.umd.edu', 80) { |http|
			http.get('/index.html').body =~ /\.\/downloads\/(tophat-\d+\.\d+\.\d+\.tar\.gz)/
			tophat = $1
			File.open(tophat, 'w'){ |file|
				file.write(http.get('/downloads/' + tophat).body)
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

desc "Install latest version of ABySS with Google Sparsehash."
task :abyssInstall do
	if !command? "ABYSS"
		Net::HTTP.start('code.google.com', 80) { |http|
			http.get('/p/google-sparsehash/').body =~ /(http:\/\/google-sparsehash\.googlecode\.com)(\/files\/sparsehash-\d+\.\d+\.tar\.gz)/
		}
		base = File.basename($1)
		gsh = $2
		puts gsh
		Net::HTTP.start(base, 80) { |http|
			File.open(File.basename(gsh), 'w'){ |file|
				file.write(http.get(gsh).body)
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

		abyss = ""
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

desc "Install latest version of BLAST+."
task :blastInstall do
	blast = ""
	if !command? "blastn"
		ftp = Net::FTP::new("ftp.ncbi.nlm.nih.gov")
		ftp.login()
		ftp.passive=true
		ftp.chdir("blast/executables/blast+/LATEST")
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

desc "Install latest version of the NT database."
task :ntInstall => :blastInstall do
	if ENV['BLASTDB'].to_s.empty? and progSettings.ntDatabase.to_s.empty?
		progSettings.ntDatabase=`#{find} "nt.nal" | head -1`.chomp('.nal')
		if !ENV['BLASTDB'].to_s.empty?
			progSettings.ntDatabase = ENV['BLASTDB'].to_s.chomp
		else
			sh %{
				mkdir /usr/share/nt;
				cd /usr/share/nt;
				/usr/bin/ncbi-blast*/c++/src/app/blast/update_blastdb.pl nt;
				for i in nt*.tar.gz; do tar -xzf $i; rm $i; done;
				`echo "export BLASTDB=/usr/share/nt/nt" >> ~/.#{shell}rc;`;
			}
		end
	end
end

desc "Install latest version of MEGAN."
task :meganInstall do
	if !command? "MEGAN"
		megan = ""
		Net::HTTP.start('www-ab.informatik.uni-tuebingen.de', 80) { |http|
			http.get('/data/software/megan/download/welcome.html').body =~ /(V\d+_\d+\/MEGAN_unix_\d+_\d+\.sh)/
			megan = $1
			File.open("#{File.basename(megan)}", 'w'){ |file|
				file.write(http.get('/data/software/megan/download/' + megan).body)
			}
		}
		megan = File.basename(megan)
		sh %{
			chmod +x #{megan};
			./#{megan};
		}
		if (arch == 64) # If CPU architecture is 64-bit, allow for more than 2GB of RAM and force 64-bit Java.
			text = File.new(`which MEGAN`.chomp).read.gsub(/"\$prg_dir\/\$progname" "-server" "-Xms\d+." "-Xmx\d+."/, "\"$prg_dir/$progname\" \"-server\" \"-d64\" \"-Xms#{memInGigs}G\" \"-Xmx#{memInGigs}G\"")
			File.open(`which MEGAN`.chomp, 'w+'){ |file| file.write(text) }
		end
	end
end

desc "Install latest version of Parallel::Iterator for perl."
task :parallelIteratorInstall do
	`perl -MParallel::Iterator -e 1`
	if $?.exitstatus != 0
		`perl -MCPAN -e 'install Parallel::Iterator'`
	end
end

task :default do
	Rake::Task[:metaGenomeAnalyzeContigs].invoke
end

desc "Install latest version of everything."
task :install do
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
	Rake::Task[:parallelIteratorInstall].invoke
end

# Reserialize object in case any changes have been made
YAML::dump(sequence, File.open(".#{seqName}", "w"))
YAML::dump(progSettings, File.open(File.expand_path("~/.#{PROG_NAME}"), "w"))