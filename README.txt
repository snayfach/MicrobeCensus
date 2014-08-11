DESCRIPTION
MicrobeCensus version 1.0.0 (28 April 2014): 
Rapidly and accurate estimate the average genome size (AGS) of a microbial community from metagenomic data. In short, AGS is estimated by aligning reads to a set of universal single-copy gene families. Because these genes occur in nearly all Bacteria and Archaea, genome size is inversely proportional to the fraction of reads which hit these genes.

AUTHORS: Stephen Nayfach (snayfach@gmail.com)

Usage: microbe_census [-options] <seqfile> <outfile>

Arguments:
  <seqfile>        Input metagenome. Can be a multi FASTA or FASTQ file.
                   Gzip (.gz) and Bzip (.bz2) file extensions recognized
  <outfile>        Tab delimited output file that includes AGS estimate.

Options:
  -h, --help       show this help message and exit
  -n NREADS        number of reads to use for AGS estimation (default = 1e6)
  -l READ_LENGTH   trim reads to this length (default = median read length)
  -f FILE_TYPE     file type: fasta or fastq (default = autodetect)
  -c QUAL_CODE     fastq quality score encoding: [sanger, solexa, illumina]
                   (default: autodetect)
  -t THREADS       number of threads to use for database search (default = 1)
  -q MIN_QUALITY   minimum base-level PHRED quality score: default = -5; no
                   filtering
  -m MEAN_QUALITY  minimum read-level PHRED quality score: default = -5; no
                   filtering
  -d               filter duplicate reads (default: False)
  -u MAX_UNKNOWN   max percent of unknown bases: default = 100%; no filtering
  -k               keep temporary files (default: False)

INSTALLATION
Simply download our software package from github:
git clone git@github.com:snayfach/MicrobeCensus.git
--OR--
wget https://github.com/snayfach/MicrobeCensus/archive/master.zip
And add the path of the src directory to your PATH.

REQUIREMENTS
x86_64 Architecture
Python 2.7.3 (Other versions of Python may be OK, but have not been tested)

RECOMMENDATIONS
* Filter duplicate reads using the -d flag.
  Be aware that this can consume large amounts of memory (>2G) when searching many reads (>20M)  
* Filter very low quality reads using -m 5 and -u 5.  
  Note that these options are only available for FASTQ files  
* Limit the number of reads searched (<nreads>) to less than 5e6.  
  We found that accurate estimates of AGS can be made using as few at 500K reads.  
* Be sure to remove potential sources of contamination from your metagenome, including  
  adaptor sequence and possibly host DNA (in the case of a host-associated metagenome).  

EXAMPLE USAGE
See the example directory.
