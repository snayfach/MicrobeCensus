DESCRIPTION
MicrobeCensus version 1.0.0 (28 April 2014): 
Rapidly and accurate estimate the average genome size (AGS) of a microbial community from metagenomic data. 
In short, AGS is estimated by aligning reads to a set of universal single-copy gene families. Because these 
genes occur in nearly all Bacteria and Archaea, genome size is inversely proportional to the fraction of 
reads which hit these genes. 

A preprint of our manuscript can be found in bioRxiv:
http://biorxiv.org/content/biorxiv/early/2014/09/11/009001.full.pdf

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
                   reads shorter than this length will be discarded
                   supported values include: 
                     50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 175, 
                     200, 225, 250, 300, 350, 400, 450, 500
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
Download our software package from github:
git clone git@github.com:snayfach/MicrobeCensus.git
--OR--
wget https://github.com/snayfach/MicrobeCensus/archive/master.zip

Add MicrobeCensus to your PATH:
export PATH=$PATH:/path/to/MicrobeCensus/src

REQUIREMENTS
Linux or Mac OSX with x86_64 Architecture
Python (tested on versions 2.6.6 and 2.7.3)

RECOMMENDATIONS
* Filter duplicate reads using the -d flag.
  Be aware that this can consume large amounts of memory (>2G) when searching many reads (>20M)  
* Filter very low quality reads using -m 5 and -u 5.  
  Note that these options are only available for FASTQ files  
* Limit the number of reads searched (<nreads>) to less than 5M.  
  We found that accurate estimates of AGS can be made using as few at 300-500K reads. 
  Using more reads will produce slightly more accurate estimates of AGS, but will take more time to run.
* Be sure to remove potential sources of contamination from your metagenome, including  
  adaptor sequence and possibly host DNA (in the case of a host-associated metagenome).  
* More reads are better than longer reads. 
  While longer reads generally produce more accurate estimates of AGS, there are diminishing returns once reads 
  are longer than ~150 bp. On the other hand, at least 300-500K reads are needed for an accurate estimate of AGS.

SOFTWARE SPEED
1 CPU: ~830   reads/second
2 CPU: ~1,300 reads/second
4 CPU: ~1,800 reads/second
8 CPU: ~2,000 reads/second

* Benchmarking was performed on a dataset containing 100 bp reads. Longer reads will have lower throughput.
* Software speed may depend on your system

OUTPUT FILE FORMAT
Tab delimited file with the following fields:
reads_sampled - the number of trimmed reads that were sampled from the input metagenome
read_length - all reads were trimmed to this length
avg_size - the size of microbial genomes present in input metagenome

EXAMPLE USAGE
The EXAMPLE.txt file contains more detailed information on running MicrobeCensus with various options.

CORRECTING FOR AGS
The NORMALIZATION.txt file contains detailed information on how to correct for AGS in your data.

