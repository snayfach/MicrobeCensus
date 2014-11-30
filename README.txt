DESCRIPTION
MicrobeCensus version 1.2.1 (11 November 2014): 
Rapidly and accurate estimate the average genome size (AGS) of a microbial 
community from metagenomic data. In short, AGS is estimated by aligning reads 
to a set of universal single-copy gene families. Because these genes occur in 
nearly all Bacteria and Archaea, genome size is inversely proportional to the 
fraction of reads which hit these genes.

AUTHORS: Stephen Nayfach (snayfach@gmail.com)

Usage: microbe_census.py [-options] <seqfile> <outfile>

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


REQUIREMENTS
Linux or Mac OSX with x86_64 Architecture
Python (tested on versions 2.6.6 and 2.7.3)
biopython
numpy


INSTALLATION
Option A. From PyPI:
Download our package
   wget https://pypi.python.org/packages/source/M/MicrobeCensus/MicrobeCensus-1.2.1.tar.gz
Unpack tarball
   tar -zxvf MicrobeCensus-1.2.1.tar.gz
Move to main directory
   cd MicrobeCensus-1.2.1
Install dependencies
   sudo python setup.py install
Add to MicrobeCensus to your path
   export PATH=$PATH:/your/install/dir/MicrobeCensus-1.2.1/src

Option B. From GitHub:
Download our package
   git clone https://github.com/snayfach/MicrobeCensus
Move to main directory
   cd MicrobeCensus
Install dependencies
   sudo python setup.py install
Add to MicrobeCensus to your path
   export PATH=$PATH:/your/install/dir/MicrobeCensus/src


USAGE RECOMMENDATIONS
- Filter duplicate reads using the -d flag.
  Be aware that this can consume large amounts of memory (>2G) when searching many reads (>20M)  
- Filter very low quality reads using -m 5 and -u 5.  
  Note that these options are only available for FASTQ files  
- Limit the number of reads searched (<nreads>) to less than 5M.  
  We found that accurate estimates of AGS can be made using as few at 300-500K reads. 
  Using more reads will produce slightly more accurate estimates of AGS, but will take more time to run.
- Be sure to remove potential sources of contamination from your metagenome, including  
  adaptor sequence and possibly host DNA (in the case of a host-associated metagenome).  
- More reads are better than longer reads. 
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
   avg_size - the average size of microbial genomes present in input metagenome


EXAMPLE USAGE
Input files: 
- MicrobeCensus/example/example.fq.gz contains 10,000 sequences in FASTQ format. Read lengths vary betweeb 60-100 bp. 
  Sequences are metagenomic reads from a stool sample.
- MicrobeCensus/example/example.fa.gz contains 10,000 500 bp sequences in FASTA format. 
  Sequences are simulated shotgun reads from the bacterial genome Treponema pallidum.
* These are toy datasets. In practice, between 300,000 to 500,000 reads are needed 
  for accurate estimates of average genome size for most metagenomes

Examples:
- Run with default options using either FASTA or FASTQ input files:
  microbe_census.py example.fq.gz fastq_example.out
  microbe_census.py example.fa.gz fasta_example.out
- Run with recommended quality filtering options:
  microbe_census.py -d -u 0 -m 5 example.fq.gz fastq_example.out
- Run with manually specified number of reads and read length:
  microbe_census.py -n 1000 -l 500 example.fa.gz fasta_example.out


DATA NORMALIZATION
We recommending using the statistic RPKG (reads per kb per genome equivalent) to quantify gene-family abundance from 
metagenomes.This is similar to the commonly used statistic RPKM, but instead of dividing by the number of mapped reads, 
we divide by the number of genome equivalents:
RPKG = (reads mapped to gene)/(gene length in kb)/(genome equivalents), and
genomes equivalent = (total DNA sequenced in bp)/(average genome size in bp), and
total DNA sequenced in bp = (read length in bp) * (reads sequenced)

Use case: 
We have two metagenomic libraries, L1 and L2, which each contain 1 million 100-bp reads:
READ_LENGTH_L1 = 100 bp
READS_SEQUENCED_L1 = 1,000,000 
TOTAL_DNA_L1 = 100,000,000 bp
READ_LENGTH_L2 = 100 bp
READS_SEQUENCED_L2 = 1,000,000 
TOTAL_DNA_L2 = 100,000,000 bp

We use MicrobeCensus to estimate the average genome size of each library: 
AGS_L1 = 2,500,000 bp
AGS_L2 = 5,000,000 bp

Next, we map reads from each library to a reference database which contains a gene of interest G. G is 1000 bp long. 
We get 100 reads mapped to gene G from each library:
LENGTH_G = 1,000 bp
MAPPED_READS_G_L1 = 100
MAPPED_READS_G_L2 = 100

Finally, we quantify RPKG for gene G in each library:
RPKG for G in L1 = (100 mapped reads)/(1 kb)/(100,000,000 bp sequenced / 2,500,000 bp AGS) = 2.5   
RPKG for G in L2 = (100 mapped reads)/(1 kb)/(100,000,000 bp sequenced / 5,000,000 bp AGS) = 5.0


TRAINING MICROBECENSUS FOR NEW MARKER GENES
The TRAINING.txt file contains detailed information on how to train MicrobeCensus using different marker genes
and/or training genomes.


TESTING MICROBECENSUS
See MicrobeCensus/test/TEST.txt for details.





