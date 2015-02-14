# MicrobeCensus
MicrobeCensus is a fast and easy to use pipeline for estimating the average genome size (AGS) of a microbial 
community from metagenomic data. 

In short, AGS is estimated by aligning reads to a set of universal single-copy gene families present in nearly all cellular microbes (Bacteria, Archaea, Fungi). 
Because these genes occur once per genome, the average genome size of a microbial community is inversely proportional to the fraction of reads which hit these genes.

### Requirements
* Python dependencies (installed via setup.py): Numpy, BioPython
* Supported platforms: Mac OSX, Unix/Linux; Windows not currently supported
* Python version 2 or 3

### Installation
Download MicrobeCensus from: https://github.com/snayfach/MicrobeCensus/archive/v1.0.2.tar.gz  

Unpack the project: 
`tar -zxvf MicrobeCensus-1.0.2.tar.gz`

Navigate to the installation directory:  
`cd /path/to/MicrobeCensus`  

Run setup.py. This will install any dependencies:  
`python setup.py install` or  
`sudo python setup.py install` to install as a superuser

Alternatively, MicrobeCensus can be installed directly from PyPI:  
`pip install MicrobeCensus` or   
`sudo pip install MicrobeCensus` to install as a superuser

### Using MicrobeCensus without installing
Although this is not recommended, users may wish to run MicrobeCensus without running setup.py.  

Both BioPython and Numpy will both need to be already installed.
You should be able to enter the following command in the python interpreter without getting an error:  
`>>> import Bio.SeqIO`  
`>>> import numpy`

Next, add the MicrobeCensus module to your PYTHONPATH environmental variable:  
`export PYTHONPATH=$PYTHONPATH:/path/to/MicrobeCensus` and  
`echo export PYTHONPATH='$PYTHONPATH':/path/to/MicrobeCensus >> ~/.bash_profile` to avoid entering the command in the future

Finally, add the scripts directory to your PATH environmental variable:  
`export PATH=$PATH:/path/to/MicrobeCensus/scripts` and  
`echo export PATH='$PATH':/path/to/MicrobeCensus/scripts >> ~/.bash_profile` to avoid entering the command in the future

Now, you should be able to enter the command into your terminal without getting an error:  
`run_microbe_census.py -h`

### Testing the software
After installing MicrobeCensus, we recommend testing the software:  
`cd /path/to/MicrobeCensus/test`  
`python test_microbe_census.py`

### Running MicrobeCensus
MicrobeCensus can either be run as a command-line script or imported to python as a module.

#### Command-line usage
**run_microbe_census.py [-options] seqfile outfile**

arguments:
* **seqfile**: path to input metagenome. can be FASTA/FASTQ formatted. gzip (.gz) and bzip2 (.bz2) compression supported.
* **outfile**: path to output file

options: 
* **-h, --help**: show this help message and exit 
* **-n NREADS**: number of reads to use for AGS estimation (default = 1e6)  
* **-l READ_LENGTH**: trim reads from 3' to this length (default = median read length of seqfile)  
* **-f FILE_TYPE {fasta,fastq}**: FASTA or FASTQ formatted seqfile (default = autodetect)
* **-c QUAL_ENCODE {fastq-sanger,fastq-solexa,fastq-illumina}**: Quality encoding for FASTQ files (default = autodetect)
* **-t THREADS**: number of threads to use for database search (default = 1)  
* **-q MIN_QUALITY**: minimum base-level PHRED quality score (default = -5)  
* **-m MEAN_QUALITY**: minimum read-level PHRED quality score (default = -5)  
* **-d**: filter duplicate reads (default = False)  
* **-u MAX_UNKNOWN**: max percent of unknown bases (default = 100)  
* **-k**: keep temporary files (default = False)  
* **-s**: suppress printing program's progress to stdout (default = False)

#### Module usage

First, import the module:  
`>>> from microbe_census import microbe_census`

Next, setup your options and arguments, formatted as a dictionary. The path to your metagenome is the only requirement (default values will be used for all other options):  
`>>> args = {'seqfile':'MicrobeCensus/microbe_census/example/example.fq.gz'}`

Alternatively, other options can be specified:
```
>>> args = {
  'seqfile':'MicrobeCensus/microbe_census/example/example.fq.gz',
  'nreads':100000,
  'read_length':100,
  'file_type':'fastq',
  'quality_type':'fastq-sanger',
  'threads':1,
  'min_quality':10,
  'mean_quality':10,
  'filter_dups':False,
  'max_unknown':0,
  'quiet':False}
```

Finally, the entire pipeline can be run by passing your arguments to the run_pipeline function. MicrobeCensus returns the estimated AGS of your metagenome, along with a dictionary of used arguments:
`average_genome_size, args = microbe_census.run_pipeline(args)`

#### Recommended options
* Use -n to limit the number of reads searched. Suggested values are between 500,000 and 1 million. Using more reads may result in slightly more accurate estimates of AGS, but will take more time to run.
* Remove potential sources of contamination from your metagenome. This may include: adaptor sequences, host DNA, or viral DNA.  
* Generally it is better to use more reads rather than longer reads.
* Filter very low quality reads using -m 5 and -u 5.

### Software speed
* Run times are for a 150 bp library. Expect longer/shorter runtimes depending on read length.

Threads (-t)  | Reads/Second
------------- | -------------
1  | 830
2  | 1,300
4  | 1,800
8  | 2,000


### Normalization
AGS estimated from MicrobeCensus can be used to normalize functional abundance profiles. We recommending using the statistic RPKG (reads per kb per genome equivalent) to quantify gene-family abundance from shotgun metagenomes. This is similar to the commonly used statistic RPKM, but instead of dividing by the number of total mapped reads, we divide by the number of genome equivalents:

>RPKG = (reads mapped to gene)/(gene length in kb)/(genome equivalents), where  
>genomes equivalent = (total DNA sequenced in bp)/(average genome size in bp), and  
>total DNA sequenced in bp = (read length in bp) * (reads sequenced)

Use case: 
We have two metagenomic libraries, L1 and L2, which each contain 1 million 100-bp reads:
>READ_LENGTH_L1 = 100 bp  
>READS_SEQUENCED_L1 = 1,000,000   
>TOTAL_DNA_L1 = 100,000,000 bp  
>READ_LENGTH_L2 = 100 bp  
>READS_SEQUENCED_L2 = 1,000,000   
>TOTAL_DNA_L2 = 100,000,000 bp  

We use MicrobeCensus to estimate the average genome size of each library: 
>AGS_L1 = 2,500,000 bp  
>AGS_L2 = 5,000,000 bp

Next, we map reads from each library to a reference database which contains a gene of interest G. G is 1000 bp long. 
We get 100 reads mapped to gene G from each library:
>LENGTH_G = 1,000 bp  
>MAPPED_READS_G_L1 = 100  
>MAPPED_READS_G_L2 = 100

Finally, we quantify RPKG for gene G in each library:
>RPKG for G in L1 = (100 mapped reads)/(1 kb)/(100,000,000 bp sequenced / 2,500,000 bp AGS) = 2.5  
>RPKG for G in L2 = (100 mapped reads)/(1 kb)/(100,000,000 bp sequenced / 5,000,000 bp AGS) = 5.0  

### Training
We have included scripts and documentation for retraining MicrobeCensus, using user-supplied training genomes and gene families. Documentation and scripts can be found under: MicrobeCensus/training

### Citing
If you use MicrobeCensus, please cite:  

Average genome size estimation improves comparative metagenomics and sheds light on the functional ecology of the human microbiome. **Stephen Nayfach and Katherine Pollard**. _Genome Biology (in review)_.

