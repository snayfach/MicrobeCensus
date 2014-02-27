seqfile='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/example/example.fq.gz'
outfile='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/example/example.test'
nreads=100000
read_length=50
min_quality=-5
file_type='fastq'
fastq_code='solexa'
nboot=1
microbe_census='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/src/microbe_census.py'
python $microbe_census -n $nboot -c $fastq_code -t 1 -q $min_quality -f $file_type $seqfile $outfile $nreads $read_length

