# copy rap 2.15 database
SRC='/home/snayfach/projects/mgs_090713/markers/dwu/db/dwu_30_bact_arch.2.15.rapdb'
DEST='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/data/rapdb_2.15'
cp $SRC $DEST

SRC='/home/snayfach/projects/mgs_090713/markers/dwu/db/dwu_30_bact_arch.2.15.rapdb.info'
DEST='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/data/rapdb_2.15.info'
cp $SRC $DEST

# copy rap 2.15 src
SRC='/home/snayfach/packages/RAPSearch2.15_64bits/bin/rapsearch'
DEST='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/src/rapsearch_2.15'
cp $SRC $DEST


# run test
seqfile='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/example/example.fq'
outfile='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/example/example.test'
nreads=100000
read_length=50
min_quality=-5
max_unknown=2
file_type='fastq'
fastq_code='solexa'
microbe_census='/home/snayfach/software/MicrobeCensus/MicrobeCensus_v0.1/src/microbe_census.py'
python $microbe_census -u $max_unknown -d -c $fastq_code -t 1 -f $file_type $seqfile $outfile $nreads $read_length






