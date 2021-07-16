#!/usr/bin/bash

set -e

\ls -1d illumina/fastq/Sebastes_aleutianus/S-aleutianus_SEB-111/*.gz >sebale_sr_fastq.fofn
\ls -1d pacbio/S_aleutianus/05.28.2019/Seb111_Cell*/*fastq.gz >sebale_pb_fastq.fofn

# Run purging - we used the first version of purge dups here
module load purge_dups/1.0
pd_config.py -l ./ -s sebale_sr_fastq.fofn -n config.Fsebale1.PB.asm1.json falcon/run1/4-polish/cns-output/cns_p_ctg.fasta sebale_pb_fastq.fofn
python ~/local_modules_sw/purge_dups/scripts/run_purge_dups.py config.Fsebale1.PB.asm1.json ~/local_modules_sw/purge_dups/bin sebale1_FLC -p bash -w 1

# Create alternative assembly with the alternate haplotigs from FALCON and the redundant sequences from purgingf
cat falcon/run1/4-polish/cns-output/cns_h_ctg.fasta purging/cns_p_ctg/seqs/cns_p_ctg.red.fa >squirrel2.hap.fasta

