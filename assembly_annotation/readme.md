Scripts used for rockfish genome assembly

A) High quality whole genomes (we refer to them as PacBio based genomes or reference genomes throughout)
1) Basic preprocessing was performed for de novo PacBio assembly (check pacbio folder)
 - FALCON/FALCON_unzip : for five species aka Sebastes aleutianus, Sebastes miniatus, Sebastes pinniger, Sebastes rosaceus, Sebastes umbrosus
 - WTDBG2 : Performed for Sebastes entomelas (lower PacBio coverage), Sebastolobus alascanus (high heterozygosity)
2) Purging and polishing was performed for all the assemblies using PacBio and Illumina PE data
3) Bionano was used for only Sebastes umbrosus (please check Vertebrate Genome Project for another version of assembly)
4) HiC scaffolding was performed for five genomes - Sebastes aleutianus, Sebastes entomelas, Sebastes miniatus, Sebastes rosaceus, Sebastes umbrosus

B) Illumina only whole genomes
1) De novo assembly was performed by masurca
2) All 102 samples were de novo assembled
3) Scaffolding was performed using reference using RagTag - uploaded on Zenodo at https://zenodo.org/record/5534983

C) Mitogenome assembly and annotation
1) Illumina PE data was error corrected using musket
2) The error corrected Illumina PE data was used to asssemble 102 samples mitogenomes using MitoZ
3) Mitogenomes were annotated using MITOS - https://gitlab.com/Bernt/MITOS/

D) Quality checks
1) BUSCO was used for quality checks for all whole genome assemblies - different modes too, V3.1 and V4
2) Merqury was used for QC of PacBio based genomes

E) Annotation
1) PacBio based genomes were annotated for coding sequences using Funannotate. RNAseq was used for all samples except Sebastes aleutianus 
2) Braker2 annotation was performed for illumina assemblies (check the illumina/masurca_Snakefile) - but not used since we lifted annotations from reference genomes
3) Repeat annotation was performed using RepeatMasker
4) Noncoding annotations using tRNAscan-SE and infernal
5) Interproscan was used to search for functional domains and gene families in the reference genomes
6) For Illumina only genomes, LiftOff was used for annotations based on the reference genomes
7) The genome annotation GFF files are on Zenodo - https://zenodo.org/record/5534983
8) Annotation files for aleutianus and Sebastolobus are missing on zenodo, hence uploaded here - aleutianus_alascanus.gff3.zip
