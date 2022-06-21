# Transcriptome assembly pipeline

This pipeline describes the process of assembling transcriptomes from RNA-seq data, as well as the evaluation and optimisation of these transcriptomes. It additionally describes the annotation and quantification of these transcriptomes. Finally, it generates a DESeq object as an input for DESeq2, which may be used for differential expression analysis. 

List of tools utilised for this pipeline: 
1. SPAdes-rna
2. Trinity
3. Salmon 
4. BBmap
5. rnaQUAST
6. BUSCO
7. Dammit
8. DESeq2 in R

## 1. Assembling the transcriptomes
Little information was available regarding the best assembler to use, therefore the two top candidates were tested. These assemblers were SPAdes-rna and Trinity

### SPAdes-rna
`for i in path/to/samples/*.fastq; 
do
spades.py –rna -s $i -o “$i”_sp
done
`

### Trinity
`for i in /path/to/samples/*.fastq; 
do
Trinity --seqType fq --max_memory 200G --single $i --output /dev/shm/trinity_out_dir/"$i"_trinity --CPU 28 --full_cleanup
done`

## 2. Evaluation of transcriptomes
Transcriptomes were evaluated using rnaQUAST, BUSCO, and tools native to the Trinity suite

### rnaQUAST
`/path/to/rnaQUAST/rnaQUAST.py --transcripts /path/to/folder/*.fasta \
--reference /path/to/reference/transcriptome.fa \
--gtf /path/to/annotation/file.gtf \
-o /path/to/output/dir --busco /path/to/busco/database \
--gene_mark --disable_infer_genes --disable_infer_transcripts
`

### BUSCO
`ls *.fasta | while read -r line; 
do;
docker run -u $(id -u) -v /Users/student/:/busco_wd ezlabgva/busco:v5.2.2_cv1 busco -i genome.fna -m transcriptome -i path/to/dir/$line -o  $line -l \ 
eudicots_odb10 --cpu 6 -f;
done`

### Trinity suite tools
I.	Get quant.sf files using Trinity:
```for i in *.fasta; 
do 
base=`basename $i .fasta`; 
~/trinityrnaseq-v2.12.0/util/align_and_estimate_abundance.pl\
--transcripts /path/to/transcripts/directory \
--seqType fa \ 
--single $i --est_method salmon --thread_count 16 --trinity_mode\
--output_dir ./"$base"; 
done

II.	Abundance estimate to matrix:

 ~/trinityrnaseq-v2.12.0/util/abundance_estimates_matrix.pl\
 --est_method salmon --gene_trans_map none --name_sample_by_basedir \
 --out_prefix Tr $(find . -name "*.sf")

III.	Ex90N50: 

 ~/trinityrnaseq-v2.12.0/util/misc/contig_ExN50_statistic.pl mock.isoform.TMM.EXPR.matrix\
/path/to/transcripts |\
 tee ExN50.stats  ExN50.txt


## 3. Transcriptome optimisation

1. Combining the results of two assemblers
`for i in ~/NAS/output/spades_transcripts/*.fasta;
do
cat $i ~/NAS/output/trinity_transcripts/$i >> "$i"_combined.fasta
done`

Alternative: CD-HIT-est

2. Deduplication and clumping of reads to optimise combined transcriptome
`for i in path/to/combined/files/*.fasta;
do
~/bbmap/dedupe.sh in=“$i” out= “$i”_co.fasta ac=f
~/bbmap/clumpify.sh in=“$i”_co.fasta out=“$i”_clumped.fasta 
done`

## 4. Annotation of transcriptome

`dammit annotate ~/NAS/output/transcripts/co_supertranscript.fasta
--busco-group embryophyta --user-databases
~/NAS/databases/eukaryota_odb10/refseq_db.faa --n_threads 20 --force`

## 5. Quantification

1. Indexing

`salmon index -t transcriptome.fa -i {out_name}`

2. Quantification

`salmon quant -i index -l A -r transcript.fasta --validateMappings -o path/to/out/dir`

## 6. Preparing input for differential expression analysis

`R
library(tximport)
library(GenomicFeatures)
library(readr)`

`txdb <- makeTxDbFromGFF("/path/to/annotation/file.gff")
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- select(txdb, keys = k, keytype = "TXNAME", columns = "GENEID")
samples <- read.csv(‘/path/to/samples/file.csv’, header=TRUE)
files <- file.path(samples$Sample, "quant.sf")`

`names(files) <- paste0(samples$Sample)
txi <- tximport(files, type="salmon", tx2gene = tx2gene)
DS2 <- DESeqDataSetFromTximport(txi, samples, ~Condition)
`


```python

```
