# **ERVmap**
ERVmap is one part curated database of human proviral ERV loci and one part a stringent algorithm to determine which ERVs are transcribed in their RNA seq data.

## Citation 
Tokuyama M. et. al., ERVmap analysis reveals genome-wide transcription of human endogenous retroviruses. Proc Natl Acad Sci USA 2018 Dec 11;115(50):12565-12572. [doi: 10.1073/pnas.1814589115](http:/doi.org/10.1073/pnas.1814589115).

## **How to use it**

### Install
This version of the tool consists on 2 steps: 1. alignment to the human genome (GRC38) and 2. quantification of the ERV regions. To download and install ERVmap latest version provided as docker image, simply type:
```
docker pull eipm/ervmap:latest
```
**NOTE**: for a specific version replace `latest` with the release version. 

### **How to run ERVmap**
To run ERVmap, you'd need: 1. an indexed genome reference for STAR; 2. A bed file with the curated ERV regions on the human genome (see `ERVmap.bed`); 3. the input FASTQ data (gzipped).  Assuming that your sample is called `SAMPLE`, and has 2 FASTQ files (one per read) in the folder `/path/to/input/data`; the reference genome is in `/path/to/genome` and the ERV bed file is in `/path/to/erv/file` here is the command:
```
docker run --rm  \
    -u $(id -u):$(id -g) \
    -v /path/to/input/data:/data:ro \
    -v /path/to/genome:/genome:ro \
    -v /path/to/erv/file:/resources:ro \
    -v /path/to/output:/results \
    ervmap \
    --read1 /data/SAMPLE_1.fastq.gz \
    --read2 /data/SAMPLE_2.fastq.gz \
    --output SAMPLE/SAMPLE. \
    --mode ALL
```
This command will generate the alignment files (BAMs) in the `/path/to/output/SAMPLE/` folder and all files will have the prefix `SAMPLE.`. The generated files will be:
```
SAMPLE.Aligned.sortedByCoord.out.bam
SAMPLE.Aligned.sortedByCoord.out.bam.bai
SAMPLE.ERVresults.txt
SAMPLE.Log.final.out
SAMPLE.Log.out
SAMPLE.Log.progress.out
SAMPLE.SJ.out.tab
```
(See [STAR documentation](https://github.com/alexdobin/STAR) for the description of the output files of the STAR aligner ). 
The results of ERV quantification will be in the `SAMPLE.ERVresults.txt` file. This is a tab-delimited file with 7 columns from  [bedtools](https://bedtools.readthedocs.io/en/latest/). For example:
```
1       896176  898458  5803    500     +       70
1       1412251 1418852 5804    500     +       36
1       3801730 3806808 5807    500     +       6
1       4178468 4187573 5808    500     +       1
```

## The **`--mode`** option
This option can only have 3 values: { `ALL`, `STAR`, `BED` }:
* `ALL` to run both the STAR aligner and the ERV quantification from start to finish; 
* `STAR` to only perform the alignment;
* `BED` to only run the ERV quantification.


### Optional parameters (recommended) 
There are a few parameters that can be added to the ERVmap image to make the process more efficient.
* `--cpus 20`: if you have a multi-core system (and you should have one), you can specify the number of CPUs to use (e.g. 20);
* `--limit-ram 48000000000`: this limits the amount of RAM used to avoid overusing the resources 
You can see the full set of parameters by typing: `docker run --rm ervmap`.

There are also other parameters from Docker that should be included before `ervmap` in the command line, e.g. 
```    
    --memory 50G \
    --memory-swap 100G
``` 

----

# Published version 

Please note that the instructions hereafter refer to the orignal published version (see [ERVmap on GitHub](https://github.com/mtokuyama/ERVmap))

## **Installing**

### Install dependencies
``` 
bedtools2
cufflinks
bwa-0.7.17
cufflinks-2.2.1.Linux_x86_64
python
samtools-1.8
tophat-2.1.1.Linux_x86_64
tophat2
trim (http://graphics.med.yale.edu/trim/)
```

### Install .pl and r files
```
erv_genome.pl
interleaved.pl
run_clean_htseq.pl
clean_htseq.pl
merge_count.pl
normalize_with_file.pl
normalize_deseq.r
```

## **Map data to human genome (hg38)**

This step will yield raw counts for cellular genes and ERVmap loci as separate files.

### For single-end sequences:
```
erv_genome.pl -stage 1 -stage2 6 -fastq /${i}_SS.fastq.gz
```

### For pair-end sequences:
```
interleaved.pl --read1  ${i}_R1.fastq.gz  --read2 ${i}_R2.fastq.gz > ${i}.fastq.gz
erv_genome.pl -stage 1 -stage2 6 -fastq /${i}.fastq.gz
```

### Store output files
```
mkdir -p output
mv ./sample/herv_coverage_GRCh38_genome.txt ./output/erv/${i}.e
mv ./sample/GRCh38/htseq.cnt ./output/cellular/${i}.c
```

## **Clean up data, merge, and normalize**

These steps will yield normalized ERV read counts based on size factors obtained through DESeq2 analysis. 
Use the output files from above. 

```
run_clean_htseq.pl ./output/cellular c c2 __
merge_count.pl 3 6 e ./output/erv > ./output/erv/merged_erv.txt
merge_count.pl 0 1 c2 ./output/cellular > ./output/cellular/merged_cellular.txt
normalize_deseq.r  ./output/cellular/merged_cellular.txt ./output/cellular/normalized_cellular ./output/cellular/normalized_factors
normalize_with_file.pl ./output/cellular/normalized_factors ./output/erv/merged_erv.txt > ./output/$folder_name.txt
```

## Authors

* Maria Tokuyama
* Yong Kong



