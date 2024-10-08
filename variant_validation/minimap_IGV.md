# Map reads using minimap2

Please make sure that you already installed necessary packages.

## Making custom references

In order to map the reads, you need to prepare references using [this notebook](make_custom_GTF_FASTA.ipynb)

## Map reads using minimap2

This script will combine fastq files for each sample and then generate and submit separate jobs for each sample to run minimap on bulk long read RNA_seq data to sam files for each sample and then convert them to bam file, followed by sorting and indexing them to be able to visualize them in IGV: 

```bash
#!/bin/sh
#SBATCH -A model-ad_lab
#SBATCH --cpus-per-task 16
#SBATCH --output=minimap.out
#SBATCH --error=minimap.err
#SBATCH --time=5:00:00
#SBATCH -J minimap
#SBATCH --mail-type=START,END
#SBATCH --partition=standard
#SBATCH --mem-per-cpu=6G

fasta_ref='path_to_reference_genome' #'/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD014_Sep_mm39/ref/MAD1.fa'
reads='path_to_reads_file' #'/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD014_Sep_mm39/fastqs'

#pseudoalign reads
COUNTER=1
IFS=$'\n'
for line in $(cat ../sample_name.txt)
do
    # Use awk to parse the columns from each line
    sample_name=$(echo $line | awk '{print $1}')
    fastqs_num=$(echo $line | awk '{print $2}')
    fastqs=()
    for i in $(seq 1 $fastqs_num);
    do
      	fastq_index=$((i+2))
        fastqs+=(${reads}/$(echo $line | awk -v idx=$fastq_index '{print $idx}'))
    done
    printf "Sample: %s\n" "$sample_name"
    printf "Fastqs: %s\n" "${fastqs[@]}"

    # Combine all fastqs into a single string with spaces
    fastqs_string=$(printf "%s " "${fastqs[@]}")

    # Remove trailing space
    fastqs_string=$(echo $fastqs_string | sed 's/ *$//')
    printf "cat $fastqs_string > ${reads}/${sample_name}.fastq.gz\n"
    cat "${fastqs[@]}" > ${reads}/${sample_name}.fastq.gz

    scriptName=${sample_name}
    curr=${sample_name}.sh
    output_sample=${output}_${sample_name}
    fastq_file=${reads}/${sample_name}.fastq.gz
    
    echo '#!/bin/bash' > ${curr}
    echo '#SBATCH -A model-ad_lab' >> ${curr}
    echo '#SBATCH --cpus-per-task 16' >> ${curr}
    echo '#SBATCH --output=minimap_%J.out' >> ${curr}
    echo '#SBATCH --error=minimap_%J.err' >> ${curr}
    echo '#SBATCH --time=05:00:00' >> ${curr}
    echo '#SBATCH -J minimap_%J' >> ${curr}
    echo '#SBATCH --mail-type=START,END' >> ${curr}
    echo '#SBATCH --partition=standard' >> ${curr}
    echo '#SBATCH --mem-per-cpu=6G' >> ${curr}


    echo "module load minimap2/2.24" >> ${curr}
    echo "module load samtools/1.10" >> ${curr}
    
    echo "fasta_ref='/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD014_Sep_mm39/ref/MAD1.fa'" >> ${curr}
    echo "reads='/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD014_Sep_mm39/fastqs'" >> ${curr}
    
    echo "minimap2 --MD -t 16 -ax splice -k14 ${fasta_ref} ${reads}/${sample_name}.fastq.gz > ${sample_name}.sam 2> ${sample_name}.log" >> ${curr}

    echo "samtools view -b ${sample_name}.sam -o ${sample_name}.bam" >> ${curr}
    echo "samtools sort ${sample_name}.bam -o ${sample_name}_sorted.bam" >> ${curr} 
    echo "samtools index ${sample_name}_sorted.bam ${sample_name}_sorted.bam.bai" >> ${curr}

    chmod +x ${curr}
    sbatch ${curr}

    COUNTER=`expr $COUNTER + 1`

done
```

`sample_name.txt` used in the script contains the name of the sample followed by the number of fastq files we have for the samples and then the corresponding fastq files in each line and it should look like this:

```
ad003_11616_lig-blk   2   ad003_11616_lig-blk_1.fastq.gz   ad003_11616_lig-blk_2.fastq.gz
ad003_11617_lig-blk   2   ad003_11617_lig-blk_1.fastq.gz   ad003_11617_lig-blk_2.fastq.gz
ad003_11625_lig-blk   2   ad003_11625_lig-blk_1.fastq.gz   ad003_11625_lig-blk_2.fastq.gz
ad003_11627_lig-blk   2   ad003_11627_lig-blk_1.fastq.gz   ad003_11627_lig-blk_2.fastq.gz
ad003_11628_lig-blk   2   ad003_11628_lig-blk_1.fastq.gz   ad003_11628_lig-blk_2.fastq.gz
ad003_11629_lig-blk   2   ad003_11629_lig-blk_1.fastq.gz   ad003_11629_lig-blk_2.fastq.gz
ad003_12517_lig-blk   2   ad003_12517_lig-blk_1.fastq.gz   ad003_12517_lig-blk_2.fastq.gz
ad003_12648_lig-blk   2   ad003_12648_lig-blk_1.fastq.gz   ad003_12648_lig-blk_2.fastq.gz
ad003_12649_lig-blk   2   ad003_12649_lig-blk_1.fastq.gz   ad003_12649_lig-blk_2.fastq.gz
ad003_12659_lig-blk   2   ad003_12659_lig-blk_1.fastq.gz   ad003_12659_lig-blk_2.fastq.gz
ad003_12660_lig-blk   2   ad003_12660_lig-blk_1.fastq.gz   ad003_12660_lig-blk_2.fastq.gz
ad003_12670_lig-blk   2   ad003_12670_lig-blk_1.fastq.gz   ad003_12670_lig-blk_2.fastq.gz
```

## Visualize reads using IGV

Once you have the sorted bam file along with BAI index, you can look at the individual reads through IGV. 
You should also use the same GTF and FASTA files that were used to map the reads.


In the end, the structure of your files would be something like this:

````
.
sample_name.txt
minimap/
├── ad014-1_26241_mux-blk.bam
├── ad014-1_26241_mux-blk.log
├── ad014-1_26241_mux-blk.sam
├── ad014-1_26241_mux-blk.sh
├── ad014-1_26241_mux-blk_sorted.bam
├── ad014-1_26241_mux-blk_sorted.bam.bai
├── ad014-1_26242_mux-blk.bam
├── ad014-1_26242_mux-blk.log
├── ad014-1_26242_mux-blk.sam
├── ad014-1_26242_mux-blk.sh
├── ad014-1_26242_mux-blk_sorted.bam
├── ad014-1_26242_mux-blk_sorted.bam.bai
├── ad014-1_26243_mux-blk.bam
├── ad014-1_26243_mux-blk.log
├── ad014-1_26243_mux-blk.sam
├── ad014-1_26243_mux-blk.sh
├── ad014-1_26243_mux-blk_sorted.bam
├── ad014-1_26243_mux-blk_sorted.bam.bai
├── ad014-1_26244_mux-blk.bam
├── ad014-1_26244_mux-blk.log
├── ad014-1_26244_mux-blk.sam
├── ad014-1_26244_mux-blk.sh
├── ad014-1_26244_mux-blk_sorted.bam
├── ad014-1_26244_mux-blk_sorted.bam.bai
├── ad014-1_26245_mux-blk.bam
├── ad014-1_26245_mux-blk.log
├── ad014-1_26245_mux-blk.sam
├── ad014-1_26245_mux-blk.sh
├── ad014-1_26245_mux-blk_sorted.bam
├── ad014-1_26245_mux-blk_sorted.bam.bai
├── ad014-1_26246_mux-blk.bam
├── ad014-1_26246_mux-blk.log
├── ad014-1_26246_mux-blk.sam
├── ad014-1_26246_mux-blk.sh
├── ad014-1_26246_mux-blk_sorted.bam
├── ad014-1_26246_mux-blk_sorted.bam.bai
├── ad014-1_26247_mux-blk.bam
├── ad014-1_26247_mux-blk.log
├── ad014-1_26247_mux-blk.sam
├── ad014-1_26247_mux-blk.sh
├── ad014-1_26247_mux-blk_sorted.bam
├── ad014-1_26247_mux-blk_sorted.bam.bai
├── ad014-1_26248_mux-blk.bam
├── ad014-1_26248_mux-blk.log
├── ad014-1_26248_mux-blk.sam
├── ad014-1_26248_mux-blk.sh
├── ad014-1_26248_mux-blk_sorted.bam
├── ad014-1_26248_mux-blk_sorted.bam.bai
├── ad014-1_26249_mux-blk.bam
├── ad014-1_26249_mux-blk.log
├── ad014-1_26249_mux-blk.sam
├── ad014-1_26249_mux-blk.sh
├── ad014-1_26249_mux-blk_sorted.bam
├── ad014-1_26249_mux-blk_sorted.bam.bai
├── ad014-1_26263_mux-blk.bam
├── ad014-1_26263_mux-blk.log
├── ad014-1_26263_mux-blk.sam
├── ad014-1_26263_mux-blk.sh
├── ad014-1_26263_mux-blk_sorted.bam
├── ad014-1_26263_mux-blk_sorted.bam.bai
├── ad014-1_26264_mux-blk.bam
├── ad014-1_26264_mux-blk.log
├── ad014-1_26264_mux-blk.sam
├── ad014-1_26264_mux-blk.sh
├── ad014-1_26264_mux-blk_sorted.bam
├── ad014-1_26264_mux-blk_sorted.bam.bai
├── ad014-1_26265_mux-blk.bam
├── ad014-1_26265_mux-blk.log
├── ad014-1_26265_mux-blk.sam
├── ad014-1_26265_mux-blk.sh
├── ad014-1_26265_mux-blk_sorted.bam
├── ad014-1_26265_mux-blk_sorted.bam.bai
├── minimap.err
├── minimap.out
└── run_minimap.sh
ref/
├── APOE.fa
├── APOE.fa.fai
├── gencode.v46.chr_patch_hapl_scaff.annotation.gtf
├── gencode.v46.chr_patch_hapl_scaff.annotation.gtf.gz
├── gencode.vM32.chr_patch_hapl_scaff.annotation.gtf
├── gencode.vM32.chr_patch_hapl_scaff.annotation.gtf.gz
├── hg38.fa
├── hg38.fa.fai
├── hg38.fa.gz
├── kallisto_index_LR.err
├── kallisto_index_LR.out
├── MAD1.cdna.fa
├── MAD1.fa
├── MAD1.fa.fai
├── MAD1.gtf
├── MAD1_k-63.idx
├── MAD1.t2g.txt
├── make_custom_GTF_FASTA.ipynb
├── MAPT.fa
└── mm39.fa
````
