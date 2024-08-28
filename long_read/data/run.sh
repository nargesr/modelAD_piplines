#!/bin/sh
#SBATCH -A model-ad_lab
#SBATCH --cpus-per-task 1
#SBATCH --output=kallisto.out
#SBATCH --error=kallisto.err
#SBATCH --time=1:00:00
#SBATCH -J kallisto_lr
#SBATCH --mail-type=START,END
#SBATCH --partition=standard

ref='/share/crsp/lab/model-ad/nargesr/kallisto_lr/ref_aug24/mm10'
ref_genome='/share/crsp/lab/model-ad/nargesr/kallisto_lr/ref_aug24/mm10.fa.gz'
ref_annot='/share/crsp/lab/model-ad/nargesr/kallisto_lr/ref_aug24/gencode.vM21.primary_assembly.annotation_UCSC_names.gtf.gz'
reads="/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD003_Aug24/rawfastq/fastqs"

path_to_kallisto='/share/crsp/lab/model-ad/nargesr/kallisto/build/src/kallisto'
path_to_bustools='/share/crsp/lab/model-ad/nargesr/bustools/build/src/bustools'

output='/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD003_Aug24/rawfastq/output'

#build index
#kb ref --kallisto ${path_to_lr_kallisto} -i ${ref}_k-63.idx -k 63 -f1 ${ref}.cdna.fa -g ${ref}.t2g.txt ${ref_genome} ${ref_annot}

#pseudoalign reads
COUNTER=1
IFS=$'\n'
for line in $(cat sample_name.txt)
do
    # Use awk to parse the columns from each line
    sample_name=$(echo $line | awk '{print $1}')
    fastq1=${reads}/$(echo $line | awk '{print $2}')
    fastq2=${reads}/$(echo $line | awk '{print $3}')
    printf "Sample: %s\n" "$sample_name"
    printf "Fastq1: %s\n" "$fastq1"
    printf "Fastq2: %s\n" "$fastq2"
    scriptName=${sample_name}
    curr=${sample_name}.sh
    output_sample=${output}_${sample_name}
    #fastq_file=${reads}/${sample}.fastq
    
    echo '#!/bin/bash' > ${curr}
    echo '#SBATCH -A model-ad_lab' >> ${curr}
    echo '#SBATCH --cpus-per-task 16' >> ${curr}
    echo '#SBATCH --output=kallisto_%J.out' >> ${curr}
    echo '#SBATCH --error=kallisto_%J.err' >> ${curr}
    echo '#SBATCH --time=02:00:00' >> ${curr}
    echo '#SBATCH -J kallisto_%J' >> ${curr}
    echo '#SBATCH --mail-type=START,END' >> ${curr}
    echo '#SBATCH --partition=standard' >> ${curr}

    echo "ref='/share/crsp/lab/model-ad/nargesr/kallisto_lr/ref_aug24/mm10'" >> ${curr}
    echo "reads='/share/crsp/lab/model-ad/nargesr/kallisto_lr/AD003_Aug24/rawfastq/fastqs'" >> ${curr}
    echo "path_to_kallisto='/share/crsp/lab/model-ad/nargesr/kallisto/build/src/kallisto'" >> ${curr}
    echo "path_to_bustools='/share/crsp/lab/model-ad/nargesr/bustools/build/src/bustools'" >> ${curr}
    echo "output='${output_sample}'" >> ${curr}
    
    echo "${path_to_kallisto} bus --long --threshold 0.8 -x bulk -i ${ref}_k-63.idx -o ${output_sample} ${fastq1} ${fastq1} ${name} -t 8" >> ${curr}

    echo "${path_to_bustools} sort -t 32 ${output_sample}/output.bus -o ${output_sample}/sorted.bus" >> ${curr} 
    echo "${path_to_bustools} count ${output_sample}/sorted.bus -t ${output_sample}/transcripts.txt  -e ${output_sample}/matrix.ec  -o ${output_sample}/count --cm -m -g ${ref}.t2g.txt" >> ${curr}

    echo "${path_to_kallisto} quant-tcc -t 32 --long -P ONT ${output_sample}/count.mtx -i ${ref}_k-63.idx -f ${output_sample}/flens.txt -e ${output_sample}/count.ec.txt -o ${output_sample}" >> ${curr}
    
    chmod +x ${curr}
    sbatch ${curr}

    COUNTER=`expr $COUNTER + 1`
    
done

