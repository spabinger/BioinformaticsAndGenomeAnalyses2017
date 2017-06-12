# Practical

In this practical you will get to know basic tools for SAM/BAM manipulation and call variants using different programs. Please check out the [Useful information](#useful-information) section.

## Required tools

* Bedtools - http://bedtools.readthedocs.io
* Freebayes - https://github.com/ekg/freebayes
* GATK - https://www.broadinstitute.org/gatk/download
* IGV - https://www.broadinstitute.org/igv/home
* Java (7) - http://www.oracle.com/technetwork/java/javase/downloads/index.html?ssSourceSiteId=otnjp
* Picard - http://broadinstitute.github.io/picard/
* SAMtools - http://samtools.sourceforge.net/‎
* VCFtools - https://vcftools.github.io/index.html
* VCFlib - https://github.com/vcflib/vcflib
* Tabix for vcf-tools - https://sourceforge.net/projects/samtools/files/tabix/



## Information
* Original data from [1000 genomes](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data/NA10847/exome_alignment/).
* Stop GATK "calling home" - [information](http://gatkforums.broadinstitute.org/gatk/discussion/1250/what-is-phone-home-and-how-does-it-affect-me)  


## Exercise

We assume that we have a properly mapped BAM file from quality checked reads.
For some variant callers we use a [target file](target.bed) to shorten variant calling time.

#### Important

* After each step inspect the generated output (cat, less, head, grep, ...).
* Organize your data and results in folders.
* Check out the specified parameters. If you don't know the exact meaning, consult the help pages (-h, --help).



#### SAMtools

__(*)__ Unzip all files for the practicals
    tar xyfz pabinger_files_for_pracicals.tar.gz -C ~


__(*)__ How big is the BAM file

    ls -lah <file.bam>
    

__(*)__ Inspect the header of the BAM file

    module add samtools
    samtools ...
    samtools view -H aln.bam


__(*)__ View the BAM file
    
    samtools view <bam.file> | less
    
__(*)__ How many reads are in the BAM file? <br/>
Is there another way to count the reads (check the samtools flagstat)
   
    samtools view <file.bam> | grep -v "^#" | wc -l
    
    
__(*)__ Answer the following questions by investigating the SAM file
* Print only the SAM header
* What version of the human assembly was used to perform the alignments?
* What version of bwa was used to align the reads?
* What is the name of the first read?
* At what position does the alignment of the read start?
```
    Use SAMtools for these questions
    samtools view -H aln.bam | less
    samtools view aln.bam | less
```
    
__(*)__ Sort the BAM file

    samtools sort -o sorted.bam -@ <THREADS> <file.bam>
    
__(*)__ Extract small file

    samtools view -h sorted.bam | head -n 50000 | samtools view -S -b -0 small_sorted.bam 
    
__(*)__ Index the bam file
    
    samtools index <sorted.bam>
        

#### Alignment stats
    samtools flagstat sorted.bam
    samtools idxstats sorted.bam


  
#### Prepare reference genome
__(*)__ Prepare dict index
    
    module add picard
    java -jar -Xmx2G /bcga2016/picard-tools-2.2.1/picard.jar CreateSequenceDictionary R=hg19.fasta O=hg19.dict

__(*)__ Prepare fai index
    
    samtools faidx hg19.fasta 


#### BAM file preparations
__(*)__ Sort with Picard
    
    java -Xmx2g -jar /bcga2016/picard-tools-2.2.1/picard.jar SortSam I=aln.bam O=sorted_picard.bam SORT_ORDER=coordinate


__(*)__ Mark duplicates
     
    java -Xmx2g -jar /bcga2016/picard-tools-2.2.1/picard.jar MarkDuplicates I=sorted_picard.bam O=dedup.bam M=metrics.txt


__(*)__ Add ReadGroup
    
    java -Xmx2g -jar /bcga2016/picard-tools-2.2.1/picard.jar AddOrReplaceReadGroups I=dedup.bam O=deduprg.bam RGID=group1 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=sample1


__(*)__ Index with Picard
    
    java -Xmx2g -jar /bcga2016/picard-tools-2.2.1/picard.jar BuildBamIndex I=deduprg.bam O=deduprg.bam.bai VALIDATION_STRINGENCY=SILENT

__(*)__ Collect insert size metrics
    
    module add R-3.2.4
    java -Xmx2g -jar /bcga2016/picard-tools-2.2.1/picard.jar CollectInsertSizeMetrics I=deduprg.bam O=insertSizeHistogram.txt H=insertSizeHistogram.pdf
    
__(*)__ View the PDF
    evince insertSizeHistogram.pdf


__(*)__ Questions
* How many reads were marked as duplicated? Look for flags.
* What are the other sorting possibilities for SortSam?
* Inspect the insert size metrics histogram.




#### GATK variant calling

__(*)__ Known indel sites are here specified as variables - either copy the whole path or use variables as well

    KNOWN_INDELS_1="1000G_phase1.indels.hg19.vcf"
    KNOWN_INDELS_2="Mills_and_1000G_gold_standard.indels.hg19.vcf"


__(*)__ Realignment target creator

    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T RealignerTargetCreator -R hg19.fasta -nt 8 -L target.bed \
    -I deduprg.bam -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o target_intervals.list

__(*)__ Perform realignment
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T IndelRealigner -R hg19.fasta -I deduprg.bam \
    -targetIntervals target_intervals.list -known ${KNOWN_INDELS_1} -known ${KNOWN_INDELS_2} -o dedup_rg_real.bam


__(*)__ Base quality recalibration
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam \
    -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o recal_data_table.txt -L target.bed 
    --maximum_cycle_value 800


__(*)__ Second pass of recalibration
     
     java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T BaseRecalibrator -R hg19.fasta -I dedup_rg_real.bam \
     -knownSites ${KNOWN_INDELS_1} -knownSites ${KNOWN_INDELS_2} -o post_recal_data_table.txt \
     -L target.bed --maximum_cycle_value 800 -BQSR recal_data_table.txt 


__(*)__ Generate before after plots (requires R and ggplot2)

    Fix missing R packages
    R
    Inside R call
    install.packages(c('reshape','gplots','gsalib'))
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T AnalyzeCovariates -R hg19.fasta -L target.bed 
    -before recal_data_table.txt -after post_recal_data_table.txt -plots recalibration_plots.pdf


__(*)__ Print recalibrated reads
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T PrintReads -R hg19.fasta -L target.bed -I dedup_rg_real.bam 
    -BQSR recal_data_table.txt -o dedup_rg_real_recal.bam


__(*)__ Now do variant calling
    
    java -Xmx8g -jar /bcga2016/GATK-3.5/GenomeAnalysisTK.jar -T HaplotypeCaller -R hg19.fasta -nct 8 -L target.bed 
    -I dedup_rg_real_recal.bam --genotyping_mode DISCOVERY -o gatk.vcf

__(*)__ Questions
* Check out the before and after plots.
* How many variants were called?



#### Display files in IGV

    (Download and open) IGV
    Load the BAM file and the VCF files into IGV
    Look at the mapping on Chr 11
    Check out the results of the different variant calling programs.



#### BONUS: FreeBayes variant calling

__(*)__ Call

     module add freebayes-1.0.2
     freebayes -f hg19.fasta -t target.bed -v freebayes.vcf deduprg.bam

__(*)__ Investigate result
  
    #Perform the same procedures as done for samtools
    #Do you notice differences?
    grep -v "^#" freebayes.vcf | wc -l



#### BONUS: VCF stats

__(*)__ VCFlib - stats - shown here for one VCF file

    vcfstats freebayes.vcf > freeb_call_vcf_lib_stats.txt



    
#### Useful information

__(*)__ Determine number of cores

    cat /proc/cpuinfo  
    
__(*)__ Determine memory size
    
    cat /proc/meminfo

__(*)__ Make file executable

    chmod +x <file.name>
    
__(*)__ Extract information from a file

    grep -i "info" <file.name>
    
__(*)__ Extract information from a file excluding the pattern and display the first 100 lines

    grep -v "^#" <file.name> | head -n 100










