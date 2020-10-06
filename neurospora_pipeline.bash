########################SNP CALLING#############################
# popolation_1.2.2 package copied to home directory without installation
# bwa-0.7.17 installed locally as in the manual and added it to $PATH, so it could be called by "bwa"
# samtools 1.9 installed locally (./configure --prefix="$HOME") as in the manual and added it to $PATH, so it could be called by "samtools"

# copied all files in alexg/neurospora

# removing spaces from names in the reference genome and putting nucl&mito together
awk '{print $1}' ~/neurospora/neurospora_crassa_or74a_12_supercontigs.fasta > ~/neurospora/neurospora_crassa_or74a_12_supercontigs_shortheader.fasta
awk '{print $1}' ~/neurospora/neurospora_crassa_or74a_mito_10_supercontigs.fasta > ~/neurospora/neurospora_crassa_or74a_mito_10_supercontigs_shortheader.fasta
cat ~/neurospora/neurospora_crassa_or74a_12_supercontigs_shortheader.fasta ~/neurospora/neurospora_crassa_or74a_mito_10_supercontigs_shortheader.fasta > ~/neurospora/or74a_full.fasta

#counting nt per scaffold in the reference neurospora or74a genome version NC12 + mitochondrial supercontig_10.21
awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next;} {seqlen += length($0)} END {print seqlen}' ~/neurospora/or74a_full.fasta

# >Supercontig_12.1
# 9798893
# >Supercontig_12.2
# 4478683
# >Supercontig_12.3
# 5274802
# >Supercontig_12.4
# 6000761
# >Supercontig_12.5
# 6436246
# >Supercontig_12.6
# 4218384
# >Supercontig_12.7
# 4255303
# >Supercontig_12.8
# 192308
# >Supercontig_12.9
# 142473
# >Supercontig_12.10
# 125404
# >Supercontig_12.11
# 31696
# >Supercontig_12.12
# 19714
# >Supercontig_12.13
# 13515
# >Supercontig_12.14
# 11565
# >Supercontig_12.15
# 9397
# >Supercontig_12.16
# 8983
# >Supercontig_12.17
# 6701
# >Supercontig_12.18
# 6309
# >Supercontig_12.19
# 4755
# >Supercontig_12.20
# 1646
# >Supercontig_10.21
# 64840

# UNIVERSAL PIPELINE, change strain

#combining two fastq files with reads into one for both mate pairs. Unfortunately have to do manually
cat /mnt/nfs/bioinfdata/home/NIOO/alexg/neurospora/9t2/FCH3YN5BBXX-wHAXPI028005-20_L7_1.fq /mnt/nfs/bioinfdata/home/NIOO/alexg/neurospora/9t2/FCH3YN5BBXX-wHAXPI028005-20_L8_1.fq > ~/neurospora/concatenated/9t2_F.fq
cat /mnt/nfs/bioinfdata/home/NIOO/alexg/neurospora/9t2/FCH3YN5BBXX-wHAXPI028005-20_L7_2.fq /mnt/nfs/bioinfdata/home/NIOO/alexg/neurospora/9t2/FCH3YN5BBXX-wHAXPI028005-20_L8_2.fq > ~/neurospora/concatenated/9t2_R.fq

#count the number of reads before trimming and put counts in a separate file
for filename in ~/neurospora/concatenated/*.*
do
echo $filename >> ~/neurospora/stats/stats.txt
awk '(NR+3)%4==0 {sum+=1} END {print "number of raw reads: "sum}' $filename >> ~/neurospora/stats/stats.txt
awk 'NR%4==2 {sum+=length($1)} END {print "total number of raw nucleotides: "sum}' $filename >> ~/neurospora/stats/stats.txt
done

# /home/NIOO.INT/alexg/neurospora/concatenated/10t1_F.fq
# number of raw reads: 4872543
# total number of raw nucleotides: 730881450
# /home/NIOO.INT/alexg/neurospora/concatenated/10t1_R.fq
# number of raw reads: 4872543
# total number of raw nucleotides: 730881450
# ...

# trimming neurospora strains (allow minimum read length of 70 bp, with quality not lower than 30)
# a loop going through raw reads directory, selecting files with "F", extracting sample name from them, adding "_F" or "_R" in the end, and trimming PE reads in pair of files and storing them in "/trimmed" directory
for filename in ~/neurospora/concatenated/*_F.*
do
SAMPLE=$(ls $filename | awk -F'[_/]' '{print$7}')
FL1=$(ls ~/neurospora/concatenated/${SAMPLE}_F.fq)
FL2=$(ls ~/neurospora/concatenated/${SAMPLE}_R.fq)
echo $FL1
echo $FL2
perl ~/popoolation_1.2.2/basic-pipeline/trim-fastq.pl --input1 ${FL1} --input2 ${FL2} --min-length 70 --quality-threshold 30 --output ~/neurospora/trimmed/${SAMPLE}_trimmed.fq --fastq-type sanger
echo "Sample ${SAMPLE} done"
done

# select trimmed files ending with "1" and "2", count the number of reads after trimming and append to stat.txt
for filename in ~/neurospora/trimmed/*_trimmed.fq_[1,2]*
do
echo $filename >> ~/neurospora/stats/stats.txt
awk '(NR+3)%4==0 {sum+=1} END {print "number of trimmed reads: "sum}' $filename >> ~/neurospora/stats/stats.txt
awk 'NR%4==2 {sum+=length($1)} END {print "total number of trimmed nucleotides: " sum}' $filename >> ~/neurospora/stats/stats.txt
done

# /home/NIOO.INT/alexg/neurospora/trimmed/10t1_trimmed.fq_1
# number of trimmed reads: 3861750
# total number of trimmed nucleotides: 560515414
# /home/NIOO.INT/alexg/neurospora/trimmed/10t1_trimmed.fq_2
# number of trimmed reads: 3861750
# total number of trimmed nucleotides: 505961216
# ...

#indexing the reference genome
bwa index ~/neurospora/or74a_full.fasta

#extracting sample names again and a loop with aligning to the reference genome; we use default settings for bwa mem (parameter -t is the number of cores used, parameter $2 is sample name parameter $3 is mapping quality cutoff). This creates *.sam file, which is an alignment map
for filename in ~/neurospora/trimmed/*.params
do
SAMPLE=$(ls $filename | awk -F'[_/]' '{print$7}')
FL1=$(ls ~/neurospora/trimmed/${SAMPLE}_trimmed.fq_1)
FL2=$(ls ~/neurospora/trimmed/${SAMPLE}_trimmed.fq_2)
echo $FL1
echo $FL2
bwa mem -t 12 ~/neurospora/or74a_full.fasta ${FL1} ${FL2} > ~/neurospora/sams/${SAMPLE}.sam
echo "Sample ${SAMPLE} done"
done

## creating a sorted bam file (filtered alignment) with mapping qualities higher than 20 (it is a good cut off)
for filename in ~/neurospora/sams/*.sam
do
SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
echo $SAMPLE >> ~/neurospora/stats/stats.txt
samtools view -q 20 -b -S ~/neurospora/sams/${SAMPLE}.sam | samtools sort -o ~/neurospora/bams/${SAMPLE}_q20.sort.bam
echo "samtools flagstat ${SAMPLE}" >> ~/neurospora/stats/stats.txt
samtools flagstat ~/neurospora/bams/${SAMPLE}_q20.sort.bam >> ~/neurospora/stats/stats.txt
samtools index ~/neurospora/bams/${SAMPLE}_q20.sort.bam
echo "samtools depth ${SAMPLE}" >> ~/neurospora/stats/stats.txt
samtools depth ~/neurospora/bams/${SAMPLE}_q20.sort.bam | awk '{sum+=$3} END {print "average coverage bam file, coverage regions = ",sum/NR}' >> ~/neurospora/stats/stats.txt
echo "Sample ${SAMPLE} done"
done

# 10t1
# samtools flagstat 10t1
# 7659617 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 secondary
# 10656 + 0 supplementary
# 0 + 0 duplicates
# 7659617 + 0 mapped (100.00% : N/A)
# 7648961 + 0 paired in sequencing
# 3824820 + 0 read1
# 3824141 + 0 read2
# 7527550 + 0 properly paired (98.41% : N/A)
# 7648165 + 0 with itself and mate mapped
# 796 + 0 singletons (0.01% : N/A)
# 27619 + 0 with mate mapped to a different chr
# 27619 + 0 with mate mapped to a different chr (mapQ>=5)
# samtools depth 10t1
# average coverage bam file, coverage regions =  25.9672
# ...

# # NOW PAIRWISE FST ANALYSES AFTER WE PREPARED SORTED/FILTERED BAM FILES --->> fst method not used !!!
# # Now we can compare ancestor to the evolved strains
# # first we create an mpileup (between ancestor (Nc159) and evolved genomes on the reference background). -B means no BAQ correction, -f next reference file in fasta

# for filename in ~/neurospora/bams/[1-8]t*_q20.sort.bam
# do
# SAMPLE=$(ls $filename | awk -F'[_/]' '{print$7}')
# echo $SAMPLE
# samtools mpileup -Bf ~/neurospora/or74a_full.fasta ~/neurospora/bams/Nc159_q20.sort.bam $filename > ~/neurospora/mpileups/Nc159_${SAMPLE}.mpileup
# ## then a popoolation2 sync file
# java -jar /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/mpileup2sync.jar --input ~/neurospora/mpileups/Nc159_${SAMPLE}.mpileup --output ~/neurospora/syncs/Nc159_${SAMPLE}.sync --threads 12
# ## then we calculate the pairwise fst (fst value of 1 means completely divergent SNPs between the two genomes) for every nucleotides covered, play with minimum coverage, we set to 5. But also try 4, 3, and see how many fst==1 we'll get. Compare.
# perl /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/fst-sliding.pl --input ~/neurospora/syncs/Nc159_${SAMPLE}.sync --output ~/neurospora/fsts/Nc159_${SAMPLE}_mincov5.fst --suppress-noninformative --min-count 2 --min-coverage 5 --max-coverage 2500 --window-size 1 --step-size 1 --pool-size 30
# ## then we make a output that sorts the fst values
# sed 's/1:2=//g' ~/neurospora/fsts/Nc159_${SAMPLE}_mincov5.fst | sort -rk6,6 > ~/neurospora/fsts/Nc159_${SAMPLE}_mincov5.fst.txt
# echo "Sample ${SAMPLE} done"
# done

# # doing the same for wild-type ancestor (Nc152) and respective evolved genomes (except 9)
# for filename in ~/neurospora/bams/[1][0-6]t*_q20.sort.bam
# do
# SAMPLE=$(ls $filename | awk -F'[_/]' '{print$7}')
# samtools mpileup -Bf ~/neurospora/or74a_full.fasta ~/neurospora/bams/Nc152_q20.sort.bam $filename > ~/neurospora/mpileups/Nc152_${SAMPLE}.mpileup
# ## then a popoolation2 sync file
# java -jar /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/mpileup2sync.jar --input ~/neurospora/mpileups/Nc152_${SAMPLE}.mpileup --output ~/neurospora/syncs/Nc152_${SAMPLE}.sync --threads 12
# ## then we calculate the pairwise fst (fst value of 1 means completely divergent SNPs between the two genomes) for every nucleotides covered, play with minimum coverage, we set to 5. But also try 4, 3, and see how many fst==1 we'll get. Compare.
# perl /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/fst-sliding.pl --input ~/neurospora/syncs/Nc152_${SAMPLE}.sync --output ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst --suppress-noninformative --min-count 2 --min-coverage 5 --max-coverage 2500 --window-size 1 --step-size 1 --pool-size 30
# ## then we make a output that sorts the fst values
# sed 's/1:2=//g' ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst | sort -rk6,6 > ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst.txt
# echo "Sample ${SAMPLE} done"
# done

# # same for the line 9 (because I could not code it in one go with the previous for loop) 
# for filename in ~/neurospora/bams/[9]t*_q20.sort.bam
# do
# SAMPLE=$(ls $filename | awk -F'[_/]' '{print$7}')
# samtools mpileup -Bf ~/neurospora/or74a_full.fasta ~/neurospora/bams/Nc152_q20.sort.bam $filename > ~/neurospora/mpileups/Nc152_${SAMPLE}.mpileup
# ## then a popoolation2 sync file
# java -jar /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/mpileup2sync.jar --input ~/neurospora/mpileups/Nc152_${SAMPLE}.mpileup --output ~/neurospora/syncs/Nc152_${SAMPLE}.sync --threads 12
# ## then we calculate the pairwise fst (fst value of 1 means completely divergent SNPs between the two genomes) for every nucleotides covered, play with minimum coverage, we set to 5. But also try 4, 3, and see how many fst==1 we'll get. Compare.
# perl /mnt/nfs/bioinfdata/home/NIOO/alexg/popoolation2_1201/fst-sliding.pl --input ~/neurospora/syncs/Nc152_${SAMPLE}.sync --output ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst --suppress-noninformative --min-count 2 --min-coverage 5 --max-coverage 2500 --window-size 1 --step-size 1 --pool-size 30
# ## then we make a output that sorts the fst values
# sed 's/1:2=//g' ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst | sort -rk6,6 > ~/neurospora/fsts/Nc152_${SAMPLE}_mincov5.fst.txt
# echo "Sample ${SAMPLE} done"
# done

# # a complicated sorting loop for FSTs
# for filename in ~/neurospora/fsts/*_mincov5.fst.txt #take each filename
# do
# awk '{print $6}' $filename | sort -ru > ${filename}.col6 #extract column 6, take only uniqe values from 1 to 0 and store in a temp file
	# for row in `cat ${filename}.col6` #take each row (each fst value) of that temp file
	# do
	# grep $row $filename | sort -k1 >> $filename.sorted #return lines having a particular fst value in the original file, and sort them by 1 column, and append
	# done
# rm ${filename}.col6 #remove temp file each cycle
# done

# #once we have a nice total FST file, we can pick the SNPs we want: e.g. FST 1.00000000 & COVERAGE >=10
# for filename in ~/neurospora/fsts/*.sorted
# do
# echo $filename
# awk '$5 >= 10 && $6 >= 0.85' $filename
# done

#indexing original bams just in case (why not?)
for filename in ~/neurospora/bams/*.bam
do
	[[ $filename == *_q20.sort.bam ]] && continue # if $filename matches *_q20.sort.bam
	[ -e "$filename" ] || continue #exclude it
	SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $filename
echo $SAMPLE
samtools index ~/neurospora/bams/${SAMPLE}.sort.bam
done



##########VARSCAN FOR SNPS AND SMALL INDELS##########################
# since fst analysis gave the same results for SNPs as varscan, we can only use varscan for SNPs and small indels
# using prevously generated pairvise mpileupfiles and feed it to VarScan v 2.4.4
# generating two vcf files: with --min-var-freq 0.5 and --min-var-freq 0.8 --> using more stringent cutoff 0.8 in further analysis
for filename in ~/neurospora/mpileups/*.mpileup
do
SAMPLE=$(ls $filename | awk -F '[./]' '{print$8}')
echo $SAMPLE
java -jar ~/varscanv2.4.4/VarScan.v2.4.4.jar mpileup2cns $filename --min-coverage 10 --min-var-freq 0.5 --p-value 0.005 --variants --strand-filter 0 --output-vcf 1 >  ~/neurospora/vcf/${SAMPLE}_freq05.vcf
echo "${SAMPLE} done"
done

# R script to process vcf files
library("vcfR") ## from bioconductor

   # *****       ***   vcfR   ***       *****
   # This is vcfR 1.8.0 
     # browseVignettes('vcfR') # Documentation
     # citation('vcfR') # Citation
   # *****       *****      *****       *****
   
setwd("P:/!!!postdoc/vcf/")
path = "P:/!!!postdoc/vcf/"
file.names<-dir(path, pattern ="9t1_freq08.vcf$") #select vcf files
#par(mfrow=c(6,6)) #creat area for the plots

for(i in 1:length(file.names)){
NAME<-as.matrix(unlist(strsplit(file.names[i],"_")))
NAME<-NAME[2,1]
vcf<-read.vcfR(file.names[i])
AD<-extract.gt(vcf, element='AD', as.numeric=TRUE)
DP<-extract.gt(vcf, element='DP', as.numeric=TRUE) 

#AD contains of all samples the alternative depth, those reads that are different from reference
#DP is the total depth

AF<-AD/DP

#AF is essentially the allele frequency, filtering on coverage >=10

head(AF)
L<-which(is.na(DP[,1])==FALSE & is.na(DP[,2])==FALSE & DP[,1]>=10 & DP[,2]>=10)
dAF<-abs(AF[,1]-AF[,2])
plot(dAF[L])

#filtering on >0.8 ratio
write.table(cbind(NAME,vcf@fix[L[which(dAF[L]>0.8)],],AF[L[which(dAF[L]>0.8)],],DP[L[which(dAF[L]>0.8)],],AD[L[which(dAF[L]>0.8)],]),sep='\t', paste0(getwd(),"/snps1.txt"),col.names=F, row.names=F, quote=F, append=TRUE)
}

p<-read.table("P:/!!!postdoc/vcf/snps.txt",header=F)
plot(p$V1)




##################BIG INDELS SEARCH, script from Joost##########
### here I assume there is a sort.bam, which is completely raw, without any filtering
#making unfiltered bams from previous sams (no filtering by -q20 flag)
for filename in ~/neurospora/sams/*.sam
do
SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
samtools view -b -S ~/neurospora/sams/${SAMPLE}.sam | samtools sort -o ~/neurospora/bams/${SAMPLE}.sort.bam
done

#select unfiltered bams and pick weird reads (with various flags)
for filename in ~/neurospora/bams/*.sort.bam
do 
	[[ $filename == *_q20.sort.bam ]] && continue # if $filename matches *_q20.sort.bam
	[ -e "$filename" ] || continue #exclude it
SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
samtools view -h $filename | awk '$6 ~ /I/ || $6 ~ /D/ || $6 ~ /S/ || $6 ~ /H/ || $2==67 || $2==131 || $2==115 || $2==179 || $2==81 || $2==161 || $2==97 || $2== 145 || $2==65 || $2==129 || $2==113 || $2==177 ||$1 ~ /^@/' | samtools view -bS - > ~/neurospora/struct_var/$SAMPLE.structvar.bam
done

#indexing bams
for filename in ~/neurospora/struct_var/*.structvar.bam
do 
SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
samtools index ~/neurospora/struct_var/${SAMPLE}.structvar.bam
done

# then we have to quantify the coverage at each base for the raw bams and the structural variation bams
#relative to Nc152 (manually remove later unnecessary combinations)
for filename in ~/neurospora/bams/*.sort.bam
do
	[[ $filename == *_q20.sort.bam ]] && continue # if $filename matches *_q20.sort.bam
	[ -e "$filename" ] || continue #exclude it
	SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
samtools depth -a ~/neurospora/bams/Nc152.sort.bam ~/neurospora/struct_var/Nc152.structvar.bam  ~/neurospora/bams/$SAMPLE.sort.bam ~/neurospora/struct_var/$SAMPLE.structvar.bam > ~/neurospora/struct_var/Nc152_$SAMPLE.depthfile.txt
done
#relative to Nc159 (manually remove later unnecessary combinations)
for filename in ~/neurospora/bams/*.sort.bam
do
	[[ $filename == *_q20.sort.bam ]] && continue # if $filename matches *_q20.sort.bam
	[ -e "$filename" ] || continue #exclude it
	SAMPLE=$(ls $filename | awk -F'[./]' '{print$8}')
echo $SAMPLE
samtools depth -a ~/neurospora/bams/Nc159.sort.bam ~/neurospora/struct_var/Nc159.structvar.bam  ~/neurospora/bams/$SAMPLE.sort.bam ~/neurospora/struct_var/$SAMPLE.structvar.bam > ~/neurospora/struct_var/Nc159_$SAMPLE.depthfile.txt
done

#R start
# now the depth file .txt is a file that contains the read depth for all reads and then for all suspicious reads
# so therefore we use this loop R script to calculate the depth in bins (*0.01) indicate the binsize is 1/0.01 = 100 bps which is a good cutoff, as reads are roughly 100 bp
# smaller bin will indicate more interesting stuff.

setwd("P:/!!!postdoc/struct_var/")
path = "P:/!!!postdoc/struct_var/"
out.file<-""
file.names<-dir(path, pattern =".depthfile.txt")
for(i in 1:length(file.names))
		{
		dat<-read.table(file.names[i], header=F)

# dat<-data[grep("^Supercontig_12.1$", data$V1), ]
# head(dat)
# tail(dat)

		WINDOWS<-round(dat$V2*0.01)
		D<-dim(dat)[[2]]
		CHRS<-unique(dat$V1)
		for (j in 1:length(CHRS)){
			L<-which(dat$V1==CHRS[j])
			OUT<-unique(WINDOWS[L])
			for (y in 3:D){
			OUT<-cbind(OUT,tapply(dat[L,y], WINDOWS[L], mean))
			}
		M<-length(unique(WINDOWS[L]))
		write.table(cbind(rep(as.character(dat$V1)[L[1]],M),round(OUT,2)), paste0(getwd(),"/",file.names[i],".binned100bp.txt"), sep='\t', col.names=F, row.names=F, quote=F, append=TRUE)
			}
		}


for(i in 1:length(file.names))
		{
		dat<-read.table(file.names[i], header=F)

# dat<-data[grep("^Supercontig_12.1$", data$V1), ]
# head(dat)
# tail(dat)

#now shift the window by 50 bp
		WINDOWS<-round((dat$V2+50)*0.01)
		D<-dim(dat)[[2]]
		CHRS<-unique(dat$V1)
		for (j in 1:length(CHRS)){
			L<-which(dat$V1==CHRS[j])
			OUT<-unique(WINDOWS[L])
			for (y in 3:D){
			OUT<-cbind(OUT,tapply(dat[L,y], WINDOWS[L], mean))
			}
		M<-length(unique(WINDOWS[L]))
		write.table(cbind(rep(as.character(dat$V1)[L[1]],M),round(OUT,2)), paste0(getwd(),"/",file.names[i],".binned100bp.txt"), sep='\t', col.names=F, row.names=F, quote=F, append=TRUE)
			}
		}
		
# once this is done you can get for each bin the relative frequency of suspicious reads to all
# then yo can caluclate the relative frequency differences using some cutoff for coverage as well 
# at this point there is no standard script so I need to help you!
# but it would be nice to visit the nioo so that is fine
# another loop that processes the binned data txt files

file.names.binned<-dir(path, pattern =".binned100bp.txt")
	par(mfrow=c(6,6))
for (j in 1:length(file.names.binned))
	{
	binned<-read.table(file.names.binned[j], header=F)

	AF_ref<-binned$V4/binned$V3
	AF_evo<-binned$V6/binned$V5

	Call_ref<-round(binned$V3)
	Calt_ref<-round(binned$V4)
	Call_evo<-round(binned$V5)
	Calt_evo<-round(binned$V6)

	LL_evo<-dbinom(Calt_evo, Call_evo, 0.5*(AF_ref+AF_evo))
	LL_ref<-dbinom(Calt_ref, Call_ref, 0.5*(AF_ref+AF_evo))
	totLL<-LL_evo*LL_ref
	Y<-0-log10(totLL)
	binned[which(Y>9),]

	write.table(binned[which(Y>9),], paste0(getwd(),"/",file.names.binned[j],".interesting.stuff.txt"),
	col.names=c("chromos","bin","tot_ref","var_ref","tot_sample","var_sample"),
	row.names=F, sep='\t', quote=F)

	cols<-as.numeric(binned$V1) %% 2 +1
	plot(Y, pch=16, cex=0.6, col=c("black","grey40")[cols],
	xlab="relative place on the genome", 
	ylab="0-log10 p value")
	}

