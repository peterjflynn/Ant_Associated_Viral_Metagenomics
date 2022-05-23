# Ant Associated Viral Metagenomics Workflow
## Dependencies 
[Java](https://www.java.com/en/) <br>
[Perl](https://www.perl.org/) <br>
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) <br>
[Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) <br>
[Samtools](http://samtools.sourceforge.net/) <br>
[Trimmomatic-0.35](http://www.usadellab.org/cms/?page=trimmomatic) <br>
[Spades 3.14.0](https://github.com/ablab/spades/releases) <br>
[BBMap](https://sourceforge.net/projects/bbmap/) <br>
[Diamond](https://github.com/bbuchfink/diamond) <br>
[Blast 2.6.0](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) <br>
[CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) <br>
[CD-Hit v4.8.1](http://weizhong-lab.ucsd.edu/cd-hit/) <br>
[VirSorter2](https://github.com/jiarong/VirSorter2) Note: used [CyVerse](https://de.cyverse.org/) version of VirSorter2 <br>

Phylogenetics Pipeline <br>
[ProtTest v3.4.2](https://github.com/ddarriba/prottest3) <br>
[EMBOSS 6.6.0](http://emboss.sourceforge.net/download/) <br>
[MAFFT v.7.309](https://mafft.cbrc.jp/alignment/software/) <br>
[Geneious 10.2.3](https://www.geneious.com/) <br>
[RAxML v8.2.11](https://cme.h-its.org/exelixis/web/software/raxml/) <br>
[R studio]() <br>
[BaTS v.]() <br>
## Bioinformatics Pipeline
This workflow starts with raw paired-end HiSeq data in demultiplexed FASTQ formated assumed to be located within a folder called raw_seq. For this pipeline I will be illustrating with a single sample (*Cephalotes varians*) 
1. Concatenate the forward and the reverse reads from multiple lanes of sequencing together. 
```sh
cat /raw/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L004_R1.fastq.gz /raw/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L005_R1.fastq.gz /raw/15A_TTGCCTAG-TAAGTGGT-AHWYVLBBXX_L005_R1.fastq.gz > /raw/concat_data/15A_concatenated_R1.fastq.gz

cat /raw/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L004_R2.fastq.gz /raw/15A_TTGCCTAG-TAAGTGGT-AHNFFJBBXX_L005_R2.fastq.gz /raw/15A_TTGCCTAG-TAAGTGGT-AHWYVLBBXX_L005_R2.fastq.gz > /raw/concat_data/15A_concatenated_R2.fastq.gz
```
2. Run fastqc for manual inspection of the quality of the sequences 
```sh
mkdir fastqc_out
fastqc -t 4 raw/data/* -o fastqc_out/
```
###![Happy Christmas](FastQC.png)

3. Trim the contactenated reads to remove adapters and low quality reads using Trimmomatic
```sh
java -jar Trimmomatic-0.35/trimmomatic-0.35.jar PE -phred33 \
raw/concat_data/15A_concatenated_R1.fastq.gz raw/concat_data/15A_concatenated_R2.fastq.gz  \
/data/15_pe1.fq /data/15_unp1.fq /data/15_pe2.fq /data/15_unp2.fq \
ILLUMINACLIP:Trimmomatic-0.35/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

 4. Builds ant species specific genome database in bowtie2. For this particular sample it was Cephalotes varians.
```sh
bowtie2-build /data/D_quadriceps.fna /data/Dquad

```
5. Using bowtie2, map ant genome to paired reads, unmapped reads are carried through to the assembly
```sh
bowtie2 -t -x /data/Dquad -1 15_pe1.fq  -2 15_pe2.fq --un-conc /data/15_conc_unmapped.fastq --al-conc /data/15_conc_mapped.sam

```
6. Use SPAdes assembly program to cross-assemble all unmapped reads from step 5. The single cell option (--sc) worked best for this data, but the --meta option may work better depending on your dataset. 
```sh
SPAdes-3.14.0-Linux/bin/spades.py --sc --pe1-1 /data/15_conc_unmapped.1.fastq --pe1-2 /data/15_conc_unmapped.2.fastq -k 21,33,55,77,99,127  -o 15_spades
```

7. Using removesmalls.pl, only retain contigs 300 bp or larger
```sh
perl /code/removesmalls.pl 300 scaffolds.fasta > 15_scaffolds_300.fasta
```

8. Using Bowtie2, map these contigs back to reads for that sample to assess read coverage of each individual contig.
```sh
bowtie2-build 15_scaffolds_300.fasta 15_scaffolds_300.fasta

bowtie2 -p 12 -x 15_scaffolds_300.fasta -1 15_pe1.fq  -2 15_pe2.fq -S 15_reads.map.sam

samtools faidx 15_scaffolds_300.fasta
```

All data from here on out is in files contigs

9. Run cd-ht on each seperate sample in your dataset. Cd-hit filters the contigs for redundancy, in this case at 95% sequence similarity.
```sh
work_dir="./contigs/"
read_dir="../cdhit_output"

cd "${work_dir}"

for i in *.fasta
do

/cd-hit-v4.8.1-2019-0228/cd-hit -i "${i}" -o "${read_dir}/${i//}_cdhit.fasta" -aS 0.95 -c 0.95 -n 5 -d 0

done

```

10. Concatenate all contigs from together into fasta file.
```sh
cat /data/cdhit_output/*.fasta > /contigs/data/all_contigs_cdhit_decon.fasta
```
11. Sort contigs by length in descending order using BBMap
```sh
/bbmap/sortbyname.sh in=/contigs/data/all_contigs_cdhit.fasta out=/contigs/data/all_contigs_cdhit_decon_sorted.fasta length descending
```

12. Creat single line fasta file from all contigs.
```sh
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' /contigs/data/all_contigs_cdhit_decon_sorted.fasta > /contigs/data/all_contigs_cdhit_decon_sorted_single.fasta
```

**VirSorter2**  
13. Use VirSorter2 version 2.1.0 to identify further viral contigs from your cross-assembled samples. I found that the CyVerse version of VirSorter2 worked better than the command line versions. I used VirSorter2 pre-set parameters for this analysis. Output file from VirSorter2 is: **virsorter_contigs.fa**

Decontaminate all samples with control sample.

14. database with taxonomy for decontamination, these files are very large so I have not included them in this workflow, but you can download them to your server from NCBI. This tutorial is helpful: https://andreirozanski.com/2020/01/03/building-a-diamond-db-using-refseq-protein/

```sh
/diamond/diamond makedb --in /decontamination/23_scaffolds_300.fasta_cdhit.fasta --db /decontamination/decontamination_db --taxonmap /nr/prot.accession2taxid.gz --taxonnodes /nr/nodes.dmp --taxonnames  /nr/names.dmp --threads 20 &
```
15. Create a database of "contaminant" samples using the contigs from the control sample.
```sh
/diamond/diamond makedb --in /data/decontamination/23_scaffolds_300.fasta_cdhit.fasta -d /data/decontamination/decontamination_db
makeblastdb -in /data/decontamination/23_scaffolds_300.fasta_cdhit.fasta -out /data/decontamination/Decon -dbtype nucl -input_type fasta
```

16. Using BLAST, blast all discovered contigs against the control sample. I used BLAST instead of DIAMOND here since it is a nucleotide-nucleotide search and DIAMOND is only functional for protein and translated protein searches.

```sh
blastn -num_threads "40" -db /data/decontamination/Decon -outfmt "6" -max_target_seqs "1" -evalue "1e-5" -max_hsps 1  -out /data/decontamination/contaminated_contigs_300.out -query /contigs/data/all_contigs_cdhit_decon_sorted.fasta &
```
17. **MAKE A CONTAMINATED TEXT FILE FROM BLAST OUTPUT**

18. delete contaminated sequences from contig file
```sh
awk 'BEGIN{while((getline<"contam.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' /contigs/data/all_contigs_cdhit_decon_sorted_single.fasta > /contigs/data/all_contigs_300_decontam_cdhit_single.fasta
```
19. blastx with contig list on RefSeq database. RefSeq Viral Protein database can be downloaded here: [RefSeq Viral Protein Database from NCBI](https://www.ncbi.nlm.nih.gov/protein?term=%28%22Viruses%22%5BOrganism%5D%20AND%20srcdb_refseq%5BPROP%5D%20AND%20viruses%5Bfilter%5D&cmd=DetailsSearch) and converted into a .dmnd file.
#
```sh
/diamond/diamond blastx -p "40" -d /RefSeq_protein_viral_database.dmnd -f "6" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle qtitle qstrand -k "1" --evalue "1e-3" --max-hsps 1 --sensitive -o /contigs/data/RefSeq_blastx_contigs_300.out -q /contigs/data/all_contigs_300_decontam_cdhit_single.fasta &
```
20. extract viral fasta sequences from  blastx refseq results (to make txt file need to do this in text wrangler)
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' /contigs/data/RefSeq_viral_contigs.txt /contigs/data/all_contigs_300_decontam_cdhit_single.fasta > /contigs/data/viral_contigs_refseq_300.fasta
```
21. Delete anything after | in fasta header to make VirSorter output file compatible with RefSeq viral contig output file.
```sh
cut -d'|' -f1 /contigs/data/virsorter_contigs.fa > /contigs/data/virsorter_contigs_1.fa
```
22. Finds viral contigs which were in common between RefSeq and VirSorter searches. 
```sh
awk '/^>/{if (a[$1]>=1){print $1}a[$1]++}' /contigs/data/virsorter_contigs_1.fa /contigs/data/viral_contigs_refseq_300.fasta > /contigs/data/common_viral.txt
```
23. Takes away > from each line 
```sh
perl -pe 's,.*>,,' /contigs/data/common_viral.txt > /contigs/data/common_viral2.txt
```

24. delete same sequences from virsorter contig file
```sh
awk 'BEGIN{while((getline<"common_viral2.txt")>0)l[">"$1]=1}/^>/{f=!l[$1]}f' /contigs/data/virsorter_contigs_1.fa > /contigs/data/virsorter_unique_viruses.fa
```
25. combine viral RefSeq contigs and Virsorter2 contigs
```sh
cat virsorter_unique_viruses.fa viral_contigs_refseq_300.fasta > /contigs/data/Refseq_virsorter_contigs.fasta
```
26. Blastx viral RefSeq and VirSorter contigs on nr
```sh
~/diamond/diamond blastx -p "50" -d ~/nr/nr_diamond.dmnd -f "6" qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames sskingdoms skingdoms sphylums stitle qtitle qstrand -k "1" --evalue "1e-3" --max-hsps 1 --sensitive -o /contigs/data/NR_blastx_contigs_300.out -q /contigs/data/Refseq_virsorter_contigs.fasta &
```
27. Extract viral fasta sequences from  blastx NR results (to make txt file need to do this in text wrangler)
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' /contigs/data/viruses_NR_300.txt /contigs/data/Refseq_virsorter_contigs.fasta > /contigs/data/final_NR_viral_contigs.fasta
```
28. Use CHECKV to look for proviral contamination (on home desktop)
```sh
export CHECKVDB=/checkv-db-v1.0
checkv contamination /contigs/data/final_NR_viral_contigs.fasta  /contigs/data/checkv_output2
```
29. extract viral fasta sequences from  blastx NR results with proviruses and retroviruses and endogenous viruses removed (to make txt file need to do this in text wrangler) and 500 bp
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' /contigs/data/final_viruses.txt /contigs/data/final_NR_viral_contigs.fasta > /contigs/data/final_viruses_aftertaxonomy.fasta
```

## Phylogenetics Pipeline
This portion of the workflow takes only the viral contigs which had most similarity to the viral phylum [Cressdnaviricota](https://talk.ictvonline.org/taxonomy/p/taxonomy-history?taxnode_id=202107372) to illustrate the phylogenetic workflow for each viral clade analyzed. 

1. CRESS viruses using EMBOSS for circular viruses 

2. blastp these 
2. extract Rep sequences for CRESS viruses, only included viruses with complete CRESS genomes (n=170)
```sh
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' Cruciviridae_Rep.txt Cruciviridae_contigs_emboss.fasta > Cruciviridae_Rep.fasta
```
THESE ARE RIGHT
CRESS_viruses_proteins_blastp_emboss_circular1.out
CRESS_rep_proteins_emboss_final.fasta

3. Download all CRESSDNA virus Rep sequences from NCBI 

3. MAFFT

4. Geneious for manual alignment (pictures)

5. Use ProtTest3 to assess best fit model of protein evolution for the CRESS phylogeny.
```sh
java -jar prottest-3.4.2.jar -i viral_alignments/CRESS_refseq_alignment.phy -all-matrices -all-distributions -o viral_alignments/CRESS_refseq.log -threads 30 &
```
6. RAxML 
```sh
raxmlHPC-PTHREADS -n CRESS_refseq -s prottest-3.4.2/viral_alignments/CRESS_refseq_alignment.phy -m PROTGAMMALG -f a -p 194955 -x 12345 -# 500 -T 50 &
```
7. ITOL (pictures) with RAxML output

How to prune tree to just CRESS samples in ITOL 

8. Co-phylo tanglegram in R
```R
setwd(Phylogenetics_data)
library(phytools)
library(ape)
#cophylo CRESS
tr1 <- read.tree("prunedtree_CRESS_JANE_branch.tree")
tr2 <- read.tree("CRESS_nona/CRESS_noboot1.tree")
assoc <- read.csv("CRESS_associations.csv")
obj<-cophylo(tr1,tr2,assoc=assoc, rotate=TRUE, plot=TRUE)
summary(obj)
## add tip labels
#tiplabels.cophylo(pch=21,frame="none",bg="grey",cex=1.5)
#tiplabels.cophylo(pch=21,frame="none",bg="grey",which="right",cex=1.5)
ant_names <- as.factor(obj$assoc[,1])
color_assoc <- rainbow(length(unique(ant_names)))
color_assoc <- c("deepskyblue1", "green", "red3", "darkorange1", "darkviolet", "khaki1", "blue", "red1", "plum1", "goldenrod2", "green4", "black", "goldenrod", "red", "gray4", "deepskyblue2")  
col <- color_assoc[ant_names]

plot(obj, link.col=col, link.lty= "solid", tip.len = 0.02, tip.lty="dashed", part=0.3, pts=FALSE)
plot(obj, ftype= "off", link.col=col, tip.len = 0.1, tip.lty="blank" , part=0.2, pts=FALSE)
tiplabels.cophylo(text=c("Pseudomyrmex", "Camponotus","Odontomachus", "Anochetus", "Ectatomma", "Neoponera", "Azteca", "Cephalotes", "Dolichoderus", "Solenopsis", "Crematogaster", "Labidus", "Eciton", "Gigantiops", "Paraponera", "Daceton"), adj = 0.03, frame= "none", cex=1.2, font=4)
```

9. Procrustes application to cophylogenetic analysis (PACo) in R
```R
setwd(Phylogenetics_data)
library (ape)
library(vegan)
library(paco)
### 1. PACo FUNCTION: adjustemt prior to Procrustes analysis
PACo <- function (H.dist, P.dist, HP.bin)
{ 
  HP.bin <- which(HP.bin > 0, arr.in=TRUE)
  H.PCo <- pcoa(H.dist, correction="cailliez") #Performs PCo of Host distances 
  P.PCo <- pcoa(P.dist, correction="cailliez") #Performs PCo of Parasite distances
  if (is.null(H.PCo$vectors.cor)==TRUE) H.PCo <- H.PCo$vectors else
    H.PCo <- H.PCo$vectors.cor      # returns corrected pcoord 
  if (is.null(P.PCo$vectors.cor)==TRUE) P.PCo <- P.PCo$vectors else
    P.PCo <- P.PCo$vectors.cor
  H.PCo <- H.PCo[HP.bin[,1],]  #adjust Host PCo vectors 
  P.PCo <- P.PCo[HP.bin[,2],]  #adjust Parasite PCo vectors
  list (H.PCo = H.PCo, P.PCo = P.PCo)
}
### 2. DATA INPUT
#2.1 Host and parasite phylogenetic data (should be one of the following):
#2.1.1 Phylogenetic trees:
TreeH <- read.tree(file.choose()) #this function reads Newick trees TRY CRESS_Refseq folder in analysis
TreeP <- read.tree(file.choose()) #for Nexus trees, use read.nexus(file.choose())
#Compute patristic distances:
host.D <- cophenetic (TreeH)
para.D <- cophenetic (TreeP)

#2.2 ## Read HP: host-parasite association matrix
#Hosts in rows, parasites in columns. Taxa names are included in the file and should match those in tree, sequence or distance files. 
HP <- as.matrix(read.table(file.choose(), header=TRUE)) 
#Sort host and parasite taxa in distance matrices to match the HP matrix:
host.D <- host.D[rownames(HP), rownames(HP)]
para.D <- para.D[colnames(HP), colnames(HP)]
# 
### 3. APPLY PACo FUNCTION  
PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
NLinks = sum(HP) #Number of H-P links; needed for further computations
#
#3.1 Plot of host and parasite ordination:
HostX <- HP.proc$X #host ordination matrix
ParY <- HP.proc$Yrot #parasite ordination matrix, scaled and rotated to fit HostX
#Plotting host and parasite ordinations
plot(HostX, asp=1, pch=46) 
points(ParY, pch=1)
arrows(ParY[,1], ParY[,2], HostX[,1], HostX[,2], length=0.12, angle=15, xpd=FALSE)
HostX <- unique(HP.proc$X) 
ParY <- unique(HP.proc$Yrot) #unique() removes duplicated points - convenient for labelling of points below
identify(ParY[,1], ParY[,2], rownames(ParY), offset=0.3, xpd=FALSE, cex=0.8) #interactive labelling
identify(HostX[,1], HostX[,2], rownames(HostX),offset=0.3, xpd=TRUE, cex= 0.8)
#
#3.2 Goodness-of-fit-test
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 10000 #set number of permutations for testing
P.value = 0
seed <-.Random.seed[trunc(runif(1,1,626))]
set.seed(seed)
#set.seed(5) ### use this option to obtain reproducible randomizations
for (n in c(1:N.perm))
{ 
  if (NLinks <= nrow(HP) | NLinks <= ncol(HP)) 	#control statement to avoid all parasites being associated to a single host 
  {	flag2 <- TRUE 
  while (flag2 == TRUE)	{ 
    HP.perm <- t(apply(HP,1,sample))
    if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
  }  
  } else { HP.perm <- t(apply(HP,1,sample))} #permutes each HP row independently
  PACo.perm <- PACo(host.D, para.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
  #write (m2.perm, file = "D:/m2_perm.txt", sep ="\t", append =TRUE) #option to save m2 from each permutation
  if (m2.perm <= m2.obs)
  {P.value = P.value + 1} 
}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")
#
#3.3 Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo) 
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
} 
SQres.jackn <- SQres.jackn**2 #Jackknifed residuals are squared
SQres <- (residuals (HP.proc))**2 # Vector of original square residuals
#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) #standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) #upper 95% confidence interval
#barplot of squared jackknifed residuals
pat.bar <- barplot(phi.mean, names.arg = " ", space = 0.25, col="white", xlab= "Host-parasite links", ylab= "Squared residuals", ylim=c(0, max(phi.UCI)), cex.lab=1.2)
text(pat.bar, par("usr")[3] - 0.001, srt = 330, adj = 0, labels = colnames(SQres.jackn), xpd = TRUE, font = 1, cex=0.6)
arrows(pat.bar, phi.mean, pat.bar, phi.UCI, length= 0.05, angle=90)
abline(a=median(phi.mean), b=0, lty=2, xpd=FALSE) #draws a line across the median residual value
```
10. JANE (in picture format )

11. Ecological Analysis using BaTS
BaTS analysis for diet
```sh
java -jar BaTS_beta.jar single BATS_Picorna_species.tree 1000 22
```
BaTS analysis for bacterial load 
```sh
java -jar BaTS_beta.jar single BATS_Picorna_species.tree 1000 22
```

BaTS analysis for habitat 
```sh
java -jar BaTS_beta.jar single BATS_Picorna_species.tree 1000 22
```
BaTS analysis for nest type 
```sh
java -jar BaTS_beta.jar single BATS_Picorna_species.tree 1000 22
```