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
[Blast v.]()
[CheckV](https://bitbucket.org/berkeleylab/checkv/src/master/) <br>
[CD-Hit v4.8.1](http://weizhong-lab.ucsd.edu/cd-hit/) <br>
[VirSorter2](https://github.com/jiarong/VirSorter2) Note: used [CyVerse](https://de.cyverse.org/) version of VirSorter2 <br>

[ProtTest v3.4.2](https://github.com/ddarriba/prottest3) <br>
[EMBOSS 6.6.0](http://emboss.sourceforge.net/download/) <br>
[MAFFT v.7.309](https://mafft.cbrc.jp/alignment/software/) <br>
[Geneious 10.2.3](https://www.geneious.com/) <br>
[RAxML v8.2.11](https://cme.h-its.org/exelixis/web/software/raxml/) <br>

## Bioinformatics Pipeline
This workflow starts with raw paired-end HiSeq data in demultiplexed FASTQ formated assumed to be located within a folder called raw_seq. For this pipeline I will be illustrating with a single sample (Cephalotes varians) 
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

VirSorter2  
13. Use VirSorter2 version 2.1.0 to identify further viral contigs from your cross-assembled samples. I found that the CyVerse version of VirSorter2 worked better than the command line versions. I used VirSorter2 pre-set parameters for this analysis. Output file from VirSorter2 is: virsorter_contigs.fa

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
17. MAKE A CONTAMINATED TEXT FILE FROM BLAST OUTPUT

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
This portion of the workflow takes only the viral contigs which had most similarity to the viral phylum CRESSDnaviricota (include link here) to illustrate the phylogenetic workflow for each viral clade analyzed. 


3. Use ProtTest3 to assess best fit model of protein evolution for the CRESS phylogeny.
```sh
java -jar prottest-3.4.2.jar -i viral_alignments/CRESS_refseq_alignment.phy -all-matrices -all-distributions -o viral_alignments/CRESS_refseq.log -threads 30 &
```
4. RAxML 
```sh
raxmlHPC-PTHREADS -n CRESS_refseq -s prottest-3.4.2/viral_alignments/CRESS_refseq_alignment.phy -m PROTGAMMALG -f a -p 194955 -x 12345 -# 500 -T 50 &
```

### Cruciviridae Rep ####
awk -F'>' 'NR==FNR{ids[$0]; next} NF>1{f=($2 in ids)} f' Cruciviridae_Rep.txt Cruciviridae_contigs_emboss.fasta > Cruciviridae_Rep.fasta

MAFFT
Geneious for manual alignment 

Co-phylo in R and PACO


JANE (in picture format )


BaTS analysis
```sh
java -jar BaTS_beta.jar single BATS_Picorna_species.tree 1000 22
```