RAW intensity files to VCF can be converted with dragena (command-line) or genome studio 2.0 (app). 

# Pipeline for dragena:
```bash
Non-Kazakh individuals were removed based on metadata to leave 224 filtered ones
# Downloaded reference file for hg19 from https://support.illumina.com/downloads/genome-fasta-files.html

array-analysis-cli genotype call \
    --bpm-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1_hg19.bpm \
    --cluster-file /home/user/biostar/gwas/redo_october/making_vcf/clustering_file/GSA-24v2-0_A1_ClusterFile.egt \
    --idat-sample-sheet /home/user/biostar/gwas/redo_october/making_vcf/sample_sheet/SampleSheet_224_aap_cli.csv \
    --output-folder /home/user/biostar/gwas/redo_october/making_vcf/gtc

dragena genotype gtc-to-vcf \
    --bpm-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1_hg19.bpm \
    --genome-fasta-file /home/user/biostar/gwas/redo_october/making_vcf/GRCh37_genome/GRCh37_genome.fa \
    --gtc-sample-sheet /home/user/biostar/gwas/redo_october/making_vcf/sample_sheet/SampleSheet_224_aap_cli_with_gtc.csv \
    --csv-manifest /home/user/biostar/gwas/redo_october/making_vcf/hg19_manifest/GSA-24v2-0_A1.csv \
    --output-folder /home/user/biostar/gwas/redo_october/making_vcf/output_dragena_vcf

bcftools merge ./*.vcf.gz -o merged_output.vcf
# I renamed sentrix names to regular sample names
```

# Quality control of kaz data:
```bash
plink --vcf merged_output.vcf --make-bed --out custom_kaz
awk -F'\t' '!seen[$1]++' GSA-24v2-0_A1_b150_rsids.txt | awk -F'\t' '!seen[$2]++' | awk -F'\t' '$2 !~ /,/' > GSA-dictionary.txt

plink --bfile custom_kaz --update-name GSA-dictionary.txt --make-bed --out custom_kaz1
awk '$2 !~ /^rs/' custom_kaz1.bim | sort -k2,2 > custom_non_rs_SNP.txt
plink --bfile custom_kaz1 --exclude custom_non_rs_SNP.txt --make-bed --out custom_kaz2

# Make FID same as IID
cat custom_kaz2.fam | grep -wf list.txt | cut -d " " -f 1
# I manually changed 1 of 122 and 1 of 123 into 122_1 and 123_1 because ther were 2 of each of them

plink --bfile custom_kaz2 --remove "reproducibility - to_remove.tsv" --make-bed --out custom_kaz3
plink --bfile custom_kaz3 --update-ids "reproducibility - changing_fid.tsv" --make-bed --out custom_kaz4

plink --bfile custom_kaz4 --impute-sex --make-bed --out custom_kaz5
plink --bfile custom_kaz5 --geno 0.02 --make-bed --out custom_kaz6
plink --bfile custom_kaz6 --mind 0.02 --make-bed --out custom_kaz7
plink --bfile custom_kaz7 --maf 0.001 --make-bed --out custom_kaz8

At this stage I removed all non-nucleotide SNPs (indels). But they are already not here. They must have gotten filtered out at some point.
plink --bfile custom_kaz8 --snps-only 'just-acgt' --make-bed --out custom_kaz9

# ROH
plink --bfile custom_kaz9_autosomal \
--homozyg-density 60 \
--homozyg-gap 500 \
--homozyg-window-snp 100 \
--homozyg-window-het 1 \
--out kaz
```

# Working with reference data:
```bash
# Merging with HGDP:
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bed ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.bim ./
cp ../redo_july/working_with_ref_data/hgdp_estonian/HGDP9.fam ./

plink --bfile custom_kaz10 --bmerge HGDP9.bed HGDP9.bim HGDP9.fam --make-bed --out merged1
plink --bfile HGDP9 --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out HGDP10
plink --bfile custom_kaz10 --exclude merged1-merge.missnp --biallelic-only strict --make-bed --out custom_kaz11
plink --bfile custom_kaz11 --bmerge HGDP10.bed HGDP10.bim HGDP10.fam --make-bed --out merged2
plink --bfile merged2 --geno 0.02 --make-bed --out merged3
plink --bfile merged3 --mind 0.02 --make-bed --out merged4
echo -e "rs9267522\nrs11229\nrs75412302\nrs12660719" > duplicates.txt
plink --bfile merged4 --exclude duplicates.txt --make-bed --out merged6

# only 54000 SNPs remaining
wget https://evolbio.ut.ee/turkic/turkic.fam
wget https://evolbio.ut.ee/turkic/turkic.bim
wget https://evolbio.ut.ee/turkic/turkic.bed
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/caucasus/caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/sakha/sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bim
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.fam
wget https://evolbio.ut.ee/jew/jew_paper_data_dbSNP-b131_pos-b37_1KG_strand.bed

plink --bfile caucasus_paper_data_dbSNP-b131_pos-b37_1KG_strand --bmerge jew_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all1
plink --bfile all1 --bmerge sakha_paper_data_dbSNP-b131_pos-b37_1KG_strand --make-bed --out all2
plink --bfile all2 --bmerge turkic --make-bed --out all3

cp ../redo_july/working_with_ref_data/hgdp_estonian/metadata.txt ./

cat metadata.txt | awk '{print $1 "\t" $1}' > ethnic.txt 
plink --bfile all3 --keep ethnic.txt --make-bed --out all4

comm -12 <(awk '{print $2}' merged6.bim | sort) <(awk '{print $2}' all4.bim | sort) > common_snps.txt
cat merged6.bim | awk '{print $2"\t" $1}' > dictionary_chr
cat merged6.bim | cut -f 2,4 > dictionary_pos
plink --bfile all4 --extract common_snps.txt --make-bed --out all5
plink2 --bfile all5 --update-chr dictionary_chr --update-map dictionary_pos --sort-vars --make-pgen --out all6
plink2 --pfile all6 --make-bed --out all7
plink --bfile merged6 --bmerge all7 --make-bed --out all8
# this line above is for making a missnp file

plink --bfile all7 --exclude all8-merge.missnp --make-bed --out all8
plink --bfile merged6 --exclude all8-merge.missnp --make-bed --out merged7

plink --bfile merged7 --bmerge all8 --make-bed --out all9
cp ../redo_july/working_with_ref_data/hgdp_estonian/metadata_kaz.txt ./
cat metadata.txt metadata_kaz.txt > ethnic1.txt

plink --bfile all9 --geno 0.02 --make-bed --out all10
plink --bfile all10 --mind 0.02 --make-bed --out all11
plink --bfile all11 --maf 0.001 --make-bed --out all12
plink --bfile all12 --genome --min 0.2 --out pihat_min0.2
plink --bfile all12 --missing --out missing_report
awk '$10 > 0.2 {print $1, $2, $3, $4}' pihat_min0.2.genome > related_pairs.txt
plink --bfile all12 --snps-only 'just-acgt' --make-bed --out all13
awk '{print $1}' all13.fam | grep -Fwf - ethnic1.txt > ethnic2.txt
cat ethnic_final.tsv | awk '{print $1"\t"$1}' > to_keep_from_metadata.tsv
plink --bfile all13 --keep to_keep_from_metadata.tsv --make-bed --out all14
plink2 --bfile all14 --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnic_final.tsv
cat ethnic_final.tsv | awk '{print $1 "\t" $1 "\t" $2 "\t" $3}' > ethnic2_fst.txt
plink2 --bfile all13 --fst CATPHENO --within ethnic2_fst.txt --double-id --out fst_output
./plot_fst_heatmap.py fst_output.fst.summary sorting_order.tsv

# admixture
plink --bfile all14 --indep-pairwise 1000 150 0.4 --out pruned_data
plink --bfile all14 --extract pruned_data.prune.in --make-bed --out all15
awk '{print $1}' all15.fam | grep -Fwf - ethnic2.txt > ethnic3.txt
for K in 8; do admixture --cv all15.bed -j8 $K | tee log${K}.out; done
ln -s ~/tools/AncestryPainter_v5/AncestryPainter.pl ./
cat ethnic3.txt | awk '{print $2"\t" $1}'  > ethnic3.ind 
perl AncestryPainter.pl -i ethnic3.ind -q ./all15.8.Q -t Kazakh -o Kazakh -f png
```

# Getting docs and graphs for publication:
```bash
# autosomal, mitochnodrial, y-chr
awk '{ if ($1 >= 1 && $1 <= 22) print $2 }' custom_kaz9.bim > snp_1_22.txt
plink --bfile custom_kaz9 --extract snp_1_22.txt --make-bed --out custom_kaz9_autosomal

awk '{ if ($1 == 26) print $2 }' custom_kaz9.bim > snp_mitoch.txt
plink --bfile custom_kaz9 --extract snp_mitoch.txt --make-bed --out custom_kaz9_mitoch

awk '{ if ($1 == 24) print $2 }' custom_kaz9.bim > snp_y.txt
plink --bfile custom_kaz9 --extract snp_y.txt --make-bed --out custom_kaz9_y_chr

#mafs
plink2 --bfile custom_kaz9  --freq --out maf_custom_kaz9
plink2 --bfile custom_kaz9_autosomal  --freq --out maf_custom_kaz9_autosomal
plink2 --bfile custom_kaz9_mitoch  --freq --out maf_custom_kaz9_mitoch
plink2 --bfile custom_kaz9_y_chr --freq --out maf_custom_kaz9_y_chr
sed -i 's/ALT_FREQS/kaz_alt_frq/g' maf_custom_kaz9_autosomal.afreq  maf_custom_kaz9_mitoch.afreq maf_custom_kaz9_y_chr.afreq maf_custom_kaz9.afreq
sed -i 's/OBS_CT/kaz_allele_count/g' maf_custom_kaz9_autosomal.afreq  maf_custom_kaz9_mitoch.afreq maf_custom_kaz9_y_chr.afreq maf_custom_kaz9.afreq

#making vcfs
plink --bfile custom_kaz9_autosomal --recode vcf --out custom_kaz9_autosomal
plink --bfile custom_kaz9_mitoch --recode vcf --out custom_kaz9_mitoch
plink --bfile custom_kaz9_y_chr --recode vcf --out custom_kaz9_y_chr
plink --bfile custom_kaz9 --recode vcf --out custom_kaz9_all

# making extended vcfs
bcftools view -h custom_kaz9_autosomal.vcf > kaz_a1.vcf && bcftools view -H custom_kaz9_autosomal.vcf > kaz_a2.vcf && tail -n +2 maf_custom_kaz9_autosomal.afreq | cut -f 6,7 > added_info.txt && head -n 1 maf_custom_kaz9_autosomal.afreq | cut -f 6,7 > added_header.txt && (sed '$d' kaz_a1.vcf; paste <(tail -n 1 kaz_a1.vcf) added_header.txt) > kaz_a4.vcf && paste kaz_a2.vcf added_info.txt > kaz_a3.vcf && cat kaz_a4.vcf kaz_a3.vcf > kaz_a5.vcf && mv kaz_a5.vcf ./autosomal_ext_for_annovar.vcf

# IBD between kazakhs only and between different ethnicities
plink --bfile custom_kaz9_autosomal --genome --out ibd_kaz
plink --bfile all14 --genome --out ibd_all
./average_ibd.py ibd_all.genome ethnic_final.tsv
cat ibd_all_average.txt | grep Kazakh | sort -gk 3,3 > average_ibd_kazakh.tsv

# ROH
plink --bfile custom_kaz9_autosomal --exclude merged-merge.missnp --bmerge hapmap3.hg19.bed hapmap3.hg19.bim hapmap3.hg19.fam --make-bed --out merged
plink --bfile custom_kaz9_autosomal --bmerge all3.bed all3.bim all3.fam --make-bed --out merged
plink --bfile custom_kaz9_autosomal --bmerge hgdp.hg19.bed hgdp.hg19.bim hgdp.hg19.fam --make-bed --out merged
```




Problem - too few SNPs with all my referent datasets. Potentially need WGS referent datasets to retain 300000 SNPs

Solution: found WGS HGDP dataset in PLINK format in PLINK resources website. Workflow:
```bash
# 1) un-archive files and rename for conversion (I will use phased dataset)
mv hgdp.psam hgdp_statphase.psam
unzstd *hgdp_statphase*
# 2) turn into PLINK binary files (split multi-allelic variants)
plink2 --pfile hgdp_statphase --make-bed --out hgdp_hg38 --allow-extra-chr --max-alleles 2
# 3) liftover my data to hg38
plink2 --bfile custom_kaz9_autosomal --recode ped --out kaz_hg19
awk '{print "chr"$1,"\t",$4,"\t",$4+1,"\t",$2}' kaz_hg19.map > kaz_hg19.bed
~/tools/liftover/liftOver kaz_hg19.bed ~/tools/liftover/hg19ToHg38.over.chain.gz kaz_hg38.bed unmap
awk '!/^#/ {print $4}' unmap > unmapped.txt
plink --file kaz_hg19 --exclude unmapped.txt --recode --out kaz_hg38
awk '{print $1,"\t",$4,"\t",0,"\t",$2}' kaz_hg38.bed > kaz_hg38.map
perl -p -i -e 's/chr//g' kaz_hg38.map
plink --file kaz_hg38 --make-bed --allow-extra-chr --out kaz_hg38
# 4) merge with HGDP dataset
cat kaz_hg38.bim | cut -f 2 > to_keep.tsv
plink2 --bfile hgdp_hg38 --extract to_keep.tsv --make-bed --out hgdp_for_merge
plink --bfile kaz_hg38 --bmerge hgdp_for_merge  --threads 32 --allow-extra-chr --make-bed --out merged
plink --bfile hgdp_for_merge --flip merged-merge.missnp --make-bed --out hgdp_for_merge2
plink --bfile kaz_hg38 --bmerge hgdp_for_merge2  --threads 32 --allow-extra-chr --make-bed --out merged
plink --bfile hgdp_for_merge2 --exclude merged-merge.missnp --make-bed --out hgdp_for_merge3
plink --bfile kaz_hg38 --bmerge hgdp_for_merge3  --threads 32 --allow-extra-chr --make-bed --out merged

# 5) QC of merged data
plink --bfile merged --allow-extra-chr --chr 1-22 --make-bed --out merged1
plink --bfile merged1 --geno 0.02 --make-bed --out merged2
plink --bfile merged2 --mind 0.02 --make-bed --out merged3
plink --bfile merged3 --maf 0.001 --make-bed --out merged4
awk '{$6 = -9; print}' merged4.fam > merged4_modified.fam && mv merged4_modified.fam merged4.fam
# modifying FID to be IID
cat merged4.fam | awk '{print $2"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' > merged5.fam
cp merged4.bim ./merged5.bim
cp merged4.bed ./merged5.bed
plink --bfile merged5 --keep to_keep_samples.tsv --make-bed --out merged6
plink2 --bfile merged6 --pca 10 --out all_pca
python plot_eigenvec.py all_pca.eigenvec ethnic.tsv
#PCA looks as expected

# 6) run ROH calculations for each individual
plink --bfile merged6 \
--homozyg \
--homozyg-snp 50 \
--homozyg-kb 1500 \
--homozyg-window-snp 100 \
--homozyg-density 60 \
--homozyg-gap 500 \
--homozyg-window-het 1 \
--homozyg-window-missing 1 \
--out kaz_hgdp_roh

# 7) calculate F(ROH) for each individuals, caclulate number of segments of different length for each individual, and calculate average parameters for all inidividuals in each population
# in python script FROH_script.ipynb
```





LSBL test for statistically significant differences between populations:
# 1) using HGDP WGS data get 3 datasets (separate) with Kazakh (224) + East Asian (152) + European (121)
```bash
cat ethnic.tsv | grep -e French -e Basque -e Sardinian -e Russian -e Orcadian -e Tuscan -e Bergamo | awk '{print $1"\t"$1}' > samples_EUR.txt
cat ethnic.tsv | grep -e Han -e Japanese -e Dai -e Hezhen -e Miao -e Naxi -e Oroqen -e She -e Tujia -e Tu -we Xibo -e Yi -e Cambodian -e Lahu -e Daur -e Mongolian | awk '{print $1"\t"$1}' > samples_EAS.txt
cat ethnic.tsv | grep Kazakh | awk '{print $1"\t"$1}' > samples_KAZ.txt
plink --bfile merged6 --keep samples_EUR.txt --make-bed --out only_EUR
plink --bfile merged6 --keep samples_EAS.txt --make-bed --out only_EAS
plink --bfile merged6 --keep samples_KAZ.txt --make-bed --out only_KAZ

# 2) get MAFs for all SNPs for each of these 3 datasets
plink2 --bfile only_EUR  --freq --out only_EUR
plink2 --bfile only_EAS  --freq --out only_EAS
plink2 --bfile only_KAZ  --freq --out only_KAZ

# 3) merge 3 MAF tables together
join -1 1 -2 1 \
  <(awk 'NR>1 {print $2, $6}' only_EUR.afreq | sort -k1,1) \
  <(awk 'NR>1 {print $2, $6}' only_EAS.afreq | sort -k1,1) | \
join -1 1 -2 1 - \
  <(awk 'NR>1 {print $2, $6}' only_KAZ.afreq | sort -k1,1) | \
awk 'BEGIN{OFS="\t"; print "rsID\tMAF_EUR\tMAF_EAS\tMAF_KAZ"} {print $1, $2, $3, $4}' > cat_3pops_mafs.txt

# 4) Get pairwise Fst values between KAZ, EUR, EAS
Via PLINK (not using MAFs calculated earlier):
cat only_EUR.fam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""EUR"}' > only_EUR1.fam
cat only_EAS.fam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""EAS"}' > only_EAS1.fam
cat only_KAZ.fam | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t""KAZ"}' > only_KAZ1.fam

echo -e "FID\tIID\tPOP" > pheno_file.txt
cat only_EUR1.fam | cut -f 1,2,6 >> pheno_file.txt
cat only_EAS1.fam | cut -f 1,2,6 >> pheno_file.txt
cat only_KAZ1.fam | cut -f 1,2,6 >> pheno_file.txt
plink --bfile only_KAZ --bmerge only_EUR.bed only_EUR.bim only_EUR1.fam  --make-bed --out temp1
plink --bfile temp1 --bmerge only_EAS.bed only_EAS.bim only_EAS1.fam  --make-bed --out temp2

plink2 --bfile temp2 \
  --pheno pheno_file.txt \
  --pheno-name POP \
  --fst POP report-variants \
  --out pairwise_fst

join -1 1 -2 1 \
  <(awk 'NR>1 {print $3, $5}' pairwise_fst.EAS.EUR.fst.var | sort -k1,1) \
  <(awk 'NR>1 {print $3, $5}' pairwise_fst.EAS.KAZ.fst.var | sort -k1,1) | \
join -1 1 -2 1 - \
  <(awk 'NR>1 {print $3, $5}' pairwise_fst.EUR.KAZ.fst.var | sort -k1,1) | \
awk 'BEGIN{OFS="\t"; print "rsID\tFST_EAS_EUR\tFST_EAS_KAZ\tFST_EUR_KAZ"} {print $1, $2, $3, $4}' > FSTs.txt
```


After obtaining the test results, I want to find new SNPs to discuss in the paper:
First I made a list of SNPs with empirical P-value < 0.01 (4969 entries).
Among those we will search for those with ANNOVAR charactertistics:

Among those SNPs I will find here let's search if they are present in my 5000 list:
```bash
cat final_annovared_extended.tsv | grep -e "LCT" -e "ExonicFunc.knownGene" -e "rs4988235" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > lactose.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -e "ALDH2" -e "ADH" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > alcohol.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -w -e "IFNL3" -e "NUDT15" -e "SLCO1B1" -e "TPMT" -e "UGT1A1" -e "CFTR" -e "CYP2B6" -e "CYP2C19" -e "CYP2C9" -e "CYP2D6" -e "CYP3A5" -e "CYP4F2" -e "DPYD" -e "VKORC1" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > pharm_idda.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | grep -w -f druggable_genome.tsv  | awk '$14 == "ExonicFunc.knownGene" || $70 == "D" && $76 == "D"' | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > druggable_mafs.tsv
cat final_annovared_extended.tsv | grep -e "exonic" -e "ExonicFunc.knownGene" | $70 == "D" && $76 == "D"' | grep -w -e "thyroid" -e "ExonicFunc.knownGene" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > thyroid.tsv
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $70 == "D" || $76 == "D"' final_annovared_extended.tsv | grep "thyroid" | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > thyroid.tsv
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $70 == "D" && $76 == "D"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > deleterious.tsv
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $333 == "drug_response"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > clinvar_drug_response.tsv
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $333 == "pathogenic"|| $333 == "risk_factor"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > clinvar_pathogenic_riskfactor.tsv

# merging them all together:
cat lactose.tsv alcohol.tsv pharm_idda.tsv druggable_mafs.tsv thyroid.tsv deleterious.tsv clinvar_drug_response.tsv clinvar_pathogenic_riskfactor.tsv > snps_of_interest.tsv

# find
cat snps_of_interest.tsv | grep -f final_list.txt

# 5) Calculate LSBL = FST(KAZ,EUR)+FST(KAZ,EAS)-FST(EUR,EAS)/2
# 6) Rank LSBLs for all SNPs in the genome and calculate where your SNPs of interest locate in it; provide percentile (empirical P-value) - genome-wide empirical distribution, using the formula PE (x)=(number of loci>x)/(total number loci)
# In Python script (notebook)
```

Calculating number of deleterious variants per individual:
```bash
awk -F'\t' '$14 == "ExonicFunc.knownGene" || $70 == "D" && $76 == "D"' final_annovared_extended.tsv | cut -f 7,9,12,14,17,19,27,303,306,307,456,687 > deleterious.tsv

cat deleterious.tsv | cut -f 12 | awk 'NR > 1 {sum += $1 * 224} END {print sum / 224}' 
# There was an average of 81.2555 likely deleterious nonsynonymous single nucleotide variants (SNVs) (those which are classified as deleterious by both SIFT and PolyPhen databases) per individual

awk -F'\t' '$70 == "D" && $76 == "D" {for (i=463; i<=686; i++) if ($i == "1/1") count[i]++} END {for (i=463; i<=686; i++) print "Individual " i-462 ": " count[i]}' final_annovared_extended.tsv 
awk -F'\t' '$70 == "D" && $76 == "D" {for (i=463; i<=686; i++) if ($i == "1/1") count[i]++} END {sum=0; for (i=463; i<=686; i++) sum+=count[i]; print "Average: ", sum/(686-463+1)}' final_annovared_extended.tsv
# An average of 19 variants per individual were homozygotes, therefore representing the presence of a considerable amount of potentially deleterious variation in the Kazakh population. 
```

```bash
# first - turning EAS and EUR to vcf files:
plink --bfile only_EAS --recode vcf --out only_EAS
plink --bfile only_EUR --recode vcf --out only_EUR
plink --bfile only_KAZ --recode vcf --out only_KAZ
# for those rows where only_EAS.vcf column 3 is present in column 23 of final_annovared_extended.tsv, add columns 70 and 76 from final_annovared_extended.tsv as the last 2 columns in the annotated_only_EAS.vcf
cat final_annovared_extended.tsv | cut -f 23,70,76 > rs_sift_polyphen.txt

awk -F'\t' 'FNR==NR {map[$1]=$2"\t"$3; next} 
            /^#/ {print; next} 
            ($3 in map) {print $0"\t"map[$3]}' OFS='\t' rs_sift_polyphen.txt only_EAS.vcf > annotated_only_EAS.vcf

awk -F'\t' 'FNR==NR {map[$1]=$2"\t"$3; next} 
            /^#/ {print; next} 
            ($3 in map) {print $0"\t"map[$3]}' OFS='\t' rs_sift_polyphen.txt only_EUR.vcf > annotated_only_EUR.vcf

awk -F'\t' 'FNR==NR {map[$1]=$2"\t"$3; next} 
            /^#/ {print; next} 
            ($3 in map) {print $0"\t"map[$3]}' OFS='\t' rs_sift_polyphen.txt only_KAZ.vcf > annotated_only_KAZ.vcf


awk -F'\t' '$162 == "D" && $163 == "D" {for (i=10; i<=161; i++) if ($i == "1/1") count[i]++} END {sum=0; for (i=10; i<=161; i++) sum+=count[i]; print "Average EAS: ", sum/(161-10+1)}' annotated_only_EAS.vcf

awk -F'\t' '$131 == "D" && $132 == "D" {for (i=10; i<=130; i++) if ($i == "1/1") count[i]++} END {sum=0; for (i=10; i<=130; i++) sum+=count[i]; print "Average EUR: ", sum/(130-10+1)}' annotated_only_EUR.vcf

awk -F'\t' '$234 == "D" && $235 == "D" {for (i=10; i<=233; i++) if ($i == "1/1") count[i]++} END {sum=0; for (i=10; i<=233; i++) sum+=count[i]; print "Average KAZ: ", sum/(233-10+1)}' annotated_only_KAZ.vcf

plink --bfile only_EAS --freq --out only_EAS
plink --bfile only_EUR --freq --out only_EUR
plink --bfile only_KAZ --freq --out only_KAZ
cat only_KAZ.frq | awk '{print $5}' 
cat only_EUR.frq | awk '{print $5}'
cat only_EAS.frq | awk '{print $5}'

awk -F'\t' '$162 == "D" && $163 == "D" {for (i=10; i<=161; i++) if ($i == "1/1") count[i]++} END {for (i=10; i<=161; i++) print "Individual " i-9 ": " count[i]}' annotated_only_EAS.vcf

awk -F'\t' '$131 == "D" && $132 == "D" {for (i=10; i<=130; i++) if ($i == "1/1") count[i]++} END {for (i=10; i<=130; i++) print "Individual " i-9 ": " count[i]}' annotated_only_EUR.vcf

awk -F'\t' '$234 == "D" && $235 == "D" {for (i=10; i<=233; i++) if ($i == "1/1") count[i]++} END {for (i=10; i<=233; i++) print "Individual " i-9 ": " count[i]}' annotated_only_KAZ.vcf

#Average EAS:  12.0855
#Average EUR:  10.7851
#Average KAZ:  11.2857


# using not only 1/1 but also 0/1,1/0 but it shows 1297 for all. why?
awk -F'\t' '$162 == "D" && $163 == "D" {for (i=10; i<=161; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=161; i++) sum+=count[i]; print "Average EAS: ", sum/(161-10+1)}' annotated_only_EAS.vcf

awk -F'\t' '$131 == "D" && $132 == "D" {for (i=10; i<=130; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=130; i++) sum+=count[i]; print "Average EUR: ", sum/(130-10+1)}' annotated_only_EUR.vcf

awk -F'\t' '$234 == "D" && $235 == "D" {for (i=10; i<=233; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=233; i++) sum+=count[i]; print "Average KAZ: ", sum/(233-10+1)}' annotated_only_KAZ.vcf


# not dividing
awk -F'\t' '$162 == "D" && $163 == "D" {for (i=10; i<=161; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=161; i++) sum+=count[i]; print " EAS: ", sum}' annotated_only_EAS.vcf

awk -F'\t' '$131 == "D" && $132 == "D" {for (i=10; i<=130; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=130; i++) sum+=count[i]; print " EUR: ", sum}' annotated_only_EUR.vcf

awk -F'\t' '$234 == "D" && $235 == "D" {for (i=10; i<=233; i++) if ($i == "1/1" || $i = "0/1" || $i == "1/0") count[i]++} END {sum=0; for (i=10; i<=233; i++) sum+=count[i]; print " KAZ: ", sum}' annotated_only_KAZ.vcf


# only 0/0
awk -F'\t' '$162 == "D" && $163 == "D" {for (i=10; i<=161; i++) if ($i == "0/0") count[i]++} END {sum=0; for (i=10; i<=161; i++) sum+=count[i]; print "Average EAS: ", sum/(161-10+1)}' annotated_only_EAS.vcf
awk -F'\t' '$131 == "D" && $132 == "D" {for (i=10; i<=130; i++) if ($i == "0/0") count[i]++} END {sum=0; for (i=10; i<=130; i++) sum+=count[i]; print "Average EUR: ", sum/(130-10+1)}' annotated_only_EUR.vcf
awk -F'\t' '$234 == "D" && $235 == "D" {for (i=10; i<=233; i++) if ($i == "0/0") count[i]++} END {sum=0; for (i=10; i<=233; i++) sum+=count[i]; print "Average KAZ: ", sum/(233-10+1)}' annotated_only_KAZ.vcf
