
--------------shell script for merging samples before running SNPRelate.R-----------------

#### REMEMBER TO ADD YOUR SAMPLES INFO TO THE FILE Population_Sample_Info.txt

module load vcftools
module load samtools/0.1.19
module load R/3.2.2 

for i in `cat IDs`
do
cd ${i}
bgzip ${i}_raw.vcf
tabix -p vcf ${i}_raw.vcf.gz
cp ${i}_raw.vcf.gz* ../
cd ..
done

vcf-merge $(ls -1 *.vcf.gz | perl -pe 's/\n/ /g') > merge.vcf


for i in `cat IDs`
do
rm ${i}_raw.vcf.gz*
done

bgzip merge.vcf
tabix -p vcf merge.vcf.gz
#vcf-merge /scratch/EXOME_DATA/Plinkseq/CEU_JPT_CHB_YRISMIRH.vcf.gz merge.vcf.gz > temp.vcf

wget ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20100804/ALL.2of4intersection.20100804.genotypes.vcf.gz
tabix -p vcf ALL.2of4intersection.20100804.genotypes.vcf.gz
vcf-merge ALL.2of4intersection.20100804.genotypes.vcf.gz merge.vcf.gz > temp.vcf

R CMD SNPRelate.R

################## SNPRelate.R script below   #############


### Install the Bioconductor R packages gdsfmt and SNPRelate ###

source("http://bioconductor.org/biocLite.R")
biocLite("gdsfmt")
biocLite("SNPRelate")


#Load the pakcages gdsfmt and SNPRelate into your R session #

library(gdsfmt)
library(SNPRelate)


## Set your working directory to where the VCF files are located - it takes also plink files
##setwd("/Users/tathornt/SISG2015_Association_Mapping/Data/")


## Convert PLINK FILES INTO GDS FORMAT FOR R
#snpgdsBED2GDS(bedfile, famfile, bimfile, "IPCgeno.gds")


##Convert merged VCF file to GDS format 

vcf.fn <- "merge.vcf.gz"
snpgdsVCF2GDS(vcf.fn, "test3.gds", method="biallelic.only")


## Open the GDS file

genofile <- snpgdsOpen("test3.gds")
head(genofile)

head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

## LD-based SNP pruning
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)


#-------------------------IBS---------------------------

RT <- c ("Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Unrelated","Twins_duplicated","Twins_duplicated","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Parent_offspring","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","Other_related","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS","FS")
Z0 <- c( 0.9498, 0.949, 0.9486, 0.9483, 0.9483, 0.9464, 0.9456, 0.9364, 0.9338, 0.9324, 0.9308, 0.9289, 0.9279, 0.9264, 0.9257, 0.9249, 0.9211, 0.9207, 0.9167, 0.9163, 0.9147, 0.91, 0.907, 0.8928, 0.8899, 0.8898, 0.8881, 0.8839, 0.8836, 0.882, 0.874, 0.8722, 0.8703, 0.8689, 0.8682, 0.8666, 0.8659, 0.8637, 0.8583, 0.857, 0.854, 0.8531, 0.8529, 0.8452, 0.8421, 0.8415, 0.8415, 0.8405, 0.8396, 0.8394, 0.8349, 0.8334, 0.8333, 0.8216, 0.8199, 0.8184, 0.8169, 0.8153, 0.8126, 0.8081, 0.8072, 0.8009, 0.8004, 0.8, 0, 0, 0.0024, 0.0096, 0, 0.0024, 0.0096, 0.0072, 0.0048, 0.0071, 0.0027, 0.0072, 0.0048, 0.0024, 0.0095, 0.0024, 0.0048, 0.0024, 0.0024, 0, 0.012, 0.0024, 0.0024, 0, 0.0072, 0.012, 0.0024, 0.0095, 0.0168, 0.0048, 0.0048, 0, 0, 0.0072, 0.0048, 0.0048, 0.0024, 0.0055, 0.0024, 0.0143, 0.0027, 0.0048, 0.0028, 0.0096, 0, 0.0049, 0.0096, 0.0024, 0.0024, 0.0074, 0.0096, 0.0194, 0.0097, 0, 0.0123, 0.0048, 0.0072, 0.0049, 0.0074, 0.0072, 0.0098, 0.0024, 0.0049, 0.0028, 0.0024, 0.0072, 0.0024, 0.0024, 0.0024, 0.0121, 0.0072, 0.0122, 0.0049, 0.0098, 0.0049, 0.0049, 0.0024, 0.0048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0097, 0.0024, 0.0073, 0.0049, 0.012, 0.0073, 0.0072, 0.0073, 0.0073, 0.0082, 0.0097, 0.0096, 0.0121, 0.7914, 0.7902, 0.7861, 0.7812, 0.779, 0.7753, 0.7712, 0.7707, 0.77, 0.7686, 0.7631, 0.7628, 0.7592, 0.7537, 0.753, 0.7527, 0.7524, 0.7496, 0.7476, 0.7459, 0.7402, 0.7382, 0.7262, 0.7089, 0.7082, 0.6892, 0.6876, 0.6569, 0.6473, 0.6083, 0.4522, 0.4055, 0.4058, 0.4077, 0.4122, 0.4186, 0.4189, 0.421, 0.426, 0.4277, 0.432, 0.4336, 0.4342, 0.4389, 0.4533, 0.4535, 0.4628, 0.4654, 0.466, 0.4697, 0.4791, 0.4797, 0.4816, 0.4837, 0.484, 0.4865, 0.4878, 0.4888, 0.4892, 0.4904, 0.4931, 0.4937, 0.4948, 0.4996, 0.5012, 0.503, 0.5063, 0.5077, 0.5102, 0.5135, 0.5215, 0.5234, 0.5236, 0.5256, 0.527, 0.5269, 0.529, 0.5299, 0.5321, 0.5339, 0.5344, 0.5797, 0.5397, 0.5455, 0.5472, 0.5543, 0.5566, 0.561, 0.561, 0.5655, 0.566, 0.5669, 0.5668, 0.5712, 0.5742, 0.5809, 0.5837, 0.584, 0.5851, 0.5864, 0.5867, 0.5876, 0.5879, 0.5882, 0.5885, 0.5946, 0.5953)
Z1 <- c( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.0211, 0, 0.0339, 0.0135, 0, 0, 0.0543, 0, 0.0428, 0.0334, 0.0573, 0.0822, 0.1072, 0.1047, 0.1102, 0.0647, 0.1161, 0.1164, 0.1135, 0.1021, 0.0983, 0.1297, 0.1311, 0.1318, 0.1334, 0.1143, 0.13, 0.1296, 0.1247, 0.146, 0.1469, 0.1471, 0.1137, 0.1464, 0.1521, 0.1585, 0.1595, 0.1604, 0.1606, 0.1651, 0.1385, 0.1576, 0.1784, 0.1561, 0.1816, 0.1364, 0.1847, 0.1874, 0.1919, 0.1928, 0.166, 0.1783, 0.2, 0.0157, 0.0435, 0.8857, 0.8721, 0.8951, 0.8905, 0.8771, 0.8841, 0.8915, 0.8882, 0.8998, 0.8935, 0.8994, 0.9049, 0.8944, 0.91, 0.906, 0.9115, 0.9123, 0.9177, 0.8972, 0.9204, 0.9256, 0.9313, 0.9191, 0.911, 0.9317, 0.9179, 0.9062, 0.9326, 0.9366, 0.9482, 0.9494, 0.9364, 0.942, 0.9424, 0.9477, 0.9419, 0.9496, 0.9299, 0.9531, 0.9493, 0.9561, 0.9448, 0.9651, 0.9556, 0.9482, 0.9627, 0.966, 0.9569, 0.9528, 0.9336, 0.9537, 0.9742, 0.9523, 0.9708, 0.9664, 0.9719, 0.9679, 0.9692, 0.9642, 0.9792, 0.9751, 0.9802, 0.982, 0.9729, 0.9835, 0.9839, 0.9851, 0.9676, 0.9794, 0.9698, 0.9851, 0.9758, 0.9856, 0.986, 0.9919, 0.9883, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0.9822, 0.9976, 0.9884, 0.9939, 0.9803, 0.9915, 0.9928, 0.9927, 0.9927, 0.9918, 0.9903, 0.9904, 0.9879, 0.2086, 0.2098, 0.2139, 0.1963, 0.221, 0.2247, 0.2172, 0.2293, 0.23, 0.2314, 0.2369, 0.1837, 0.2408, 0.2463, 0.247, 0.2473, 0.2476, 0.2504, 0.2524, 0.2541, 0.2598, 0.2269, 0.2738, 0.2911, 0.2918, 0.3108, 0.309, 0.3431, 0.2963, 0.3917, 0.4941, 0.5945, 0.5942, 0.5923, 0.5878, 0.5814, 0.5811, 0.579, 0.574, 0.5723, 0.568, 0.5664, 0.5658, 0.5611, 0.5467, 0.5465, 0.5372, 0.5346, 0.534, 0.5303, 0.5209, 0.5203, 0.5184, 0.5163, 0.516, 0.5135, 0.5122, 0.5112, 0.5108, 0.5096, 0.5069, 0.5063, 0.5052, 0.5004, 0.4988, 0.497, 0.4937, 0.4923, 0.4898, 0.4865, 0.4785, 0.4766, 0.4764, 0.4744, 0.473, 0.4731, 0.471, 0.4701, 0.4679, 0.4661, 0.4656, 0.3794, 0.4603, 0.4545, 0.4528, 0.4457, 0.4434, 0.439, 0.439, 0.4345, 0.434, 0.4331, 0.4332, 0.4288, 0.4258, 0.4191, 0.4163, 0.416, 0.4149, 0.4136, 0.4133, 0.4124, 0.4121, 0.4118, 0.4115, 0.4054, 0.4047)
 
datafile_contaminated = data.frame (RT=RT,Z0=Z0,Z1=Z1)
colfunc <- colorRampPalette(c("green", "red", "yellow"))


sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
ibd <- snpgdsIBDMoM(genofile, sample.id=sample.id, snp.id=snpset.id,
maf=0.05, missing.rate=0.05, num.thread=2)
    
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

jpeg('IBS.jpg')
plot(Z0, Z1,main="IBS PLOT", col = colfunc(length(RT))[rank(RT)],cex=0.75, pch=16)
## Add grid behind the plot
axis(1, tck=1, col.ticks="light gray")
axis(1, tck=-0.015, col.ticks="black")
axis(2, tck=1, col.ticks="light gray", lwd.ticks="1")
axis(2, tck=-0.015)

legend(1,1, xjust=1, yjust=1, legend=levels(datafile_contaminated$RT), pch=16, col=colfunc(5))
##### Plot MY DATA

points(ibd.coeff$k0, ibd.coeff$k1, pch=17,cex=1)
dev.off()

#-----------------PCA/MDS---------------------------------------

#### REMEMBER TO ADD YOUR SAMPLES IDS TO THE FILE Population_Sample_Info.txt

##Convert merged VCF file to GDS format 

vcf.fn <- "temp.vcf"
snpgdsVCF2GDS(vcf.fn, "test2.gds", method="biallelic.only")


## Open the GDS file

genofile <- snpgdsOpen("test2.gds")
head(genofile)

head(read.gdsn(index.gdsn(genofile, "sample.id")))
head(read.gdsn(index.gdsn(genofile, "snp.id")))

## LD-based SNP pruning
set.seed(1000)

# Try different LD thresholds for sensitivity analysis
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)
snpset.id <- unlist(snpset)

######


pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)

# make a data.frame
tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)


## NOW INCORPORATE THE POPULATION MEMBERSHIPS IN THE PCA PLOT
POPINFO=read.table(file="Population_Sample_Info.txt",header=TRUE)

## Obtain the number of sample individuals from each population
table(POPINFO$Population)

##  Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
population=POPINFO$Population


## PLOT the data
tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(population)[match(pca$sample.id, sample.id)],
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

jpeg('PCA.jpg')
plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), xlab="eigenvector 2", ylab="eigenvector 1", main="PCA using HAPMAP PHASE3 DATA")
legend("bottomright", legend=levels(tab$pop), pch="o", col=1:(nlevels(tab$pop)))
dev.off()


