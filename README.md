# wgs
This pipeline was designed for the wheat re-sequence projects, only java source code are posted here. You can use NetBeans to open this project.
## Usage:
You can only download dist/wgs.jar and run the command, you can also rebuild the jar file using NetBeans.
## 1. Pipeline for construction wheat variant map I (VMapI)
### Filtering SNPs
#### 1.1 Filtering by quality
java -jar fold-to-wgs/wgs.jar --model vcf --type quality --file fold-to-read/file.vcf --out fold-to-write/out.vcf [--MQ num] [--FS num] [--MQRankSum num] [--SOR num]
#### 1.2 Filtering by depth(min, max and SD, SD/depth)
java -jar 
#### 1.3 Filtering by Segeragation test

#### 1.4 Filtering by LD

#### 1.5 Filtering by IBD

#### 1.6 Filtering by maf (using vcftools)

### Spliting and merging vcf files
### Estimate pair-wise genetic distance from vcf file
##### firstly, estimate the heterozygousity using vcftools
##### then, estimate the pair-wise genetic divergency using wgs
    java -jar wgs.jar --model vcf --type GeneticDivergency --file fold-to-read/file.vcf --out fold-to-write/GeneticDivergency.txt
##### finaly, calculate pair-wise distance using R
  D = 1 - 0.25*(Pi1+Pi2)/Pi12
  
