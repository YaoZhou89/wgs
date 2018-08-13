# wgs
This pipeline was designed for the wheat re-sequence projects, only java source code are posted here. You can use NetBeans to open this project.
## Usage:
You can only download dist/wgs.jar and run the command, you can also rebuild the jar file using NetBeans.
### Filtering SNPs
#### Filtering by quality
java -jar fold-to-wgs/wgs.jar --model vcf --type quality --file fold-to-read/file.vcf --out fold-to-write/out.vcf [--MQ num] [--FS num] [--MQRankSum num] [--SOR num]
#### Filtering by depth(min, max and SD, SD/depth)
java -jar 
#### Filtering by Segeragation test

#### Filtering by LD

#### Filtering by IBD

### Spliting and merging vcf files
### Estimate pair-wise genetic distance from vcf file
##### firstly, estimate the heterozygousity using vcftools
##### then, estimate the pair-wise genetic divergency using wgs
    java -jar wgs.jar --model vcf --type GeneticDivergency --file fold-to-read/file.vcf --out fold-to-write/GeneticDivergency.txt
##### finaly, calculate pair-wise distance using R
  D = 1 - 0.25*(Pi1+Pi2)/Pi12
  
