/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package script;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class GenerateScripts {
    public GenerateScripts(String inFile,String outFile,String model){
        if(model.equals("model1")){
            this.getRScripts();
        }
        if(model.equals("model2")){
            runbgzip(inFile,outFile);
        }
        if(model.equals("model3")){
            this.getMacscript();
        }
        if(model.equals("model4")){
            this.getMapping(inFile,outFile);
        }
        if(model.equals("model5")){
            this.getSh(inFile);
        }
        if(model.equals("model6")){
            this.getVar(inFile,outFile);
        }
        if(model.equals("model7")){
            this.getBamHeader(inFile, outFile);
        }
        if(model.equals("model8")){
            this.getBam(inFile,outFile);
        }
        if(model.equals("model9")){
            this.getBamIndex(inFile, outFile);
        }
        if(model.equals("model10")){
            this.getS();
        }
        if(model.equals("model11")){
            this.makeDir(inFile);
        }
        if(model.equals("model12")){
            this.getConf(inFile,outFile);
        }
    }
    public void getRScripts(){
         BufferedWriter bw = IOUtils.getTextWriter("scripts.sh");
            try {
                for (int i = 10; i < 15;i++){
                   bw.write("nohup Rscript FullModle_AD.R "+ i+ " > ./log/AD_BLUPres_"+ i +".log &\n");
                }
                bw.flush();
                bw.close();
           
            } catch (IOException ex) {
                ex.printStackTrace();
            };
    }
   
    public void getMapping(String inFile,String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        File dir = new File(outFile+"/Scripts_ref");
        if(!dir.exists()) dir.mkdir();
        File Bamdir = new File(outFile+"/ref_bam");
        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, "q.gz");
        File[] reFs = getPath(subFs);   
        try {
            for(File entry : reFs){
                File[] fq = IOUtils.listRecursiveFiles(entry);
                if(fq.length!=2) {
                    System.out.println("more than two files!");
                }
                String[] pa = entry.toString().split("/");
                String name = pa[pa.length-1];
                BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+name+".sh");
                bws.write(header.toString());
                bws.write("yhrun -n 1 -c 24 speedseq align -R\t"+"\"@RG\\tID:"+name+"\\tSM:Seq01\\tPL:ILLUMINA\\tPI:330\""+
                        "\t-o\t"+ Bamdir.toString()+"/"+name  +"\t-t 24 /BIGDATA/pp811/data/tomato/ref/S_lycopersicum_chromosomes.3.00.fa\t"
                        +fq[0].getAbsolutePath()+"\t"+fq[1].getAbsolutePath());
                bws.flush();
                bws.close();
            }
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getVar(String inFile,String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        File dir = new File(outFile+"/Scripts_ref");
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".bam");
//        File[] reFs = getPath(subFs);   
        StringBuilder command =new StringBuilder("yhrun -n 1 -c 24 speedseq var ");
        try {
            for(File entry : subFs){
                if(entry.length()/1024/1024 < 100) continue;
//                String[] pa = entry.toString().split("/");
//                String name = pa[pa.length-1];
                command.append(entry.toString());
                command.append("\t");
            }
            BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+"var.sh");
            bws.write(header.toString());
            bws.write(command.toString());
            bws.newLine();
            bws.flush();
            bws.close();
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getBamHeader(String inFile, String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        File dir = new File(outFile+"/Scripts_lyc");
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".bam");
//        File[] reFs = getPath(subFs);  
        String outbam ="";
        StringBuilder command =new StringBuilder();
        int num = 0,f=0;
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+"change"+f+".sh");
        try {
            for(File entry : subFs){
                if(entry.length()/1024/1024 < 100) continue;
//                String[] pa = entry.toString().split("/");
//                String name = pa[pa.length-1];
                num++;
                
                if(num % 24 == 0){
                    f++;
                    bws.write(header.toString());
                    bws.write(command.toString());
                    bws.write("wait");
                    bws.newLine();
                    bws.flush();
                    bws.close();
                    bws = IOUtils.getTextWriter(dir.toString()+"/"+"change"+f+".sh");
                    command =new StringBuilder();
                }
                command.append("yhrun -n 1 -c 1 samtools view -L /BIGDATA/pp811/data/tomato/lycopersicoides/lyc.bed -o\t");
                outbam = entry.toString().replace(".bam", ".chr.bam");
                command.append(outbam);
                command.append("\t");
                command.append(entry.toString());
                command.append("\t&\n");
                
            }
            bws.write(header.toString());
            bws.write(command.toString());
            bws.write("wait");
            bws.newLine();
            bws.flush();
            bws.close();
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getBamIndex(String inFile, String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        File dir = new File(outFile+"/Scripts_lyc");
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".bam");
//        File[] reFs = getPath(subFs);  
        String outbam ="";
        StringBuilder command =new StringBuilder();
        int num = 0,f=0;
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+"index"+f+".sh");
        try {
            for(File entry : subFs){
                if(entry.length()/1024/1024 < 100) continue;
//                String[] pa = entry.toString().split("/");
//                String name = pa[pa.length-1];
                num++;
                
                if(num % 24 == 0){
                    f++;
                    bws.write(header.toString());
                    bws.write(command.toString());
                    bws.write("wait");
                    bws.newLine();
                    bws.flush();
                    bws.close();
                    bws = IOUtils.getTextWriter(dir.toString()+"/"+"index"+f+".sh");
                    command =new StringBuilder();
                }
                command.append("yhrun -n 1 -c 1 samtools index \t");
                outbam = entry.toString().replace(".bam", ".bam.bai");
                command.append(entry.toString());
                command.append("\t");
                command.append(outbam);
                command.append("\t&\n");
                
            }
            bws.write(header.toString());
            bws.write(command.toString());
            bws.write("wait");
            bws.newLine();
            bws.flush();
            bws.close();
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getBam(String inFile, String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        File dir = new File(outFile+"/Scripts");
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".chr.bam");
//        File[] reFs = getPath(subFs);  
        String outbam ="";
        StringBuilder command =new StringBuilder();
        int num = 0,f=0;
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+"cbam"+f+".sh");
        try {
            for(File entry : subFs){
                if(entry.length()/1024/1024 < 100) continue;
//                String[] pa = entry.toString().split("/");
//                String name = pa[pa.length-1];
                num++;
                
                if(num % 24 == 0){
                    f++;
                    bws.write(header.toString());
                    bws.write(command.toString());
                    bws.write("wait");
                    bws.newLine();
                    bws.flush();
                    bws.close();
                    bws = IOUtils.getTextWriter(dir.toString()+"/"+"cbam"+f+".sh");
                    command =new StringBuilder();
                }
                String[] pa = entry.toString().split("/");
                String name = pa[pa.length-1];
                name = name.replace(".chr.bam", "");
                outbam = entry.toString().replace(".chr.bam", ".header.bam");
                command.append("yhrun -n 1 -c 1 java -jar /BIGDATA/pp811/softwares/picard.jar AddOrReplaceReadGroups I=");
                command.append(entry.toString());
                command.append("\tO=");
                command.append(outbam);
                command.append("\tRGID=");
                command.append(name);
                command.append("\tRGSM=");
                command.append(name);
                command.append("\tRGLB=lib1 RGPL=illumina RGPU=u1");
                command.append("\t&\n");
                
            }
            bws.write(header.toString());
            bws.write(command.toString());
            bws.write("wait");
            bws.newLine();
            bws.flush();
            bws.close();
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void runbgzip(String inFile, String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        
        File dir = new File(outFile);
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".gz");
//        File[] reFs = getPath(subFs);  
        StringBuilder command =new StringBuilder();
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+"bgzip"+".sh");
        try {
            bws.write(header.toString());
            for(File entry : subFs){
                String[] names = entry.toString().split("/");
                String name = names[names.length-1].split("\\.")[0];
                for(int i = 1; i<43; i++){
                    bws.write("bgzip "+outFile+"/"+name+"/"+i+".vcf");
                    bws.newLine();
                    bws.write("tabix "+outFile+"/"+name+"/"+i+".vcf.gz");
                    bws.newLine();
                }
            }
            bws.flush();
            bws.close();
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public static void mergeLineage(String inFile, String outFile,String chr){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        header.append("gatk CombineGVCFs \\\n" +
"   -R /data1/home/yaozhou/data/ref/wheat/genome/num/ab_iwgscV1.fa \\\n");
        File dir = new File(outFile);
        if(!dir.exists()) dir.mkdir();
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".gz");
//        File[] reFs = getPath(subFs);  
        StringBuilder command =new StringBuilder();
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+chr+".sh");
        try {
            
            bws.write(header.toString());
            for(File entry : subFs){
                String[] names = entry.toString().split("/");
                String name = names[names.length-1].split("\\.")[0];
                if(!name.equals(chr)){
                    continue;
                }
                bws.write("-variant "+ entry.toString()+" \\");
                bws.newLine();
            }
            bws.write(" -O "+outFile+"/"+chr+".lineage.g.vcf.gz");
            bws.flush();
            bws.close();
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getMacscript(){
        StringBuilder header = new StringBuilder();
        File dir = new File("Scripts");
        if(!dir.exists()) dir.mkdir();
        String[] chr = {"A","AB","ABD","D"};
        String name = null;
        try {
            BufferedWriter bwr = IOUtils.getTextWriter(dir.toString()+"/wheat_step2.sh");
            for(int i = 1;i<5;i++){
                for(int j = 0;j<chr.length;j++){
                    name = "Step2_scripts_"+chr[j]+"_"+i+".script";
                    bwr.write("bsub < "+name +" &");
                    bwr.newLine();
                    bwr.newLine();
                    BufferedWriter bw = IOUtils.getTextWriter(dir.toString()+"/"+name);
                    bw.write(header.toString());
                    bw.write("vcftools --vcf ~/data/Evo/"+chr[j]+"/"+chr[j]+".all.depth.vcf "
                        + "--keep ~/data/Evo/"+chr[j]+"/group"+i+".txt"+"  --min-alleles 2 --max-alleles 2 "
                        +"--recode --recode-INFO-all --out ~/data/Evo/"+chr[j]+"/"+chr[j]+".group"+i);
                    bw.newLine();
                    bw.flush();
                    bw.close();   
                }
            }
            bwr.flush();
            bwr.close();
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public File[] getPath(File[] folder){
        TreeSet<File> fileTree = new TreeSet();
        TreeSet<File> subfile = new TreeSet();
        for (File entry : folder){
            if(entry.isHidden()) entry.delete();
            File path = new File(entry.getParent());
            fileTree.add(path);
        }
        return fileTree.toArray(new File[fileTree.size()]);
    }
    private void getS(){
        StringBuilder header = new StringBuilder("#BSUB -L /bin/bash\n" +
            "#BSUB -J barely_lastz\n" +
            "#BSUB -n 1\n" +
            "#BSUB -e %J.err\n" +
            "#BSUB -o %J.out\n" +
            "#BSUB -q vvl\n");
        String[] geno = {"A","B","D"};
        
        for(int i =1;i<8;i++){
            for(int j=0;j<3;j++){
                StringBuilder command = new StringBuilder();
                command.append(header);
                command.append("cd /public-supool/home/yaozhou/data/Genomes/barely\n");
                command.append("lastz /public-supool/home/yaozhou/data/ref/wheat/iwgsc_refseqv1.0_all_chromosomes/byChr/mask/chr"+
                        i+geno[j] +".fa "+ "*.fna --notransition --step=20 --nogapped --ambiguous=iupac " +
                        "--format=differences > barely_chr" + i+geno[j]+".fa.diff ");
                String outName = "barely_chr"+i+geno[j]+".sh";
                BufferedWriter sh = IOUtils.getTextWriter(outName);
                try {
                    sh.write(command.toString());
                    sh.newLine();
                    sh.flush();
                    sh.close();
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
                
            }
        }
    }
    private void getSh(String inFile){
        StringBuilder header = new StringBuilder();
        File dir = new File(inFile+"/Scripts");
        if(!dir.exists()) dir.mkdir();
        header.append("#!/bin/bash\n");
        if(!dir.exists()) dir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".sh");
        int i = 1;
        int j =0;
        String name ="s1";
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+name+".sh");
        try {
            bws.write(header.toString());
            for(File entry : subFs){
                j++;
                if(j%20==0){
                    bws.flush();
                    bws.close();
                    i++;
                    name = "s"+Integer.toString(i);
                    bws = IOUtils.getTextWriter(dir.toString()+"/"+name+".script");
                    bws.write(header.toString());
                }
                bws.write("sh\t"+entry);
                bws.newLine();
                
            }
            bws.flush();
            bws.close();
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void makeDir(String inFile){
        
//        File Bamdir = new File(outFile+"/ref_bam");
//        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".fa");
//        try {
            for(File entry : subFs){
                String name = entry.toString();
                String d = name.split(".fa")[0];
                File dir = new File(d);
                if(!dir.exists()) dir.mkdir();
            }
            
//        } 
    }
    private void getConf(String inFile,String outFile){
        StringBuilder header = new StringBuilder("GENETIC_CODE_TABLE=1\n" +
            "GENETIC_CODE_TABLENAME=Standard\n" +
            "MITO_GENETIC_CODE_TABLE=2\n" +
            "MITO_GENETIC_CODE_TABLENAME=Vertebrate Mitochondrial\n" +
            "\n" );
        String end ="# DBSNP_VCF_FILE=Homo_sapiens.vcf.gz\n" +
                    "#Running SIFT 4G\n" +
                    "SIFT4G_PATH=~/software/sift4g/bin/sift4g\n" +
                    "PROTEIN_DB=~/data/uniref90/uniref90.fasta\n" +
                    "COMPUTER=Lulab2\n" +
                    "\n" +
                    "# Sub-directories, don't need to change\n" +
                    "GENE_DOWNLOAD_DEST=gene-annotation-src\n" +
                    "CHR_DOWNLOAD_DEST=chr-src\n" +
                    "LOGFILE=Log.txt\n" +
                    "ZLOGFILE=Log2.txt\n" +
                    "FASTA_DIR=fasta\n" +
                    "SUBST_DIR=subst\n" +
                    "ALIGN_DIR=SIFT_alignments\n" +
                    "SIFT_SCORE_DIR=SIFT_predictions\n" +
                    "SINGLE_REC_BY_CHR_DIR=singleRecords\n" +
                    "SINGLE_REC_WITH_SIFTSCORE_DIR=singleRecords_with_scores\n" +
                    "DBSNP_DIR=dbSNP\n" +
                    "\n" +
                    "# Doesn't need to change\n" +
                    "FASTA_LOG=fasta.log\n" +
                    "INVALID_LOG=invalid.log\n" +
                    "PEPTIDE_LOG=peptide.log\n" +
                    "ENS_PATTERN=ENS\n" +
                    "SINGLE_RECORD_PATTERN=:change:_aa1valid_dbsnp.singleRecord\n";
        
        for(int i = 0; i < 43; i++){
            StringBuilder m = new StringBuilder( "PARENT_DIR=/data1/home/yaozhou/data/ref/wheat/SIFT/byChr/chr"
                    + String.format("%03d", i)+"\n"+
                "ORG=wheat\n" +
                "ORG_VERSION=iwgsc_v1.0_chr"+String.format("%03d", i));
            m.append("\n");
            String out = outFile+"/chr"+String.format("%03d", i)+"_config.txt";
            BufferedWriter bw = IOUtils.getTextWriter(out);
            try {
                bw.write(header.toString());
                bw.write(m.toString());
                bw.write(end);
                bw.flush();
                bw.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
            
        }
    }
}
