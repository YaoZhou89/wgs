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
    public GenerateScripts(){
        
    }
    public GenerateScripts(String inFile,String outFile,String model,String suffix){
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
            this.getMapping(inFile,outFile,suffix);
        }
        if(model.equals("model5")){
            this.getSh(inFile);
        }
        if(model.equals("model6")){
            System.out.println("modle6");
            this.getVar(inFile,outFile,suffix);
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
        if(model.equals("model9s")){
            this.getBamIndexSmall(inFile, outFile);
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
        if(model.equals("model13")){
            this.check(inFile);
        }
        if(model.equals("model14")){
            this.getblastp(inFile,outFile);
        }
        if(model.equals("model15")){
            this.bamIndex(inFile, outFile);
        }
        if(model.equals("model16")){
            getPara(inFile,suffix,outFile);
        }
        if(model.equals("model17")){
            checker(inFile,outFile);
        }
        if(model.equals("model18")){
            checker2file(inFile,outFile);
        }
        if(model.equals("model19")){
            CLUMPPFile(inFile,suffix,outFile);
        }
    }
    
    public GenerateScripts(String inFile, String outFile, String chrs){
        getByChr(inFile,outFile,chrs);
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
    public void check(String inFile){
        File test = new File(inFile);
        File[] t = IOUtils.listRecursiveFiles(test);
        File[] baiFile = IOUtils.listFilesEndsWith(t, ".bai");
        File[] bamFile = IOUtils.listFilesEndsWith(t, ".bam");
        Set bai = new HashSet();
        int i = 0;
        for(File f : baiFile){
            String a = f.toString();
            String b = a.replace(".bai","");
            bai.add(b);
//            if(i == 0){
//                System.out.println(b);
//                i++;
//            }
        }
        for(File f : bamFile){
//            if(i==1){
//                System.out.println(f.toString());
//                i++;
//            }
            if(bai.add(f.toString())){
                System.out.println(f.toString());
            }
        }
    }
    public void getMapping(String inFile,String outFile,String ref){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        header.append("module load bwa\n");
        File dir = new File(outFile+"/Scripts_ref");
        if(!dir.exists()) dir.mkdir();
        File Bamdir = new File(outFile+"/ref_bam");
        if(!Bamdir.exists()) Bamdir.mkdir();
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".fq.gz");
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
                bws.write("yhrun -n 1 -c 24 speedseq align -R\t"+"\"@RG\\tID:"+name+"\\tSM:"+name+"\\tPL:ILLUMINA\\tPI:330\""+
                        "\t-o\t"+ Bamdir.toString()+"/"+name  +"\t-t 24 "+ ref +"\t"
                        +fq[0].getAbsolutePath()+"\t"+fq[1].getAbsolutePath());
                bws.flush();
                bws.close();
            }
            
        } catch (IOException ex) {
                ex.printStackTrace();
        }
    }
    public void getVar(String inFile,String outFile,String suffix){
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
        command.append("-o ");
        command.append(outFile+ " ");
        command.append("-t 24 ");
        command.append(suffix);
        command.append(" ");
        
        System.out.println(command.toString());
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
        header.append("module load samtools\n");
        File dir = new File(outFile+"/bam_index_sh");
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
                if(entry.length()/1024 < 100) continue;
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
    public void getBamIndexSmall(String inFile, String outFile){
        StringBuilder header = new StringBuilder();
        header.append("#!/bin/bash\n");
        header.append("module load samtools\n");
        File dir = new File(outFile+"/bam_index_small_sh");
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
                if(entry.length()/1024/1024 >= 100) continue;
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
        BufferedWriter bws = IOUtils.getTextWriter(dir.toString()+"/"+name+".script");
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
                bws.write("yhbatch -N 1 -p paratera \t"+entry);
                bws.newLine();
                bws.write("sleep 2s\n");
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
    public void getByChr(String inFile, String outFile,String chrs){
        File test = new File(inFile);
        File dir1 = new File(outFile+"/splitSh");
        if(!dir1.exists()) dir1.mkdir();
        File dir2 = new File(outFile+"/splitBam");
        if(!dir2.exists()) dir2.mkdir();
       
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".bam");
        for (File f : subFs){
            if(f.length()/1024/1024 < 100) continue;
            String[] pa = f.getAbsolutePath().split("/");
            String name = pa[pa.length -1].replace(".bam", "");
//            System.out.println(name);
            BufferedWriter bw = IOUtils.getTextWriter(outFile+"/splitSh/" + name + ".bamSplitChr.sh");
            getBody(bw,chrs,f.getAbsolutePath(),outFile);
        }
    }
    private void getBody(BufferedWriter bw,String chrs,String name,String outFile){
        try {
            bw.write("#!/bin/bash\n");
            bw.write("module load samtools\n");
            String[] chr = chrs.split(" ");
            String[] pa = name.split("/");
            String names = pa[pa.length -1].replace(".bam", "");            
            for (int i = 0; i < chr.length; i++){
                bw.write("yhrun -n 1 -c 1 samtools view -b " + name +" "+ chr[i] + " > " + outFile +"/splitBam/"+ names + "." + chr[i] +".bam &\n");
            }
            bw.write("wait\n");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getSamtoolsCalling(String inFile,String outFile, String suffix,String chrs){
        try {
            StringBuilder header = new StringBuilder();
            header.append("#!/bin/bash\n");
            header.append("module load samtools\n");
            header.append("module load bcftools\n");
           
           
            String[] chr = chrs.split(" ");
            BufferedWriter bw = IOUtils.getTextWriter(outFile+"/vcf.sh");
            bw.write(header.toString());
            bw.newLine();
            for (int i = 0; i < chr.length; i++){
                File test = new File(inFile);
                File[] t = IOUtils.listRecursiveFiles(test);
                File[] tt = IOUtils.listFilesEndsWith(t, ".bam");
                File[] subF = IOUtils.listFilesContains(tt, chr[i]);
                StringBuilder command = new StringBuilder();
                command.append("yhrun -n 1 -c 1 samtools mpileup -uf ");
                command.append(suffix + " ");
                for (File f:subF){
                    command.append(" ");
                    command.append(f.getAbsoluteFile());
                }
                command.append(" | bcftools call -mv -Oz -o ");
                command.append(outFile);
                command.append("/");
                command.append(chr[i]+".vcf.gz &\n");
                bw.write(command.toString());
                
            }
            bw.write("wait\n");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getblastp(String inFile,String outFile){
        StringBuilder header = new StringBuilder();
        BufferedWriter bw = IOUtils.getTextWriter(outFile+"/"+1+".script");
        header.append("#!/bin/bash\n");
        File in = new File(inFile);
        File[] f = IOUtils.listRecursiveFiles(in);
        File[] fs = IOUtils.listFilesEndsWith(f,".fasta");
        int i = 0;
        StringBuilder command = new StringBuilder();
        int t = 1;
        try {
            for(File s :fs){
                i++;
                String p = s.getAbsolutePath().toString();
                String[] na = p.split("/");
                String name = na[na.length-1].split("\\.")[0];
                if(i % 23 == 0){
                    t++;
                    command.append("yhrun -n 1 -c 1 blastp -db /PARA/pp811/BIGDATA-2/data/wheat/orthologs/fasta/orthomcl -query ");
                    command.append(s.getAbsolutePath().toString());
                    command.append(" -seg yes -out ");
                    command.append(outFile+"/"+name+".blastout");
                    command.append(" -evalue 1e-5 -outfmt 7  -num_threads 1 ");
                    command.append("\n");
                    command.append("wait\n");
                    bw.write(header.toString());
                    bw.write(command.toString());
                    bw.flush();
                    bw.close();
                    command = new StringBuilder();
                    bw = IOUtils.getTextWriter(outFile+"/"+t+".script");
                }else{
                    command.append("yhrun -n 1 -c 1 blastp -db /PARA/pp811/BIGDATA-2/data/wheat/orthologs/fasta/orthomcl -query ");
                    command.append(s.getAbsolutePath().toString());
                    command.append(" -seg yes -out ");
                    command.append(outFile+"/"+name+".blastout");
                    command.append(" -evalue 1e-5 -outfmt 7  -num_threads 1 &");
                    command.append("\n");
                    command.append("sleep 5s\n");
                }
            }
            bw.write(header.toString());
            bw.write(command.toString());
            bw.flush();
            bw.close();  
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }   
    }
    public void bamIndex(String inFile,String outFile){
        File f = new File(inFile);
        File[] fs = IOUtils.listRecursiveFiles(f);
        File[] bfs = IOUtils.listFilesEndsWith(fs, ".bam");
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        int j = 0;
        try {
            for( File subFile : bfs){
                j++;
                if(j % 155 == 0){
                    bw.write("samtools index " + subFile.getAbsolutePath().toString());
                }else{
                    bw.write("samtools index " + subFile.getAbsolutePath().toString() + " &");
                }
                
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
                Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
            }
    }
    public void getPara(String inFile,String inFile2,String outFile){
        try {
            File f1 = new File(inFile);
            File f2 = new File(inFile2);
            File[] ff1 = IOUtils.listRecursiveFiles(f1);
            File[] ff2 = IOUtils.listRecursiveFiles(f2);
            File[] sf1 = IOUtils.listFilesEndsWith(ff1, ".bam");
            File[] sf2 = IOUtils.listFilesEndsWith(ff2, ".bam");
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            bw.write("Taxa\tReference\tBamPath\n");
            for(File subFile : sf1){
                String[] path = subFile.getAbsolutePath().toString().split("/");
                String name = path[path.length-1].split("\\.")[0];
                bw.write(name+"\t");
                if(name.startsWith("A0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/a_iwgscV1.fa");
                    bw.write("\t");
                }else if(name.startsWith("B0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/ab_iwgscV1.fa");
                    bw.write("\t");
                }else if(name.startsWith("B1")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/ab_iwgscV1.fa");
                    bw.write("\t");
                }
                else if(name.startsWith("D0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/d_iwgscV1.fa");
                    bw.write("\t");
                }else {
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/abd_iwgscV1.fa");
                    bw.write("\t");
                }
                bw.write(subFile.getAbsolutePath().toString());
                bw.newLine();
            }
            for(File subFile : sf2){
                if(subFile.toString().split(":").length>1) continue;
                String[] path = subFile.getAbsolutePath().toString().split("/");
                String name = path[path.length-1].split("\\.")[0];
                bw.write(name+"\t");
                if(name.startsWith("A0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/a_iwgscV1.fa");
                    bw.write("\t");
                }else if(name.startsWith("B0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/ab_iwgscV1.fa");
                    bw.write("\t");
                }else if(name.startsWith("B1")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/ab_iwgscV1.fa");
                    bw.write("\t");
                }else if(name.startsWith("D0")){
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/d_iwgscV1.fa");
                    bw.write("\t");
                }else{
                    bw.write("/data1/home/yaozhou/data/ref/wheat/genome/num/abd_iwgscV1.fa");
                    bw.write("\t");
                }
                bw.write(subFile.getAbsolutePath().toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void checker(String inFile, String outFile){
        File in = new File(inFile);
        Set <String> all = new HashSet();
        for ( int i = 1; i < 10 ; i++){
           all.add("A00"+i);
           all.add("B00"+i);
           all.add("TW00"+i);
           all.add("D00"+i);
        }
        for ( int i = 10; i < 31 ; i++){
            all.add("A0"+i);
            all.add("B0"+i);
            all.add("TW0"+i);
            all.add("D0"+i);
        }
        for ( int i = 31; i < 55 ; i++){
            all.add("A0"+i);
            all.add("B0"+i);
            all.add("TW0"+i);
        }
        for ( int i = 55; i < 92 ; i++){
            all.add("A0"+i);
            all.add("B0"+i);
        }
        for ( int i = 92; i < 99 ; i++){
            all.add("B0"+i);
        }
        for ( int i = 100; i < 126; i++){
            all.add("B"+i);
        }
        Set <String> fbam = new HashSet();
        Set <String> fbai = new HashSet();
        File[] subf = IOUtils.listRecursiveFiles(in);
        File[] bam = IOUtils.listFilesEndsWith(subf, ".bam");
        File[] bai = IOUtils.listFilesEndsWith(subf, "bai");
        for (File ff : bam){
            String[] nameS = ff.getAbsoluteFile().toString().split("/");
            String name = nameS[nameS.length-1].split("\\.")[0];
            fbam.add(name);
        }
        for (File ff : bai){
            String[] nameS = ff.getAbsoluteFile().toString().split("/");
            String name = nameS[nameS.length-1].split("\\.")[0];
            fbai.add(name);
        }
//        for(String a : fbai){
//            System.out.println(a);
//        }
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        try{
            for(String a : all){
//                System.out.println(a);
                if(fbam.add(a)){
                    bw.write("bam is missing: "+ a);
                    bw.newLine();
                }
                if(fbai.add(a)){
                    bw.write("bai is missing: " + a);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            
        }
    }
    
    public void checker2file(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inFile);
            BufferedReader br1 = IOUtils.getTextGzipReader(outFile);
            String[] te = null;
            List< String[]> geno = new ArrayList();
            String temp = "";
            int s1 = 0, s2 = 0;
            String[] taxa = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#C")) taxa = temp.split("\t");
                if(temp.startsWith("#")) continue;
                s1++;
                te = temp.split("\t");
                String[] ge = new String[te.length-9];
                for(int i = 9; i < te.length ; i++){
                    ge[i-9] = te[i];
                }
                geno.add(ge);
            }
            System.out.println("SNPs in file1 is: " + s1);
            int diff = 0;
            while((temp = br1.readLine())!=null){
                if(temp.startsWith("#")) continue;
               
                te = temp.split("\t");
                String[] ge = new String[te.length-9];
                for(int i = 9; i < te.length ; i++){
                    if(geno.get(s2)[i-9].startsWith("0/0")){
                        if(te[i].startsWith("1")) {
                            diff++;
                            System.out.println(taxa[i]);
                            System.out.println(geno.get(s2)[i-9]);
                            System.out.println(te[1]+"\t"+te[i]);
                            break;
                        }
                    }else if(geno.get(s2)[i-9].startsWith("1/1")){
                        if(te[i].startsWith("0/0")) diff++;
                    }
                }
                 s2++;
            }
            System.out.println("SNPs in file2 is: " + s2);
            System.out.println("Different sites number are : " + diff);
            
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void CLUMPPFile(String inFile, String suffix,String outFile){
        try {
            BufferedReader br1 = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            while ((temp = br1.readLine())!=null){
                if(temp.startsWith("K")){
                    bw.write("K "+suffix);
                    bw.newLine();
                }else{
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GenerateScripts.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
