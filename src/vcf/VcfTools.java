/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
//import java.lang.Object;
//import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;


/**
 *
 * @author yaozhou
 */
public class VcfTools {
    int dim = 500;
    public VcfTools(String inFile,String outFile){
        this.getABD(inFile,outFile);
    }
    public VcfTools(String inFile,Integer ExSize){
        this.splitByChr(inFile,ExSize);
    }
    public VcfTools(String inFile,String outFile,String model,int dim){
        this.dim = dim;
        if(model.equals("depth")){
            this.getDepthAll(inFile,outFile);
        }
        if(model.equals("het")){
            this.getHet(inFile,outFile);
        }
        if(model.equals("maf")){
            this.getMAF(inFile);
        }
        if(model.equals("MQ")){
            this.getMQ(inFile);
        }
        if(model.equals("GQ")){
            this.getGQ(inFile);
        }
        if(model.equals("eachDepth")){
            this.getDepthAll(inFile,outFile);
        }
        if(model.equals("vcfToStructure")){
            this.getStructure(inFile);
        }
        if(model.equals("vcfToXPCLR")){
            this.getXPCLR(inFile);
        }
        if(model.equals("DPRankSum")){
            this.filterByDPRankSum(inFile,outFile);
        }
        
    }
    public VcfTools(String inFile, String outFile, String suffix ){
        this.mergeVCF(inFile,outFile,suffix);
    }
    public VcfTools(String anchor, String inFile, String outFile, int minComp, double maxIBDDist, 
        int windowSize, int numThreads, int bestContrasts){
        ibdfilter.PfilterBasedOnIBD(anchor, inFile, outFile, minComp, maxIBDDist, windowSize, numThreads, bestContrasts);
    }
    public VcfTools(String inFile, String outFile,int windowSize, double threshold){   
        LDFilter.calLD(inFile,outFile,windowSize,threshold);
    }
     
   
    // type : filtered
    // parameters: MQ, FS, MQRankSum, and ReadPosRankSum
    public VcfTools(String inFile,String outFile,String MQ,String FS,String MQRankSum,String ReadPosRankSum,String BSQRankSum){
        this.getFilterd(inFile,outFile,MQ,FS,MQRankSum, ReadPosRankSum,BSQRankSum);
    }
 
    public VcfTools(String inFile,String outFile,String model,int size,int SNPnum,int header){
        this.getSub(inFile,outFile,size,SNPnum,header);
    }
    public VcfTools(String inFile,String outFile, double a, double b, double sd,
        int mindepth,int maxdepth,int dim){
        this.dim = dim;
        this.getDepthFilterd(inFile,outFile,a,b,sd,mindepth,maxdepth);
    }
    
    public void mergeVCF(String inFile, String outFile, String suffix){
        BufferedReader br;
        BufferedWriter bw;
        File test = new File (inFile);
        File[] fs = IOUtils.listRecursiveFiles(test);
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
        File[] subFs = IOUtils.listFilesEndsWith(fs, suffix);
        try{
            if(outFile.endsWith(".gz"))  bw = IOUtils.getTextGzipWriter(outFile);
            else bw = IOUtils.getTextWriter(outFile);
            boolean head = true;
            for (File f : subFs){
                br = IOUtils.getTextReader(f.toString());
                String temp = "";
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("#")){
                        if(head) {
                            bw.write(temp);
                            bw.newLine();
                        }
                    }else{
                       head = false;
                       bw.write(temp); 
                       bw.newLine();
                       bw.flush();
                    }
                }
                br.close();
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
    private void filterByDPRankSum(String inFile, String outFile){
        
        try {
            BufferedReader br ;
            if(inFile.endsWith(".gz")) br = IOUtils.getTextGzipReader(inFile);
            else br = IOUtils.getTextReader(inFile);
            String temp ;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private void getXPCLR(String inFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String outFile = inFile.replace(".vcf", ".geno");
            String map = inFile.replace(".vcf", ".snp");
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = null;
            while((temp = br.readLine())!=null){
                if(!temp.startsWith("#")){
                    String[] tmp = temp.split("\t");
                }
            }
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    private void getStructure(String inFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String outFile = inFile.replace("_in", "_out");
            String temp= null;
            StringBuilder A1= new StringBuilder();
            StringBuilder A2 = new StringBuilder();
            int line = 0;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while ((temp = br.readLine())!=null){
                
              if(line<2){
                  bw.write(temp);
                  bw.newLine();
              }else{
                A1= new StringBuilder();
                A2 = new StringBuilder();
                A1 = getOdd(temp,0);
                A2 = getOdd(temp,1);
//                System.out.println(A1.toString().length());
                bw.write(A1.toString());
                bw.newLine();
                bw.write(A2.toString());
                bw.newLine(); 
              }
              line++;
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private StringBuilder getOdd(String temp,int dir){
        // dir = 0 for odd and 1 for even
        String[] tmp = temp.split(" ");
        StringBuilder A = new StringBuilder(tmp[0]);
        int l = tmp.length;
//        if (l > 10000) l = 10000;
        for(int i = dir+2; i < l; i+=2 ){
            A.append("\t");
            A.append(tmp[i]);
        }
        return A;
    }
    private void getMQ(String inFile){
        try {
            System.out.println("Calculating MQ...");
            BufferedReader vcf;
            String outFile;
            if(inFile.endsWith("gz")){
                vcf = IOUtils.getTextGzipReader(inFile);
                outFile = inFile.replace(".vcf.gz","mq");
            }
            else {
                vcf = IOUtils.getTextReader(inFile);
                outFile = inFile.replace(".vcf",".mq");
            }
            String temp = null;
            String mq;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                  mq = temp.split("MQ=")[1].split(";")[0];
                  bw.write(mq);
                  bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private void getFilterd(String inFile,String outFile,String MQ, String FS, String ReadPosRankSum
    , String MQRankSum , String BSQRankSum){
        try {
            System.out.println("Filtering by FS, MQ, MQRankSum,ReadPosRankSum...");
            BufferedReader vcf;
            if(inFile.endsWith("gz")){
                vcf = IOUtils.getTextGzipReader(inFile);
//                outFile = inFile.replace(".vcf.gz",".filtered.vcf");
            }
            else {
                vcf = IOUtils.getTextReader(inFile);
//                outFile = inFile.replace(".vcf",".filtered.vcf");
            }
            String mq = "70", fs = "60", mqranksum = "-8", readposranksum = "-8",
                    bqrs = "0";
            String temp = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    try{
                      mq = temp.split("MQ=")[1].split(";")[0];
                    }catch (Exception ex){
                    };
                    if(Double.parseDouble(mq) < Double.parseDouble(MQ)) continue;
                    try{
                        fs = temp.split("FS=")[1].split(";")[0];}
                    catch (Exception ex){
                    };
                    if(Double.parseDouble(fs) > Double.parseDouble(FS)) continue;
                    try{
                        readposranksum = temp.split("ReadPosRankSum=")[1].split(";")[0];}
                    catch (Exception ex){
                    };
                    if(Double.parseDouble(readposranksum) < Double.parseDouble(ReadPosRankSum)) continue;
                    try{
                        mqranksum = temp.split("MQRankSum=")[1].split(";")[0];}
                    catch (Exception ex){
                    };
                    
                    if(Double.parseDouble(mqranksum) < Double.parseDouble(MQRankSum)) continue;
                    try{
                        bqrs = temp.split("BaseQRankSum=")[1].split(";")[0];}
                    catch (Exception ex){
                    };
                    if(Double.parseDouble(bqrs) < Double.parseDouble(BSQRankSum)) continue;
                    bw.write(temp);
                    bw.newLine();
                }else{
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private void getGQ(String inFile){
        try {
            System.out.println("Calculating GQ...");
            BufferedReader vcf;
            String outFile;
            if(inFile.endsWith("gz")){
                vcf = IOUtils.getTextGzipReader(inFile);
                outFile = inFile.replace(".vcf.gz","GQ");
            }
            else {
                vcf = IOUtils.getTextReader(inFile);
                outFile = inFile.replace(".vcf",".GQ");
            }
            String temp = null;
           
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            Integer indNum = 1;
            String[] tmp;
            String[] m=null;
            while((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    StringBuilder mq= new StringBuilder();
                    tmp = temp.split("\t");
                    indNum = tmp.length  ;
                    for(int i = 0;i<tmp.length-1;i++){
                        if(!tmp[i+9].startsWith(".")){
                             mq.append(m[3]+"\t");
                        }else{
                            m = tmp[i+9].split(":");
                            mq.append("0\t");
                        }
                    }
                    if(!tmp[tmp.length-1].startsWith(".")){
                        m = tmp[tmp.length-1].split(":");
                        mq.append(m[3]);
                    }else {
                        mq.append("0");
                    }
                    bw.write(mq.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
//    private void getDepth(String inFile,String outFile,int dim){
//        try {
//            System.out.println("Now analyzing SNP depth...");
//            BufferedReader vcf;
//            BufferedWriter VcfDepth;
//            BufferedWriter rsdFile;
//            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
//            else  vcf = IOUtils.getTextReader(inFile);
//            VcfDepth = IOUtils.getTextWriter(outFile+".all");
//            rsdFile = IOUtils.getTextWriter(outFile+".rsd");
//            String temp = null;
//            String[] tem = null;
////            Set depth = new HashSet(); // store the depth catalog
//            Integer d = 0;
//            int depthArray[][] = new int[101][300];// store the depth for each catalog
//            int snp = 0;
//            int sampleNum = 0;
//            int depthAll[] = new int[dim+1];
//            while ((temp = vcf.readLine())!=null){
//                if(!temp.startsWith("#")){
//                    snp++;
//                    if(snp%1000000 == 0) System.out.println("Anlyzing " + snp +"...");
//                    tem = temp.split("\t");
//                    sampleNum = tem.length - 9;
////                    int dp =0;
//                    double rsd = 0;
//                    if(tem[4].length()==1){
////                        for (int i = 9 ; i < tem.length; i++){
////                            String re = tem[i];
////                            if(re.split(":").length < 3){
//////                                System.out.println(temp);
////                                continue;
////                            }
////                            d = Integer.parseInt(re.split(":")[2]);
////                            if (d > 100) d = 100;
////                            depthArray[d][i-9]++;
////                        } 
//                        double avg = 0;
//                        double dep[] = new double[tem.length-9];
//                        for (int i = 9; i < tem.length; i++){
//                            String re = tem[i];
//                            if(re.split(":").length<3) continue;
//                            d = Integer.parseInt(re.split(":")[2]);
//                            dep[i-9] = d;
//                        } 
////                        dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
//                        avg = dp/sampleNum;
//                        double var=0;
//                        for (int i = 0; i<dep.length; i++){
//                            var+=(dep[i]-avg)*(dep[i]-avg);
//                        }
////                        rsd = Math.sqrt(var)/avg;
//                        rsd = Math.sqrt(var);
//                        if (dp > dim) dp = dim;
//                        depthAll[dp]++;
//                    }
//                    rsdFile.write(Integer.toString(dp)+"\t" + Double.toString(rsd));
//                    rsdFile.newLine();
//                }
//            }
//            rsdFile.flush();
//            rsdFile.close();
////            for (int i = 0; i < 101;i++){
////                VcfDepth.write(Integer.toString(i));
////                for (int j = 0; j < sampleNum;j++){
////                    VcfDepth.write("\t" + depthArray[i][j]);
////                    VcfDepth.flush();
////                }
////                VcfDepth.newLine();
////            }
////            VcfDepth.close();
//            for (int i = 0; i < dim+1; i++){
//                VcfDepth.write(Integer.toString(i));
//                VcfDepth.write("\t" + depthAll[i]);
//                VcfDepth.newLine();
//            }
//            VcfDepth.flush();
//            VcfDepth.close();
//        } catch (IOException ex) {
//           ex.printStackTrace();
//        }
//    }
    private void getABD(String inFile, String outFile) {
        BufferedReader vcf;
        BufferedWriter nvcf,bed,stat,Avcf;
        if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
        else  vcf = IOUtils.getTextReader(inFile);
        nvcf = IOUtils.getTextWriter(outFile + ".all.vcf");
        Avcf = IOUtils.getTextGzipWriter(outFile + ".A.vcf.gz");
        bed = IOUtils.getTextWriter(outFile + ".bed");
        stat = IOUtils.getTextWriter(outFile + ".stat");
        String temp = null;
        String[] temps = null;
        int a = 0, b = 0, d = 0;
        String chr = "0";
        boolean wrt = false;
        try {
            while((temp = vcf.readLine())!=null){
                if(temp.startsWith("#")){
                    nvcf.write(temp + "\n");
                    Avcf.write(temp+"\n");
                }else{
                    temps = temp.split("\t");
                    if(temps[0].substring(1,2).equals("A")){
//                        a++;
                        nvcf.write(temp + '\n');
                        Avcf.write(temp+"\n");
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }else if(temps[0].substring(1,2).equals("B")){
//                        b++;
                        nvcf.write(temp + '\n');
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }else if(temps[0].substring(1,2).equals("D")){
//                        d++;
                        nvcf.write(temp + '\n');
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }
                    if(wrt){
                        if(chr.equals("0")){
                            chr = temps[0];
                            a++;
                        }else{
                            if(chr.equals(temps[0])){
                                a++;
                            }else{
                                stat.write(chr+"\t"+a+"\n");
                                a = 1;
                                chr = temps[0];
                            }
                        }
                    }
                }
            }
            stat.write(chr+"\t"+a);
            stat.flush();
            stat.close();
            bed.flush();
            bed.close();
            nvcf.flush();
            nvcf.close();
            Avcf.flush();
            Avcf.close();
        } catch (IOException ex) {
            System.out.println("Reading vcf file failed!");
        }
    }
    private void splitByChr(String inFile, Integer ExSize) {
        BufferedReader vcf;
        BufferedWriter nvcf,bed,stat,Avcf;
        String outFile;
        String end = inFile.split("/")[inFile.split("/").length-1];
        String Path = inFile.replace(end,"");
        
        File theDir = new File(Path+"/ByChr/");
        if(inFile.endsWith("gz")) {
            vcf = IOUtils.getTextGzipReader(inFile);
            end =  end.replace(".vcf.gz","");
            outFile = Path +end;
        }else {
            vcf = IOUtils.getTextReader(inFile);
            end =  end.replace(".vcf","");
            outFile = Path + end;
        }
        
        if(!theDir.exists()) theDir.mkdir();
        stat = IOUtils.getTextWriter(theDir.toString() +"/"+ end+ ".stat");
        nvcf = null;
        String temp = null;
        String[] temps = null;
        int SNPnum = 0;
        String chr="chr";
        StringBuilder annotation = new StringBuilder();
        boolean wrt = false;
        int size = 0;
        
        try {
            while((temp = vcf.readLine())!=null){
                if(temp.startsWith("#")){
                    annotation.append(temp);
                    annotation.append("\n");
                }else{
                    temps = temp.split("\t");
                    if(temps[0].equals(chr)){
                        SNPnum++;
                        if(SNPnum%ExSize == 0) {
                            size++;
                            nvcf.flush();
                            nvcf.close();
                            
//                            stat.flush();
//                            stat.close();
                            System.out.println("Chrmosome: " + chr + ","+ size +"*"+ExSize);
                            nvcf = IOUtils.getTextGzipWriter(theDir.toString()+"/"+end+".chr"+chr+".temp"+size+".vcf.gz");
                            nvcf.write(annotation.toString());
                        }
                        nvcf.write(temp);
                        nvcf.newLine();
                    }else{
                        if(SNPnum > 0){
                            stat.write(chr+"\t"+SNPnum+"\n");
//                            stat.newLine();
                        }
                        SNPnum = 1;
                        size = 1;
                        chr = temps[0];
                        System.out.println("Chrmosome: " + chr + ","+ size +"*"+ExSize);
                        nvcf = IOUtils.getTextGzipWriter(theDir.toString()+"/"+end+".chr"+chr+".temp"+size+".vcf.gz");
                        nvcf.write(annotation.toString());
                        nvcf.write(temp);
                        nvcf.newLine();
                    }
                }
            }
            stat.write(chr+"\t"+SNPnum);
            stat.flush();
            stat.close();
            nvcf.flush();
            nvcf.close();
        } catch (IOException ex) {
            System.out.println(size);
            ex.printStackTrace();
        }
    }
    private void getDepthFilterd(String inFile,String outFile, double a, double b, double sd,
            int mindepth,int maxdepth) {
        try {
            System.out.println("Now Filtering SNP based on depth...");
            BufferedReader vcf;
            BufferedWriter VcfDepth;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            outFile = inFile.replace(".vcf",".depth.vcf");
            VcfDepth = IOUtils.getTextWriter(outFile);
            String temp = null;
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            Integer d = 0;
            int snp = 0;
            int sampleNum = 0;
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snp++;
                    if(snp%1000000 == 0) System.out.println("Anlyzing " + snp +"...");
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    int dp =0;
                    double rsd = 0;
                    int meanCUT = (int) ((sd-a)/b);
                    if(tem[4].length()==1){
//                        dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
                        double dep[] = new double[sampleNum];
                        for (int i = 9; i < tem.length; i++){
                            String re = tem[i];
                            if(re.split(":").length<3) continue;
                            if(!re.split(":")[0].startsWith(".")){
                                d = Integer.parseInt(re.split(":")[2]);
                            }else{
                                d = 0;
                            }
                            dep[i-9] = d;
                        }
                        for(int i = 0; i< dep.length;i++){
                            dp+=dep[i];
                        }
                        if(dp < mindepth ) continue;
                        if (dp > maxdepth) continue;
//                        StandardDeviation rrsd = new StandardDeviation();
                        rsd = getSD(dep);
//                        double rsd1 = getSD(dep);
                        if(rsd > sd) continue;
                        if(dp > meanCUT){
                            VcfDepth.write(temp);
                            VcfDepth.newLine();
                            VcfDepth.flush();
                        }else {
                            double cut = a+b*dp;
                            if (rsd < cut){
                                VcfDepth.write(temp);
                                VcfDepth.newLine();
                                VcfDepth.flush();
                            }
                        }
                    }
                }else{
                    VcfDepth.write(temp);
                    VcfDepth.newLine();
                    VcfDepth.flush();
                }
            }
            VcfDepth.flush();
            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    /*
    * remove SNPs by IBS
    *
    */
    private void removeByIDS(String inFile,String outFile){
        System.out.println("Now removing by IBS....");
        BufferedReader vcf;
        if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
        else  vcf = IOUtils.getTextReader(inFile);
        BufferedWriter ibs = IOUtils.getTextWriter(outFile);
        String temp = null;
        StringBuilder snp = new StringBuilder();
    }
    private double getSD(double dep[]){
        double sd=0;
        double sum = 0;
        for(int i = 0; i< dep.length;i++){
            sum+=dep[i];
        }
        double mean = sum/dep.length;
        sum =0;
        for(int i = 0; i< dep.length;i++){
           sum += ((dep[i]-mean)*(dep[i]-mean));
        }
        sd = Math.sqrt(sum/(dep.length-1));
        return sd;
    }
    private void getHet(String inFile,String outFile){
        try {
            System.out.println("Now analyzing the heterozyous....");
            BufferedReader vcf;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            String hetFile = inFile+".het";
            if(!outFile.endsWith(".vcf")) outFile = outFile+".het.vcf";
            BufferedWriter vcfhet = IOUtils.getTextWriter(outFile);
            BufferedWriter het = IOUtils.getTextWriter(hetFile);
            BufferedWriter het1 = IOUtils.getTextWriter(inFile+".perSNP.het");
//            BufferedWriter het2 = IOUtils.getTextWriter(inFile + ".IndHet");
            String temp = null;
            String tem[] = null;
            int a0 = 0, a1 = 0, a2 = 0, an = 0;
            int num = 0;
            while ((temp=vcf.readLine())!=null){
                if(temp.startsWith("#")){
                    vcfhet.write(temp);
                    vcfhet.newLine();
                    vcfhet.flush();
                }else{
                    double hetrate = 0;
                    double hetNum = 0, TotalNum = 0;
                    
                    num++;
                    tem = temp.split("\t");
                    StringBuilder te = new StringBuilder();
//                    te = null;
                    te.append(tem[0]);
                    
                    for (int i = 1;i<9; i++){
                        te.append("\t");
                        te.append(tem[i]);
                    }
                    if(tem[4].length() >1) continue;
                    for (int i = 9; i < tem.length;i++){
                        if(tem[i].startsWith(".")){
                            an++;
                            tem[i] = "./.";
                        }else if(tem[i].startsWith("0/0")){
                            a0++;
                            TotalNum++;
                        }else if(tem[i].startsWith("1/1")){
                            a2++;
                            TotalNum++;
                        }else if(tem[i].startsWith("1/0")|(tem[i].startsWith("0/1"))){
                            a1++;
                            hetNum++;
                            TotalNum++;
                            tem[i] = "./.";
                        }else{
//                            tem[i] = "./.";
                            System.out.println(tem[i]);
                            System.out.println(temp);
                        }
                        te.append("\t");
                        te.append(tem[i]);
                    }
                    hetrate = hetNum/TotalNum;
                    het1.write(Double.toString(hetrate)+"\n");
//                    vcfhet.write(te.toString());
//                    vcfhet.newLine();
                }
            }
            het1.flush();
            het1.close();
            vcfhet.flush();
            vcfhet.close();
            het.write("Number of 00 is: "+ a0);
            het.newLine();
            het.write("Number of 10 is: "+ a1);
            het.newLine();
            het.write("Number of 11 is: "+ a2);
            het.newLine();
            het.write("Number of missing is: "+ an);
            het.newLine();
            het.write("Total SNPs is: " + num);
            het.newLine();
            het.flush();
            het.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private void getHetStat(String inFile,String outFile){
        try {
            System.out.println("Now analyzing the heterozyous....");
            BufferedReader vcf;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            String hetFile = inFile+".het";
//            if(!outFile.endsWith(".vcf")) outFile = outFile+".het.vcf";
//            BufferedWriter vcfhet = IOUtils.getTextWriter(outFile);
            BufferedWriter het = IOUtils.getTextWriter(hetFile);
            String temp = null;
            String tem[] = null;
            int a0 = 0, a1 = 0, a2 = 0, an = 0;
            int num = 0;
            while ((temp=vcf.readLine())!=null){
                if(temp.startsWith("#")){
//                    vcfhet.write(temp);
//                    vcfhet.newLine();
//                    vcfhet.flush();
                }else{
                    
                    num++;
                    tem = temp.split("\t");
                    StringBuilder te = new StringBuilder();
//                    te = null;
                    te.append(tem[0]);
                    
                    for (int i = 1;i<9; i++){
                        te.append("\t");
                        te.append(tem[i]);
                    }
                    if(tem[4].length() >1) continue;
                    for (int i = 9; i < tem.length;i++){
                        if(tem[i].startsWith("./.")|tem[i].startsWith(".")){
                            an++;
                            tem[i] = "./.";
                        }else if(tem[i].startsWith("0/0")){
                            a0++;
                        }else if(tem[i].startsWith("1/1")){
                            a2++;
                        }else if(tem[i].startsWith("1/0")|(tem[i].startsWith("0/1"))){
                            a1++;
                            tem[i] = "./.";
                        }else{
//                            tem[i] = "./.";
                            System.out.println(tem[i]);
                            System.out.println(temp);
                        }
                        te.append("\t");
                        te.append(tem[i]);
                    }
//                    vcfhet.write(te.toString());
//                    vcfhet.newLine();
                }
            }
//            vcfhet.flush();
//            vcfhet.close();
            het.write("Number of 00 is: "+ a0);
            het.newLine();
            het.write("Number of 10 is: "+ a1);
            het.newLine();
            het.write("Number of 11 is: "+ a2);
            het.newLine();
            het.write("Number of missing is: "+ an);
            het.newLine();
            het.write("Total SNPs is: " + num);
            het.newLine();
            het.flush();
            het.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private void getindHet(String inFile,String outFile){
        try {
            System.out.println("Now analyzing the heterozyous for each individual....");
            BufferedReader vcf;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            String hetFile = inFile+".IndHet";
//            if(!outFile.endsWith(".vcf")) outFile = outFile+".het.vcf";
//            BufferedWriter vcfhet = IOUtils.getTextWriter(outFile);
            BufferedWriter het = IOUtils.getTextWriter(hetFile);
            String temp = null;
            String tem[] = null;
            int num = 0;
            while ((temp=vcf.readLine())!=null){
                if(temp.startsWith("#")){
//                    vcfhet.write(temp);
//                    vcfhet.newLine();
//                    vcfhet.flush();
                }else{
                    num++;
                    tem = temp.split("\t");
                    StringBuilder te = new StringBuilder();
//                    te = null;
                    te.append(tem[0]);
                    
                    for (int i = 1;i<9; i++){
                        te.append("\t");
                        te.append(tem[i]);
                    }
                    if(tem[4].length() >1) continue;
                    for (int i = 9; i < tem.length;i++){
                        if(tem[i].startsWith("./.")|tem[i].startsWith(".")){
                            
                            tem[i] = "./.";
                        }else if(tem[i].startsWith("0/0")){
//                            a0[i-9]++;
                        }else if(tem[i].startsWith("1/1")){
//                            a2[i-9]++;
                        }else if(tem[i].startsWith("1/0")|(tem[i].startsWith("0/1"))){
//                            a1[i-9]++;
                            tem[i] = "./.";
                        }else{
//                            tem[i] = "./.";
                            System.out.println(tem[i]);
                            System.out.println(temp);
                        }
                        te.append("\t");
                        te.append(tem[i]);
                    }
//                    vcfhet.write(te.toString());
//                    vcfhet.newLine();
                }
            }
//            vcfhet.flush();
//            vcfhet.close();
//            het.write("Number of 00 is: "+ a0);
//            het.newLine();
//            het.write("Number of 10 is: "+ a1);
//            het.newLine();
//            het.write("Number of 11 is: "+ a2);
//            het.newLine();
//            het.write("Number of missing is: "+ an);
//            het.newLine();
//            het.write("Total SNPs is: " + num);
//            het.newLine();
//            het.flush();
//            het.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private void getMAF(String inFile){
        try {
            System.out.println("Now analyzing the minor allele frequencing....");
            BufferedReader vcf;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            String hetFile = inFile+".het";
            String outFile = inFile + ".maf";
            BufferedWriter maf = IOUtils.getTextWriter(outFile);
            String temp = null;
            String tem[] = null;
            
            while ((temp=vcf.readLine())!=null){
                int a0 = 0, a1 = 0, a2 = 0, an = 0;
                if(!temp.startsWith("#")){
                    tem = temp.split("\t");
                    for (int i = 9; i < tem.length;i++){
                        if(tem[i].startsWith("./.")){
                            an++;
                        }else if(tem[i].startsWith("0/0")){
                            a0++;
                        }else if(tem[i].startsWith("1/1")){
                            a2++;
                        }else if(tem[i].startsWith("1/0")|(tem[i].startsWith("0/1"))){
                            a1++;
//                            tem[i] = "./.";
                        }else{
                            System.out.println(tem[i]);
                        }
                    }
                    double m1 = (double)(a1+2*a2)/(2*(tem.length-9));
                    double m2 = 1 - m1;
                    if(m2 < m1) m1 = m2;
                    maf.write(Double.toString(m1));
                    maf.newLine();
                }
            }
            maf.flush();
            maf.close();
            
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    /* 
    * Objects: get the randomly sampled vcfs
    * Inputs: size, total number to be sampled
    * statFile: file containing the information of vcf;, if not provided;
    * Outputs: file with $size number vcfs, and evenly distributed in each chr.
    */
    private void getSub(String inFile,String outFile,int size,int SNPnum,int header){
        
        try {
            BufferedReader br;
            if(inFile.endsWith("gz"))  br = IOUtils.getTextGzipReader(inFile);
            else  br = IOUtils.getTextReader(inFile);
//            String outFile = inFile.replace(".vcf", ".sub.vcf");
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = null;
            String[] tem = null;
            int i = 0;
            Set pos = new HashSet();
            double rn = 0;
            Map po = new HashMap();
            for(int j = 0; j < size;j++){
                int a = (int)(1+Math.random()*size);
                pos.add(a);
                if(pos.size()==SNPnum) break;
                po.put(a,0);
            }
            StringBuilder headerS = new StringBuilder();
            boolean w = true;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#") ){
                    headerS.append(temp+"\n");
                }else{
                    if(w){
                        bw.write(headerS.toString());
                        w = false;
                    }
                    if (po.get(i)!=null){
                        bw.write(temp);
                        bw.newLine(); 
                    }
                }
                i++; 
                if(i%1000000==0){
                    System.out.println("Analyzing snps: "+i+" ....");
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    /*
    * Checking the four groups of A/AB
    * Methods: sample from each group, and then calculate the sd without missing
    *
    */
    private void getDepthAll(String inFile,String outFile){
        try {
            System.out.println("Now Testing SNP depth...");
            BufferedReader vcf;
            BufferedWriter VcfDepth;
            BufferedWriter rsdFile;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            VcfDepth = IOUtils.getTextWriter(outFile+".eachLocus");
            String temp = null;
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            Integer d = 0;
            int snp = 0;
            int sampleNum = 0;
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snp++;
                    if(snp%1000000 == 0) System.out.println("Anlyzing " + snp +"...");
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    int dp =0;
                    double rsd = 0;
                    StringBuilder depth = new StringBuilder();
                    String geno = null;
                    double[] deptheach = new double[sampleNum];
//                    double[] depthmean0 = new double[sampleNum];
                    if(tem[4].length()==1){
//                        depth.append(Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]));
                        for (int ch = 0 ;ch < sampleNum; ch ++){
                            if(tem[ch+9].startsWith(".")) {
                                geno = "0";
                            }else{
                                geno = tem[ch+9].split(":")[2];
                            }
                            depth.append(geno);
                            depth.append("\t");
                            deptheach[ch] = Integer.parseInt(geno);
//                            depthmean0[ch] = Double.parseDouble(geno);
                        }
                        
                        depth.append(getSD(deptheach));
                        VcfDepth.write(depth.toString());
                        VcfDepth.newLine();
                    }
                }
            }
            StringBuilder DM = new StringBuilder();
//            for (int i = 0; i< sampleNum -1;i++){
//                depthmean[i] = depthmean[i]/snp;
//                DM.append(depthmean[i]);
//                DM.append("\t");
//            };
//            depthmean[sampleNum-1] = depthmean[sampleNum-1]/snp;
//            DM.append(depthmean[sampleNum]);
            VcfDepth.write(DM.toString());
            VcfDepth.flush();
            VcfDepth.close();
//            for (int i = 0; i < 101;i++){
//                VcfDepth.write(Integer.toString(i));
//                for (int j = 0; j < sampleNum;j++){
//                    VcfDepth.write("\t" + depthArray[i][j]);
//                    VcfDepth.flush();
//                }
//                VcfDepth.newLine();
//            }
//            VcfDepth.close();
//            for (int i = 0; i < dim+1; i++){
//                VcfDepth.write(Integer.toString(i));
//                VcfDepth.write("\t" + depthAll[i]);
//                VcfDepth.newLine();
//            }
//            VcfDepth.flush();
//            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    
    
}
