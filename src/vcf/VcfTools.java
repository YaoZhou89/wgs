/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import htsjdk.tribble.index.tabix.TabixIndex;
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
import math.SegeregationTest;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.variant.variantcontext.VariantContextUtils;
import java.text.DecimalFormat;
import java.util.Hashtable;

//import java.lang.Object;
//import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;


/**
 *
 * @author yaozhou
 */
public class VcfTools {
    int SNPsAll = 0;
    int SNPsFail = 0;
    public VcfTools(){
        
    }
    public VcfTools(String inFile,String outFile){
        this.getABD(inFile,outFile);
    }
    public VcfTools(String inFile,Integer ExSize){
        this.splitByChr(inFile,ExSize);
    }
    public VcfTools(String inFile,String outFile,String model,double se){
        if(model.equals("TotalDepthFilter")){
            filtTotalDepth(inFile, outFile,se);
        }else
        if(model.equals("het")){
            this.getHet(inFile,outFile);
        }else
        if(model.equals("maf")){
            this.getMAF(inFile);
        }else
        if(model.equals("MQ")){
            this.getMQ(inFile);
        }else
        if(model.equals("GQ")){
            this.getGQ(inFile);
        }else
        if(model.equals("checkDepth")){
            this.getDepthAll(inFile,outFile,se);
        }else
        if(model.equals("vcfToStructure")){
            this.getStructure(inFile);
        }else
        if(model.equals("vcfToXPCLR")){
            this.getXPCLR(inFile,outFile);
        }else if(model.equals("dpfilter")){
            this.dpFilter(inFile, outFile,se);
        }else if(model.equals("mac")){
            this.macFilter(inFile, outFile);
        }else if (model.equals("random")){
            getRandom(inFile,outFile,se);
        }else if (model.equals("phasing")){
            phaseVCF(inFile,outFile);
        }else if (model.equals("checkChr")){
            checkChr(inFile);
        }else if(model.equals("statGenotype")){
            genotypeStat(inFile,outFile);
        }else if(model.equals("phaseVCF")){
            phase(inFile,outFile);
        }else if (model.equals("vcfToEhhMap")){
            getEhhMap(inFile,outFile);
        }else if(model.equals("changeVCF")){
            changeVCF(inFile, outFile);
        }else if (model.equals("toFasta")){
            toFasta(inFile,outFile);
        }else if (model.equals("mergeABD")){
            mergeABD(inFile,outFile);
        }else if (model.equals("changeName")){
            changeName(inFile,outFile);
        }else if (model.equals("diagnostic")){
            diagSNP(inFile,outFile);
        }else if(model.equals("equalDepth")){
            equalDepth(inFile,outFile);
        
        }else if (model.equals("fastDiversity")){
            getDiversityFast(inFile,outFile);
        }else if (model.equals("getOneChr")){
            getOneChr(inFile,outFile);
        }
        
    }
    public void getOneChr(String inFile, String outFile){
        try {
            String temp = "";
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else if(temp.startsWith("Chr1\t")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    break;
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public VcfTools(String inFile, String bedFile,String outFile,Integer windowSize){
            getDiversity(inFile,bedFile,outFile,windowSize);
    }
    
    public void getDerivedSNP(String inFile,String inFile2, String outFile){
        try {
            TabixReader br1 = new TabixReader(inFile);
            TabixReader br2 = new TabixReader (inFile2);
            BufferedWriter bw1 = IOUtils.getTextWriter(outFile + ".s1.snp.txt");
            BufferedWriter bw2 = IOUtils.getTextWriter(outFile + ".s2.snp.txt");
            BufferedWriter bw12 = IOUtils.getTextWriter(outFile + ".s12.snp.txt");
            String temp = "";
            String[] te = null;
            Set<String> p1 = new HashSet();
            Set<String> pa = new HashSet();
            Set<String> p12 = new HashSet();
            System.out.println("Reading file 1...");
            while((temp = br1.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    String rs = "rs"+te[0]+"_"+te[1];
                    p1.add(rs);
                    pa.add(rs);
                }
            }
            System.out.println("Reading file 2...");
            while((temp = br2.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    String rs = "rs"+te[0]+"_"+te[1];
                    if(pa.add(rs)){
                        bw2.write(rs);
                        bw2.newLine();
                    }else{
                        p12.add(rs);
                        bw12.write(rs);
                        bw12.newLine();
                    }
                }
            }
            bw2.flush();
            bw2.close();
            bw12.flush();
            bw12.close();
            System.out.println("Writing file 1...");
            for (String a : p1){
                if(p12.add(a)){
                    bw1.write(a);
                    bw1.newLine();
                }
            }
            bw1.flush();
            bw1.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getDiversityFast(String inFile, String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            Integer[] geno = null;
            int sampleSize = 0;
            StringBuilder geo = new StringBuilder(500000);
            int SNP = 0;
            int posStart = 0;
            int posEnd = 1;
            int geneLength = 1;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            bw.write("chr\tStartPos\tendPos\tPI\n");
            String chr = "";
            boolean start = true;
            boolean end = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    if(temp.startsWith("#C")){
                        sampleSize = temp.split("\t").length-9;
                    }
                    continue;
                }
                
                SNP++;
                geno = new Integer[sampleSize];
                geo = new StringBuilder(500000);
                te = temp.split("\t");
                if(start) chr = te[0];
                
                if(!chr.equals(te[0])){
//                    double pi = estimatePi(posStart, posEnd, sampleSize, genotype,h);
//                    bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
//                    end = true;
                    start = true;
                }
                
                if(posStart ==0) posStart = Integer.parseInt(te[1]);
                posEnd = Integer.parseInt(te[1]);
                for (int i = 0;i < te.length-9; i++){
                    if(te[i+9].startsWith(".")){
                        geno[i] = -1;
                    }else{
                        String[] te1 = te[i+9].split("/");
                        geno[i] = Integer.parseInt(te1[0]) + Integer.parseInt(te1[1].split(":")[0]);
                    }
                }
                
                if(!end && (posEnd-posStart) > 5000){
//                    double pi = estimatePiFast(posStart, posEnd, sampleSize, genotype);
//                    bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
                    
                    posStart = 0;
                    start = true;
                }else if(end){
                    posStart = 0;
                    end = false;
                    start = true;
                }
            }
//            double pi = estimatePi(posStart, posEnd, sampleSize, genotype,h);
//            bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
            bw.flush();
            bw.close();
            System.out.println("Total SNPs detected:\t"+SNP);
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
//    public double estimatePiFast(int posStart, int posEnd,int sampleSize,List<Integer[]> genotype){
//        double pi = 0;
//        int sum = 0;
////            int diff  = 0;
////            double pij = 0;
//        for(int i = 0, length = genotype.size(); i < length-1;i++){
//                for(int j = i+1; j < length -1;j++){
//                    sum = 0;
//                    diff = 0;
//                    for(int z = 0; z < sampleSize;z++){
//                        if(genotype.get(i)[z]>-1 && genotype.get(j)[z]>-1){
//                            sum++;
//                            if(genotype.get(i)[z] != genotype.get(j)[z]){
//                                diff++;
//                            }
////                            System.out.println(diff/sum);
//                        }
//                    }
//                    if(sum == 0){
//                        pij = 0;
//                    }else{
//                        pij += (double)diff/sum * h.get(genotype.get(i))*h.get(genotype.get(j))/geneLength/geneLength; 
//                    }
//                }
//            }
//            pi = pij*2;
//    }
    public void getDiversity(String inFile,String bedFile,String outFile,int windowSize){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            Integer[] geno = null;
            int sampleSize = 0;
            List<Integer[]> genotype = new ArrayList();
            Set ge = new HashSet();
            Map<Integer[],Integer> h = new HashMap();
            StringBuilder geo = new StringBuilder(500000);
            int SNP = 0;
            int posStart = 0;
            int posEnd = 1;
            int geneLength = 1;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            bw.write("chr\tStartPos\tendPos\tPI\n");
            String chr = "";
            boolean start = true;
            boolean end = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    if(temp.startsWith("#C")){
                        sampleSize = temp.split("\t").length-9;
                    }
                    continue;
                }
                
                SNP++;
                geno = new Integer[sampleSize];
                geo = new StringBuilder(500000);
                te = temp.split("\t");
                if(start) chr = te[0];
                
                if(!chr.equals(te[0])){
                    double pi = estimatePi(posStart, posEnd, sampleSize, genotype,h);
                    bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
                    h = new HashMap();
                    genotype = new ArrayList();
                    ge = new HashSet();
                    end = true;
                    start = true;
                }
                
                if(posStart ==0) posStart = Integer.parseInt(te[1]);
                posEnd = Integer.parseInt(te[1]);
                for (int i = 0;i < te.length-9; i++){
                    if(te[i+9].startsWith(".")){
                        geno[i] = -1;
                    }else{
                        String[] te1 = te[i+9].split("/");
                        geno[i] = Integer.parseInt(te1[0]) + Integer.parseInt(te1[1].split(":")[0]);
                    }
                }
                if(ge.add(geno)){
                    h.put(geno,1);
                }else{
                    int t = h.get(geno)+1;
                    h.put(geno,t);
                }
                genotype.add(geno);
                if(!end && (posEnd-posStart) > windowSize){
                    double pi = estimatePi(posStart, posEnd, sampleSize, genotype,h);
                    bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
                    h = new HashMap();
                    genotype = new ArrayList();
                    ge = new HashSet();
                    posStart = 0;
                    start = true;
                }else if(end){
                    posStart = 0;
                    end = false;
                    start = true;
                }
            }
            double pi = estimatePi(posStart, posEnd, sampleSize, genotype,h);
            bw.write(chr+"\t"+posStart+"\t"+posEnd+"\t"+pi+"\n");
            bw.flush();
            bw.close();
            System.out.println("Total SNPs detected:\t"+SNP);
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private double estimatePi(int posStart, int posEnd,int sampleSize,List<Integer[]> genotype,
            Map<Integer[],Integer> h){
            int sum = 0;
            int diff  = 0;
            double pij = 0;
            double pi = 0;
//            System.out.println(genotype.get(0)[0]);
            int geneLength = posEnd - posStart +1;
            for(int i = 0, length = genotype.size(); i < length-1;i++){
                for(int j = i+1; j < length -1;j++){
                    sum = 0;
                    diff = 0;
                    for(int z = 0; z < sampleSize;z++){
                        if(genotype.get(i)[z]>-1 && genotype.get(j)[z]>-1){
                            sum++;
                            if(genotype.get(i)[z] != genotype.get(j)[z]){
                                diff++;
                            }
//                            System.out.println(diff/sum);
                        }
                    }
                    if(sum == 0){
                        pij = 0;
                    }else{
                        pij += (double)diff * h.get(genotype.get(i))*h.get(genotype.get(j))/sampleSize/sampleSize/geneLength; 
                    }
                }
            }
            pi = pij*2;
            return pi;
    }
    public void equalDepth(String inFile,String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            int[] depth = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            List<int[]> depthInfo = new ArrayList();
            int sampleSize = 0;
            boolean start = true;
            int posPre = 0;
            List genotype = new ArrayList();
            boolean pass = false;
            while ((temp = br.readLine())!= null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                    if(temp.startsWith("#C")){
                        sampleSize = temp.split("\t").length-9;
                        depth = new int[sampleSize];
                        for (int i = 0; i<sampleSize;i++){
                            depth[i] = 0;
                        }
                    }
                    continue;
                }
                te = temp.split("\t");
                te[2]= "rs" + te[0]+"_"+te[1];
                StringBuilder geno = new StringBuilder();
                for (int i =0; i < te.length-1;i++){
                    geno.append(te[i]);
                    geno.append("\t");
                }
                geno.append(te[te.length-1]);
                genotype.add(geno.toString());
                if(start){
                    posPre = Integer.parseInt(te[1]);
                    start = false;
                }
                for (int i = 0; i < te.length-9; i++){
                    if(!te[i+9].startsWith(".")){
                        try {
                            depth[i] = Integer.parseInt(te[i+9].split(":")[2]);
                        }catch (Exception e){
                            depth[i] = 0;
                        }
                    }
                }
                depthInfo.add(depth);
                if((Integer.parseInt(te[1]) - posPre) > 20000){
                    pass = checkDepth(depthInfo);
                    start = true;
                    if(pass){
                        writeGeno(genotype,bw);
                    }
                    genotype = new ArrayList();
                    depthInfo = new ArrayList();
                }
            }
            bw.flush();
            bw.close();
         } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    private boolean checkDepth(List<int[]> depthInfo){
        boolean pass = false;
        int sum = 0;
        int sumLine = 0;
        int passed = 0;
        int p = 0;
        int snp = 0;
        for (int[] line:depthInfo){
            sumLine = 0;
            sum++;
            snp = 0;
            for(int i =0; i< line.length;i++){
                sumLine++;
//                if(line[i]>0) sumLine++;
                if(line[i]>1 & line[i] <5 ) p++;
            }
            if(sumLine - p < 5) passed++;
        }
        if ((sum - passed) < 0.01*sum) pass = true;
        return pass;
    }
    private void writeGeno(List<String> geno,BufferedWriter bw){
        for (String ge:geno){
            try {
                bw.write(ge);
                bw.newLine();
            } catch (IOException ex) {
                Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
    public static void diagSNP(String inFile, String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String[] te = null;
            double sum = 0;
            double het = 0;
            double alt = 0;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")) continue;
                te = temp.split("\t");
                sum = 0;
                het = 0;
                alt = 0;
                for (int i =9 ; i < te.length; i++){
                   if(te[i].startsWith("0/0")) {
                       sum+=2;
                   }else if (te[i].startsWith("0/1")){
                       sum+=2;
                       het++;
                       alt++;
                   }else if (te[i].startsWith("1/1")){
                       sum+=2;
                       alt+=2;
                   }
                }
                bw.write(het/sum+"\t"+alt/sum+"\n");
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void changeName(String inFile,String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    te[2] = "rs"+te[0]+"_"+te[1];
                    for (int i = 0; i< te.length-1;i++){
                        bw.write(te[i]+"\t");
                    }
                    bw.write(te[te.length-1]+"\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void mergeABD(String inFile, String outFile){
        try {
            String temp = "";
            String tem = "";
            BufferedReader bf = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            int i = 0;
            int s = 0;
            String[] te = null;
            List name = new ArrayList();
            int[][] index = new int[3][300];
            for (int a = 0 ; a < index.length;a++){
                for (int b = 0; b< index[0].length;b++){
                    index[a][b] = 0;
                }
            }
            int fn = -1;
            
            while ((temp = bf.readLine() )!=null){
                fn++;
                TabixReader br = new TabixReader(temp);
                while ((tem = br.readLine())!=null){
                    i ++;
                    if(tem.startsWith("##") && i < 3){
                        bw.write(tem);
                        bw.newLine();
                    }else if (tem.startsWith("#C")){
                        te = tem.split("\t");
                        for(int j = 0; j < te.length ; j++){
                            int pos = searchArray(name,te[j]);
                            if (pos > -1){
                                index[fn][j] = pos;
                            }else{
                                name.add(te[j]);
                                index[fn][j] = name.size()-1;
                            }
                        }
                        tem = null;
                    }
                }
                br.close();
            }
            bf.close();
            for(int j = 0;j<name.size()-1;j++){
                bw.write(name.get(j)+"\t");
            }
            bw.write(name.get(name.size()-1)+"\n");
            fn = -1;
            String[] t = new String[name.size()];
            bf = IOUtils.getTextReader(inFile);
            while ((temp = bf.readLine())!=null){
                t = initialString(new String[name.size()]);
                fn++;
                TabixReader br = new TabixReader(temp);
                while ((tem = br.readLine())!=null){
                    if(tem.startsWith("#")){
                        continue;
                    }
                    te = tem.split("\t");
                    for(int j =0; j< te.length;j++){
                        t[index[fn][j]] = te[j];
                    }
                    for(int j = 0; j < t.length -1; j++){
                        bw.write(t[j]+"\t");
                    }
                    bw.write(t[t.length-1]+"\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static String[] initialString (String[] s){
        String[] r = new String[s.length];
        for (int i = 0 ; i < s.length; i++){
            r[i]  = "./.";
        }
        return r;
    }
    public static int searchArray(List A, String name){
        int index = -1;
        for (int i = 0; i < A.size();i++){
            if(A.get(i).equals(name)){
                index = i;
                break;
            }
        }
        return index;
    }
    public static void toFasta(String inFile, String outFile){
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inFile);
            if(!inFile.endsWith(".gz")) br = IOUtils.getTextReader(inFile);
            StringBuilder[] fa = null;
            String[] sample = null;
            String temp = "";
            int SN = 0;
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("##")){
                    continue;
                }else if (temp.startsWith("#")){
                    te = temp.split("\t");
                    SN = te.length - 9;
                    sample = new String[SN];
                    fa = new StringBuilder[SN];
                    for(int i = 0; i < SN; i++){
                        sample[i] = te[i+9];
                        fa[i] = new StringBuilder("");
                    }
                }else{
                    te = temp.split("\t");
                    for (int i = 0; i < SN; i++){
                        if(te[i+9].startsWith("0/0")){
                            fa[i].append(te[3]);
                        }else if (te[i+9].startsWith("1/1")){
                            fa[i].append(te[4]);
                        }else {
                            fa[i].append("-");
                        }
                    }
                }
            }
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            for (int i = 0 ; i < SN;i++){
                bw.write(">"+sample[i]);
                bw.newLine();
                bw.write(fa[i].toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void changeVCF(String inFile, String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null , p1 = null;
            int DP = 0;
            String chr = "";
            Integer pos = 0;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    p1 = te[0].split(":");
                    chr = p1[0];
                    pos = Integer.parseInt(p1[1].split("-")[0]) 
                            + Integer.parseInt(te[1]);
                    te[0] = chr;
                    te[1] = Integer.toString(pos);
                    for (int i = 0; i < te.length - 1; i++){
                        bw.write(te[i]+"\t");
                    }
                    bw.write(te[te.length - 1]);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void filtTotalDepth(String inFile, String outFile, double dp){
        try {
            TabixReader br = new TabixReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            int DP = 0;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    DP = Integer.parseInt(temp.split("DP=")[1].split(";")[0]);
                    if(DP > dp) continue;
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
    
    public static void getEhhMap(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp ="";
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("ID")) continue;
                te = temp.split(" ");
                bw.write(te[0].split("-")[0] + " " + te[0] + " " + Double.parseDouble(te[1])/100000000 + " " + te[1]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void phase(String inFile,String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String te = "";
//            String map = outFile.replace(".vcf",".map");
//            BufferedWriter bm = IOUtils.getTextWriter(map);
            System.out.print("Analyzing: "+inFile);
            String[] tem = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.replaceAll("/", "|");
                    tem = temp.split("\t");
//                    bm.write(tem[0]+" "+ tem[1]+" "+Double.parseDouble(tem[2])/100000000+" " +tem[3]);
                    bw.newLine();
                    bw.write(te);
                    bw.newLine();
                }
            }
//            bm.flush();
//            bm.close();
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public void genotypeStat(String inFile, String outFile){
        try{
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            int[][] res = null;
            String[] name = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("##")) continue;
                if(temp.startsWith("#")){
                    te = temp.split("\t");
                    res = new int[te.length-9][4];
                    name = new String[te.length-9];
                    for (int i = 0; i < te.length-9;i++){
                        name[i] = te[i+9];
                        for (int j = 0; j< 4; j ++){
                            res[i][j] = 0;
                        }
                    }
                }else{
                    te = temp.split("\t");
                    for(int i = 0; i<te.length-9;i++){
                        if(te[i+9].startsWith("0/0")){
                            res[i][0]++;
                        }else if(te[i+9].startsWith("1/1")){
                            res[i][1]++;
                        }else if(te[i+9].startsWith("0/1")){
                            res[i][2]++;
                        }else {
                            res[i][3]++;
                        }
                    }
                }
            }
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            bw.write("Hom(REF)\tHom(ALT)\tHet\tMissing\n");
            for (int i = 0;i<te.length-9;i++){
                bw.write(name[i]+"\t"+res[i][0]+"\t"+res[i][1]+"\t"+res[i][2]+"\t"+res[i][3]+"\n");
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            
        }
    }
    public void checkChr(String inFile){
        try {
            TabixReader br = new TabixReader(inFile);
            Set chr = new HashSet();
            String temp = "";
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    if(chr.add(te[0])){
                        System.out.println("Chomosome:\t"+te[0]);
                    }
                }  
            }
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static void generateNULL(String inFile, String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void mergeAll(String inFile,String outFile,String suffix){
        File inF = new File(inFile);
        File[] ins = IOUtils.listRecursiveFiles(inF);
        File[] subin = IOUtils.listFilesEndsWith(ins, suffix);
        boolean start = true;
        BufferedReader br = null;
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        try {
            for (File ff : subin){

                br = IOUtils.getTextReader(ff.getAbsolutePath().toString());
                String temp = "";
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("#")){
                        if(start){
                            bw.write(temp);
                            bw.newLine();
                            continue;
                        }else{
                            continue;
                        }
                    }
                    bw.write(temp);
                    bw.newLine();
                }
                start = false;

            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
                Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    public void mergeFiles(String inFile,String inFile2,String outFile){
        boolean start = true;
        BufferedReader br = null;
        BufferedReader br2 = IOUtils.getTextReader(inFile2);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String ff = "";
        int i = 0;
        try {
            while((ff = br2.readLine())!=null){
                br = IOUtils.getTextGzipReader(ff);
                String temp = "";
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("#")){
                        if(start){
                            bw.write(temp);
                            bw.newLine();
                            continue;
                        }else{
                            continue;
                        }
                    }
                    i++;
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();
                start = false;
            }
            System.out.println("i is: " + i);
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    
    public static void merge(String inFile,String inFile2, String outFile){
        try {
            TabixReader ref = new TabixReader(inFile);
            TabixReader alt = new TabixReader(inFile2);
            String temp1 = "";
            String temp2 = "";
            String[] te1 = null;
            String[] te2 = null;
            Set pos = new HashSet();
            StringBuilder head = new StringBuilder();
            int snp = 0;
            int read = 0;
            while ((temp2 = alt.readLine())!=null){
                if(temp2.startsWith("##")){
                    head.append(temp2);
                    head.append("\n");
                }else{
                    read ++;
                    if (read % 100000 == 0){
                        System.out.println("vcf file SNP readed:\t"+read/1000+"K...");
                    }
                    te2 = temp2.split("\t");
                    pos.add(te2[1]);
                }
            }
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            BufferedWriter bw1 = IOUtils.getTextWriter(outFile+".log");
            
            System.out.println("VCF file "+ read +" SNPs readed! Now reading gvcf file...");
            
            read = 0;
            while ((temp1 = ref.readLine())!=null){
                if(temp1.startsWith("##")) continue;
                if(temp1.startsWith("#")){
                    head.append(temp1);
                    bw.write(head.toString());
                    bw.newLine();
                    continue;
                }
                read ++;
                if (read % 100000 == 0){
                    System.out.println("gvcf file SNP readed:\t"+read/1000+"K...");
                }
                te1 = temp1.split("\t");
                if(pos.add(te1[1])) continue;
                snp++;
                bw.write(temp1);
                bw.newLine();
            }
            bw1.write("Total SNP is:\t"+snp);
            bw1.flush();
            bw1.close();
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public VcfTools(String inFile,String outFile,Integer windowSize){
        vcf a = new vcf();
        a.initialVCFread(inFile);
//        List<String[]> res = new ArrayList(a.sampleSize + 3);
        String[] temp = new String[a.sampleSize + 3];
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        try{
        while(!a.checkEnd()){
            a.readVCFBlock(windowSize);
            temp[0] = a.getChr();
            temp[1] = a.getStartPos();
            temp[2] = a.getendPos();
            int[][]aa =  blockHet(a);
            for(int i = 0;i<a.sampleSize;i++){
                temp[i+3] = Double.toString(aa[0][i]/((double)aa[1][i]));
            }
            writeString(temp,bw);
        }
        bw.flush();
        bw.close();
        } catch (Exception e){
                    
        }
        
       }
    private void writeString(String[] a,BufferedWriter bw){
        try {
            for(int i = 0; i< a.length-1; i++){
                bw.write(a[i]);
                bw.write("\t");
            }
            bw.write(a[a.length-1]);
            bw.newLine();
        } catch (IOException ex) {
                Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
            }
    }
    public VcfTools(String anchor, String inFile, String outFile, int minComp, double maxIBDDist, 
        int windowSize, int numThreads, int bestContrasts){
        ibdfilter.PfilterBasedOnIBD(anchor, inFile, outFile, minComp, maxIBDDist, windowSize, numThreads, bestContrasts);
    }
    public VcfTools(String inFile, String outFile,int windowSize, double threshold){   
        LD.calLD(inFile,outFile,windowSize,threshold);
    }
    public VcfTools(String inFile,String outFile,double se,String suffix){
        this.rsdFilter(inFile,outFile,se,suffix);
    }
   
    // type : filtered
    // parameters: MQ, FS, MQRankSum, and ReadPosRankSum
    public VcfTools(String inFile,String outFile,String MQ,String FS,String MQRankSum,String ReadPosRankSum,String BSQRankSum,String SOR){
        this.getFilterd(inFile,outFile,MQ,FS,MQRankSum, ReadPosRankSum,BSQRankSum,SOR);
    }
    private void phaseVCF(String inFile,String outFile){
        try {
            TabixReader in = new TabixReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            while((temp = in.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    for (int i = 0; i < 9 ; i++){
                        bw.write(te[i]+"\t");
                    }
                    for(int i = 0;i < te.length-10;i++){
                        if(te[i+9].startsWith("0/0")){
                            bw.write("0|0\t");
                        }else if (te[i+9].startsWith("1/1")){
                            bw.write("1|1\t");
                        }else if(te[i+9].startsWith("0/1")){
                            bw.write("0|1\t");
                        }else{
                            bw.write("?|?\t");
                        }
                    }
                    if(te[te.length-1].startsWith("0/0")){
                        bw.write("0|0\t");
                    }else if (te[te.length-1].startsWith("1/1")){
                        bw.write("1|1\t");
                    }else if(te[te.length-1].startsWith("0/1")){
                        bw.write("0|1\t");
                    }else{
                        bw.write("?|?\t");
                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
//    public VcfTools(String inFile,String outFile,String model,int size,int SNPnum,int header){
//        this.getSub(inFile,outFile,size,SNPnum);
//    }
    public VcfTools(String inFile,String outFile, double a, double b, double sd,
        int mindepth,int maxdepth){
        this.getDepthFilterd(inFile,outFile,a,b,sd,mindepth,maxdepth);
    }
    public VcfTools(String inFile, String suffix,String outsuffix,
            String MQ,String FS,String MQRankSum,String ReadPosRankSum,String BSQRankSum,
            String SOR,int maxdepth, int mindepth, double maxSD, int windowSize,
            double threshold){
        try{
            String anchor = inFile.replace(suffix,".anchor.vcf");
            BufferedWriter bw = IOUtils.getTextWriter(inFile.replace(".vcf",".txt"));
            String outFile = inFile.replace(".vcf",".quality"+MQ+".vcf");
            this.getFilterd(inFile,outFile,MQ,FS,MQRankSum, ReadPosRankSum,BSQRankSum,SOR);
            bw.write(SNPsAll+"\t"+SNPsFail);
            bw.newLine();
            bw.flush();
            inFile = outFile;
            outFile = inFile.replace(".quality"+MQ+".vcf",".depth"+MQ+".vcf");
            this.getDepthFilterd(inFile,outFile,0,0,maxSD,mindepth,maxdepth);
            bw.write(SNPsAll+"\t"+SNPsFail);
            bw.newLine();
            bw.flush();
            inFile = outFile;
            outFile = inFile.replace(".depth"+MQ+".vcf",".ST"+MQ+".vcf");
            Integer[] SNPnum = new SegeregationTest().getST(inFile,outFile);
            bw.write(SNPnum[0]+"\t"+SNPnum[1]);
            bw.newLine();
            bw.flush();
            inFile = outFile;
            outFile = inFile.replace(".ST"+MQ+".vcf",".LD"+MQ+".vcf");
            SNPnum = LD.calLD(inFile,outFile,windowSize,0.1);
            bw.write(SNPnum[0]+"\t"+SNPnum[1]);
            bw.newLine();
            bw.flush();
            inFile = outFile;
            outFile = inFile.replace(".LD"+MQ+".vcf",".IBD"+MQ+".vcf");
            
            ibdfilter.PfilterBasedOnIBD(anchor, inFile, outFile, 200, 0.03, windowSize, 4, -1);
            bw.close();
        }catch(Exception e){
            
        }
    }
    public static void mergeVCF(String inFile, String outFile){
        BufferedReader br1;
        BufferedWriter bw;
//        StringBuilder files = new StringBuilder();
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
        int i = 0;
        try{
            TabixReader br = null;
            bw = IOUtils.getTextWriter(outFile);
            String temp = "", file=null;
            String[] te = null;
            br1 = IOUtils.getTextReader(inFile);
            boolean head = true;
            while ((file = br1.readLine())!=null){
                br = new TabixReader(file);
                System.out.println("Reading "+ file);
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("#")){
                        if(head) {
                            bw.write(temp);
                            bw.newLine();
                        }
                    }else{
                       i++;
                       head = false;
                       bw.write(temp); 
                       bw.newLine();
                       bw.flush();
                    }
                }
                br.close();
            }
            System.out.println("Total SNPs is : "+i);
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
    
    public static void subChr(String inFile,String chr,String outFile){
        try {
            TabixReader br = new TabixReader(inFile) ;
            
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            Set chrs = getSplit(chr);
            String temp =null;
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    if(chrs.add(te[0])){
                        chrs.remove(te[0]);
                    }else{
                        bw.write(temp);
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            
        }
    }
    public static void statHet(String inFile,String outFile,int window){
        try {
            TabixReader br = new TabixReader(inFile) ;
            
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private static Set getSplit(String chr){
        String[] chrs = chr.split(",");
        Set a = new HashSet();
        for (int i=0;i< chrs.length;i++){
            a.add(chrs[i]);
        }
        return a;
    }
    public static void splitByChr(String inFile,String outFile,vcf a){
        a.initialVCFread(inFile);
        File byChr = new File(outFile+"/ByChr");
        if(!byChr.exists()) byChr.mkdirs();
        while(!a.checkEnd()){
            a.readVCFByChr();
            String outFiles = outFile+"/ByChr/"+a.chrom+".vcf";
            a.initialVCFwrite(outFiles);
            a.writeVCF();
            a.close();
        }
    }
    public static void splitByChrs(String inFile,String outFile){
        File input = new File(inFile);
        File[] ins = IOUtils.listRecursiveFiles(input);
        File[] subFile = IOUtils.listFilesEndsWith(ins, ".gz");
        for (File fs : subFile){
            String in = fs.getAbsolutePath().toString();
            System.out.println("Now analyzing: "+in);
            String[] names = in.split("/");
//            System.out.println(in);
            String ind = names[names.length-1].split("\\.")[0];
            String out = outFile+"/"+ind;
            splitByChr(in,out);
        }
    }
        
    private void getFilterd(String inFile,String outFile,String MQ, String FS, String MQRankSum
    , String ReadPosRankSum , String BSQRankSum,String SOR){
        try {
            System.out.printf("Filtering by FS > %s; MQ < %s; MQRankSum < %s; ReadPosRankSum < %s;"
                    + " BSQRankSum < %s; SOR > %s ...\n",FS,MQ,MQRankSum,ReadPosRankSum,
                  BSQRankSum,SOR);
            TabixReader vcf = new TabixReader(inFile);
//            if(inFile.endsWith("gz")){
//                vcf = IOUtils.getTextGzipReader(inFile);
////                outFile = inFile.replace(".vcf.gz",".filtered.vcf");
//            }
//            else {
//                vcf = IOUtils.getTextReader(inFile);
////                outFile = inFile.replace(".vcf",".filtered.vcf");
//            }
            String mq = "70", fs = "40", mqranksum = "-12.5", readposranksum = "-8",
                    bqrs = "0",sor="1";
            String temp = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            int snps = 0, rmsnp = 0;
            boolean filter = true;
            while((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snps++;
                    if(snps % 100000 == 0) System.out.println("Analyzing SNPs\t"+snps+"...........");
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
                        sor = temp.split("SOR=")[1].split("\t")[0];}
                    catch (Exception ex){
                    };
                    if(Double.parseDouble(sor) > Double.parseDouble(SOR)) continue;
                    
                    try{
                        bqrs = temp.split("BaseQRankSum=")[1].split(";")[0];}
                    catch (Exception ex){
                    };
                    if(Double.parseDouble(bqrs) < Double.parseDouble(BSQRankSum)) continue;
                    rmsnp++;
                    bw.write(temp);
                    bw.newLine();
                    bw.flush();
                }else if(temp.startsWith("##FILTER")){
                    if(filter){
                        filter = false;
                        bw.write("##FILTER=<ID=LowMQ,Description=\"MQ < "+MQ+"\">\n");
                        bw.write("##FILTER=<ID=HighFS,Description=\"FS > "+FS+"\">\n");
                        bw.write("##FILTER=<ID=HighSOR,Description=\"SOR < "+SOR+"\">\n");
                        bw.write("##FILTER=<ID=LowDP,Description=\"DP < "+2+"\">\n");
                        bw.write("##FILTER=<ID=LowQD,Description=\"QD < "+2.0+"\">\n");
                        bw.write("##FILTER=<ID=LowBaseQRankSum,Description=\"BaseQRankSum < "+bqrs+"\">\n");
                        bw.write("##FILTER=<ID=LowMQRankSum,Description=\"MQRankSum < "+MQRankSum+"\">\n");
                        bw.write("##FILTER=<ID=LowReadPosRankSum,Description=\"ReadPosRankSum < "+ReadPosRankSum+"\">\n");
                    }
                }else{
                    bw.write(temp);
                    bw.newLine();
                    bw.flush();
                }
            }
            SNPsAll = snps;
            SNPsFail = snps-rmsnp;
            System.out.println("Total SNPs:\t"+snps);
            System.out.println("SNPs removed:\t"+(snps-rmsnp));
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    private void getXPCLR(String inFile,String outFiles){
        try {
            TabixReader br = new TabixReader(inFile);
            String outFile = outFiles + ".geno";
            String map = outFiles + ".snp";
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            BufferedWriter bsnp = IOUtils.getTextWriter(map);
            String temp = null;
            int snp = 0;
            while((temp = br.readLine())!=null){
                if(!temp.startsWith("#")){
                    snp++;
                    if(snp%1000000==0){
                        System.out.println("SNPs\t"+snp/1000+"K...");
                    }
                    String[] tmp = temp.split("\t");
                    bsnp.write(" rs"+tmp[0]+"_"+tmp[1]+"\t"+tmp[0]+"\t"+Double.valueOf(Integer.parseInt(tmp[1]))/100000000+"\t"+tmp[1]+"\t"+tmp[3]+"\t"+tmp[4]+"\n");
                    
                    for (int i = 0; i<tmp.length-10;i++){
                        if(tmp[i+9].startsWith("0/0")){
                            bw.write("0 0 ");
                        }else if (tmp[i+9].startsWith("1/1")){
                            bw.write("0 1 ");
                        }else if (tmp[i+9].startsWith("0/1")){
                            bw.write("1 1 ");
                        }else{
                            bw.write("9 9 ");
                        }
                    }
                    if(tmp[tmp.length-1].startsWith("0/0")){
                        bw.write("0 0\n");
                    }else if (tmp[tmp.length-1].startsWith("1/1")){
                        bw.write("0 1\n");
                    }else if (tmp[tmp.length-1].startsWith("0/1")|tmp[tmp.length-1].startsWith("0/2")){
                        bw.write("1 1\n");
                    }else if(tmp[tmp.length-1].startsWith("2/2")){
                        bw.write("1 1\n");
                    }else {
                        bw.write("9 9\n");
                    }
                }
            }
            bw.flush();
            bw.close();
            bsnp.flush();
            bsnp.close();
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
                outFile = inFile.replace(".vcf.gz",".mq");
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
                            nvcf = IOUtils.getTextWriter(theDir.toString()+"/"+end+".chr"+chr+".temp"+size+".vcf");
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
                        nvcf = IOUtils.getTextWriter(theDir.toString()+"/"+end+".chr"+chr+".temp"+size+".vcf");
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
    private static void splitByChr(String inFile, String outFile) {
        try {
            File outf = new File(outFile);
            if(!outf.exists()) outf.mkdirs();
            TabixReader vcf = new TabixReader(inFile);
            String temp = "";
            String chr = "";
            BufferedWriter bw = null ;
            StringBuilder head = new StringBuilder();
            while((temp = vcf.readLine())!=null){
                if(temp.startsWith("#")){
                    head.append(temp);
                    head.append("\n");
                }else{
                    String[] temps = temp.split("\t");
                    if(temps[0].equals(chr)){
                        bw.write(temp);
                        bw.newLine();
                    }else{
                        if(bw!=null){
                            bw.flush();
                            bw.close(); 
                        }
                        System.out.println("Chromosome " + chr + " finished!");
                        chr = temps[0];
                        String out = outFile+"/"+chr+".vcf";
                        bw = IOUtils.getTextWriter(out);
                        System.out.println("Chrmosome: " + chr + " Starting");
                        bw.write(head.toString());
                        bw.write(temp);
                        bw.newLine();
                    }
                }
            }
            System.out.println("Chromosome " + chr + " finished!");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null,ex);
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
//            outFile = inFile.replace(".vcf",".depth.vcf");
            VcfDepth = IOUtils.getTextWriter(outFile);
            String temp = null;
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            Double d = 0.0;
            int snps = 0,rmsnps = 0;
            int sampleNum = 0;
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snps++;
                    if(snps%1000000 == 0) System.out.println("Anlyzing " + snps +"...");
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    int dp =0;
                    double rsd = 0;
                    int s=0;
                    if(tem[4].length()==1){
//                        dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
                        double dep[] = new double[sampleNum];
                        for (int i = 9; i < tem.length; i++){
                            String re = tem[i];
                            if(re.split(":").length<3){
//                                rmsnps++;
                                continue;
                            }
                            if(!re.split(":")[0].startsWith(".")){
                                d = Double.parseDouble(re.split(":")[2]);
                                s++;
                            }else{
                                d = Double.NaN;
                            }
                            dep[i-9] = d;
                        }
                        for(int i = 0; i< dep.length;i++){
                            if(!Double.isNaN(dep[i])) dp+=dep[i];
                        }
                        if(dp < mindepth ){
                            rmsnps++;
                            continue;
                        }
                        if (dp > maxdepth) {
                            rmsnps++;
                            continue;
                        }
//                        StandardDeviation rrsd = new StandardDeviation();
                        rsd = getSD(dep);
//                        double rsd1 = getSD(dep);
                        if(rsd > sd){
                            rmsnps++;
                            continue;
                        }
//                        if(dp > meanCUT){
                        VcfDepth.write(temp);
                        VcfDepth.newLine();
                        VcfDepth.flush();
                    }
                }else{
                    VcfDepth.write(temp);
                    VcfDepth.newLine();
                    VcfDepth.flush();
                }
            }
            SNPsAll = snps;
            SNPsFail = rmsnps;
            VcfDepth.flush();
            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    private void dpFilter(String inFile,String outFile,double se) {
        try {
            System.out.println("Now Filtering SNP based on depth...");
            BufferedReader vcf;
            BufferedWriter VcfDepth;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
//            outFile = inFile.replace(".vcf",".depth.vcf");
            VcfDepth = IOUtils.getTextWriter(outFile);
            String temp = null;
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            int snps = 0,rmsnps = 0;
            int sampleNum = 0;
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snps++;
                    if(snps%1000000 == 0) System.out.println("Anlyzing " + snps +"...");
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    int dp =0;
                    double d = 0;
                    double rsd = 0;
                    int s=0;
                    if(tem[4].length()==1){
//                        dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
                        double dep[] = new double[sampleNum];
                        for (int i = 9; i < tem.length; i++){
                            String re = tem[i];
                            if(re.split(":").length<3){
//                                rmsnps++;
                                continue;
                            }
                            if(!re.split(":")[0].startsWith(".")){
                                d = Double.parseDouble(re.split(":")[2]);
                                dp += d;
                                s++;
                            }else{
                                d = Double.NaN;
                            }
                            dep[i-9] = d;
                        }
//                        if(dp < s * 1.8 | dp > s* 5.5 ) continue;
//                        StandardDeviation rrsd = new StandardDeviation();
                        double rs = getSD(dep);
                        double rsd1 = rs/dp;
                        if(rsd1 > se) {
                            continue;
                        }
//                        if(dp > meanCUT){
                        VcfDepth.write(temp);
                        VcfDepth.newLine();
                        VcfDepth.flush();
                    }
                }else{
                    VcfDepth.write(temp);
                    VcfDepth.newLine();
                    VcfDepth.flush();
                }
            }
            SNPsAll = snps;
            SNPsFail = rmsnps;
            VcfDepth.flush();
            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    public  void rsdFilter(String inFile,String outFile,double se,String suffix) {
        try {
            System.out.println("Now Filtering SNP based on rsd...");
            File f = new File(inFile);
            File[] fs = IOUtils.listRecursiveFiles(f);
            File[] sub = IOUtils.listFilesEndsWith(fs, suffix);
            BufferedReader vcf;
            BufferedWriter VcfDepth;
            VcfDepth = IOUtils.getTextWriter(outFile);
            for (File fi : sub){
                vcf = IOUtils.getTextReader(fi.toString());
                String out = fi.toString().replace(".vcf",".dp.vcf");
                BufferedWriter bw = IOUtils.getTextWriter(out);
                VcfDepth.write(fi.toString()+"\t");
                String temp = null;
                String[] tem = null;
                int snps = 0,rmsnps = 0;
                int sampleNum = 0;
                while ((temp = vcf.readLine())!=null){
                    if(!temp.startsWith("#")){
                        snps++;
                        if(snps%1000000 == 0) System.out.println("Anlyzing " + snps +"...");
                        tem = temp.split("\t");
                        sampleNum = tem.length - 9;
                        int dp =0;
                        double d = 0;
                        double rsd = 0;
                        int s=0;
                        if(tem[4].length()==1){
//                        dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
                            double dep[] = new double[sampleNum];
                            for (int i = 9; i < tem.length; i++){
                                String re = tem[i];
                                if(re.split(":").length<3){
                                    continue;
                                }
                                if(!re.split(":")[0].startsWith(".")){
                                    d = Double.parseDouble(re.split(":")[2]);
                                    dp += d;
                                    s++;
                                }else{
                                    d = Double.NaN;
                                }
                                dep[i-9] = d;
                            }
                            double rs = getSD(dep);
                            double rsd1 = rs/dp;
                            if(rsd1 > se){
                                rmsnps++;
                                continue;
                            }
                            bw.write(temp);
                            bw.newLine();
                        }
                    }else{
                        bw.write(temp);
                        bw.newLine();
                    }
                    
                }
                bw.flush();
                bw.close();
                VcfDepth.write(snps+"\t"+rmsnps+"\n");
                SNPsAll += snps;
                SNPsFail += rmsnps;
            }
            VcfDepth.write("Total"+SNPsAll+"\t"+SNPsFail+"\n");
            VcfDepth.flush();
            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    private double getSD(double dep[]){
        double sd=0;
        double sum = 0;
        int s = 0;
        for(int i = 0; i< dep.length;i++){
            if(!Double.isNaN(dep[i])){
               sum+=dep[i]; 
               s++;
            }
            
        }
        double mean = sum/s;
        sum =0;
        for(int i = 0; i< dep.length;i++){
            if(!Double.isNaN(dep[i])){
                sum += ((dep[i]-mean)*(dep[i]-mean));
            }
        }
        sd = Math.sqrt(sum/(s-1));
        return sd;
    }
    private void getHet(String inFile,String outFile){
        try {
            System.out.println("Now analyzing the heterozyous....");
            TabixReader vcf = new TabixReader(inFile);
//            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
//            else  vcf = IOUtils.getTextReader(inFile);
//            String hetFile = inFile+".het";
//            if(!outFile.endsWith(".vcf")) outFile = outFile+".het.vcf";
            BufferedWriter vcfhet = IOUtils.getTextWriter(outFile);
//            BufferedWriter het = IOUtils.getTextWriter(hetFile);
//            BufferedWriter het1 = IOUtils.getTextWriter(inFile+".perSNP.het");
            String temp = null;
            String te[] = null;
           
            while ((temp=vcf.readLine())!=null){
                int a0 = 0, a1 = 0, a2 = 0, an = 0;
                int num = 0;
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    double hetrate = 0;
                    double hetNum = 0;
                    if(te[4].length() >1){
                        continue;
                    }
                    for (int i = 9; i < te.length;i++){
                        if(!te[i].startsWith(".")){
                            num++;
                            if(te[i].startsWith("0/1") || te[i].startsWith("1/0")) {
                                hetNum++;
                            }
                        }
                    }
                    hetrate = hetNum/num;
                    vcfhet.write(hetrate+"\n");
                }
            }
            
            vcfhet.flush();
            vcfhet.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private static int[][] blockHet(vcf inVCF){
        // string: het site rate per individual
        int[][] res = new int[2][inVCF.sampleSize];
        for (int i = 0; i< res[0].length;i++){
            res[1][i] = 0;
            res[0][i] = 0;
        }
//        String temp = null;
//        String tem[] = null;
        for(int i = 0;i < inVCF.genotype.size(); i++){
           res = getHetres(inVCF.genotype.get(i),res);
        }     
        
        return res;
    }
    private static int[][] getHetres(String[] geno, int[][] res){
        int[][] res0 = new int[2][res[0].length];
        for (int i = 0; i < res[0].length; i++){
            if(geno[i].startsWith("0/1")){
                res[0][i]++;
                res[1][i]++;
            }else if(geno[i].startsWith(".")){
                
            }else{
                res[1][i]++;
            }
        }
        return res;
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
    public static void getSub(String inFile,String outFile,int size,int SNPnum){
        
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
                po.put(a,0);
                if(pos.size()==SNPnum) break;
            }
//            StringBuilder headerS = new StringBuilder();
//            boolean w = true;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#") ){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    i++;
                    if(i%1000000==0){
                        System.out.println("Analyzing snps: "+i+" ....");
                    }
                    if (po.get(i)!=null){
                        bw.write(temp);
                        bw.newLine(); 
                    }
                    
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
    private void getDepthAll(String inFile,String outFile,double r){
        try {
            System.out.println("Now Testing SNP depth...");
            TabixReader vcf = new TabixReader(inFile);
            BufferedWriter VcfDepth = IOUtils.getTextWriter(outFile);
            
            String temp = "";
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            int snp = 0;
            int sampleNum = 0;
//            VcfDepth.write("chr\tposition\tdepth\tsd");
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snp++;
                    if(snp%100000 == 0) System.out.println("Anlyzing " + snp/1000 +"K...");
                    if(Math.random()>r) continue;
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    double dp =0;
                    double rsd = 0;
                    StringBuilder depth = new StringBuilder();
                    String geno = null;
                    double[] deptheach = new double[sampleNum];
                    int s = 0;
                    if(tem[4].length()==1){
//                        depth.append(Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]));
                        for (int ch = 0 ;ch < sampleNum; ch ++){
                            if(tem[ch+9].startsWith(".")) {
                                deptheach[ch] = Double.NaN;
                            }else{
                                geno = tem[ch+9].split(":")[2];
                                deptheach[ch] = Double.parseDouble(geno);
                                dp += Double.parseDouble(geno);
                                s++;
                            }
                        }
                        VcfDepth.write(tem[0]+"\t"+tem[1]+"\t"+dp+"\t"+s+"\t"+getSD(deptheach));
                        VcfDepth.newLine();
                    }
                }
            }
            StringBuilder DM = new StringBuilder();

            VcfDepth.write(DM.toString());
            VcfDepth.flush();
            VcfDepth.close();
       } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    public void macFilter(String inFile, String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while((temp = br.readLine())!=null){
                int hetNum = 0;
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    for (int i = 0; i<te.length-9;i++){
                        if(te[i].startsWith("0/1")) hetNum++;
                    }
                    if(hetNum > 1){
                        bw.write(temp);
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public void getRandom(String inFile, String outFile, double rate){
        try {
            TabixReader br = new TabixReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    if (Math.random()<rate){
                        bw.write(temp);
                        bw.newLine();
                    };
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(VcfTools.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public void changeHeader(String inFile,String outFile,String suffix){
        try {
            BufferedReader  br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#CHRO")){
                    te = temp.split("\t");
                    te[te.length-1] = suffix;
                    bw.write(te[0]);
                    for (int i = 1; i < te.length; i++){
                        bw.write("\t" + te[i]);
                    }
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
}
