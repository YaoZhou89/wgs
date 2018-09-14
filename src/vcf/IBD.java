/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.BufferedWriter;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

/**
 *
 * @author yaozhou
 */
public class IBD {
    double meanDist=0;
//    int[][] anchorGenotype;
//    int[][] vcfGenotype;
    int contrastNum;
    public IBD(){
        
    }
    /*
    anchor: anchor File
    vcf: SNPs file 
    outFile: output name
    minComp: minimum for comparsion
    maxIBDdist: distance to take as IBD
    windowSize: 
    bestContrast:
    */
    public void ibdFilter(String anchorFile, String vcfFile,String outFile,int minComp, double maxIBDDist, 
        int windowSize){
        vcf a = new vcf();
        vcf v = new vcf();
        a.initialVCFread(anchorFile);
        a.readVCFByChr();
        System.out.println("anchor readed!");
        v.initialVCFread(vcfFile);
        v.readVCFByChr();
         System.out.println("vcf readed!");
        int sampleSize = a.sampleSize, siteSizeA = a.snpNum, startA = a.pos.get(0);
        int siteSizeV = v.snpNum, startV = a.pos.get(0);
        v.initialVCFwrite(outFile);
        Integer siteEndA = a.pos.get(siteSizeA-1);
//        double[] geneticDistance = new double[sampleSize*(sampleSize-1)/2];
        int anchorStart = 0, anchorEnd = windowSize;
        BufferedWriter binfo = IOUtils.getTextWriter(outFile.replace(".vcf",".txt"));
        try {
            binfo.write("Start\tend\tmeanDist\tcontrast_number\n");
        } catch (Exception ex) {
            
        }
        for(int i = 0; i < siteSizeV; i++){
            if(v.pos.get(i) < startA||v.pos.get(i) > siteEndA) {
                v.writeVCF(i);
                continue;
            }
//            if((siteSizeA - anchorEnd) < windowSize) anchorEnd = siteSizeA - 1;
            if(anchorEnd > siteSizeA -1) anchorEnd = siteSizeA-1;
            int siteV = v.posMap.get(a.pos.get(anchorEnd));
            boolean[] pass = new boolean[sampleSize*(sampleSize-1)/2];
            int[][] anchorGenotype = getGenotype(a.genotype,anchorStart,anchorEnd);
            pass = getGeneticDistance(anchorGenotype,minComp,maxIBDDist);
           
            try{
                binfo.write(a.pos.get(anchorStart)+"\t");
                binfo.write(a.pos.get(anchorEnd)+"\t");
                binfo.write(meanDist+"\t");
                binfo.write(contrastNum+"\n");
            }catch (Exception e){
                
            }
            
            int[][] vcfGenotype = getGenotypeTrans(v.genotype,i,siteV);
//            gint[] test = getCol(vcfGenotype,1);
            while(i < siteV){
                int m = 0;
                int[] count = new int[2];
                count[0] = 0; count[1] = 0;
                for(int s1 = 0; s1< sampleSize-1;s1++){
                    int ir = sampleSize/25;
                    for(int s2 = s1+1; s2 < sampleSize;s2++){
                        if(pass[m]){
                            int[] counts = getCount(vcfGenotype[s1],vcfGenotype[s2]);
//                            int[] counts = getCount(test,test);
                            count[0] += counts[0];
                            count[1] += counts[1];
                        }
                        m++;
                    }
                }
                
                if(count[1]> 2* count[0]){
                    v.writeVCF(i);
                }else if((count[0]+count[1])==0){
                    v.writeVCF(i);
                };
                i++;
                if(i%100000 ==0){
                    System.out.println(i+"\tSNPs procceeded!");
                }
            }
            anchorStart = anchorEnd;
            anchorEnd +=windowSize;
        }
        try{
            binfo.flush();
            binfo.close();
            v.close();
        } catch (Exception e){
            
        }
    }
    public void ibdFilterparrallel(String anchorFile, String vcfFile,String outFile,int minComp, double maxIBDDist, 
        int windowSize){
        vcf a = new vcf();
        vcf v = new vcf();
        a.initialVCFread(anchorFile);
        a.readVCFByChr();
        System.out.println("anchor readed!");
        v.initialVCFread(vcfFile);
        v.readVCFByChr();
         System.out.println("vcf readed!");
        int sampleSize = a.sampleSize, siteSizeA = a.snpNum, startA = a.pos.get(0);
        int siteSizeV = v.snpNum, startV = a.pos.get(0);
        v.initialVCFwrite(outFile);
        Integer siteEndA = a.pos.get(siteSizeA-1);
//        double[] geneticDistance = new double[sampleSize*(sampleSize-1)/2];
        int anchorStart = 0, anchorEnd = windowSize;
        BufferedWriter binfo = IOUtils.getTextWriter(outFile.replace(".vcf",".txt"));
        try {
            binfo.write("Start\tend\tmeanDist\tcontrast_number\n");
        } catch (Exception ex) {
            
        }
        for(int i = 0; i < siteSizeV; i++){
            if(v.pos.get(i) < startA||v.pos.get(i) > siteEndA) {
                v.writeVCF(i);
                continue;
            }
//            if((siteSizeA - anchorEnd) < windowSize) anchorEnd = siteSizeA - 1;
            if(anchorEnd > siteSizeA -1) anchorEnd = siteSizeA-1;
            int siteV = v.posMap.get(a.pos.get(anchorEnd));
            boolean[] pass = new boolean[sampleSize*(sampleSize-1)/2];
            int[][] anchorGenotype = getGenotype(a.genotype,anchorStart,anchorEnd);
            pass = getGeneticDistance(anchorGenotype,minComp,maxIBDDist);
           
            try{
                binfo.write(a.pos.get(anchorStart)+"\t");
                binfo.write(a.pos.get(anchorEnd)+"\t");
                binfo.write(meanDist+"\t");
                binfo.write(contrastNum+"\n");
            }catch (Exception e){
                
            }
            
            int[][] vcfGenotype = getGenotypeTrans(v.genotype,i,siteV);
//            gint[] test = getCol(vcfGenotype,1);
            while(i < siteV){
                int m = 0;
                int[] count = new int[2];
                count[0] = 0; count[1] = 0;
                ExecutorService executor = Executors.newFixedThreadPool(30);
                Count aaa = new Count(vcfGenotype,pass);
                executor.execute(aaa);
                int[] counts = aaa.getCount();
                count[0] += counts[0];
                count[1] += counts[1];
                executor.shutdown();
                if(count[1]> 2* count[0]){
                    v.writeVCF(i);
                }else if((count[0]+count[1])==0){
                    v.writeVCF(i);
                };
                i++;
                if(i%100000 ==0){
                    System.out.println(i+"\tSNPs procceeded!");
                }
            }
            anchorStart = anchorEnd;
            anchorEnd +=windowSize;
        }
        try{
            binfo.flush();
            binfo.close();
            v.close();
        } catch (Exception e){
            
        }
        
    }
    private int[] getCol(int[][] a,int index){
        int[] res = new int[a.length];
        for(int i = 0; i< a.length;i++){
            res[i] = a[i][index];
        }
        return res;
    }
    private boolean[] getGeneticDistance(int[][]anchorGenotype, int miniComp,double maxIBDDist){
        int sampleSize = anchorGenotype[0].length;
        boolean[] pass = new boolean[sampleSize*(sampleSize-1)/2];
        int m = 0;
        int mn = 0;
        contrastNum = 0;
        for(int i =0; i<sampleSize-1;i++){
            for(int j=i+1; j < sampleSize;j++){
                double dist = getTwoDistance(getCol(anchorGenotype,i),getCol(anchorGenotype,j),miniComp);
                if(Double.isNaN(dist)){
                    pass[m] = false;
                }else if(dist< maxIBDDist){
                    meanDist +=dist;
                    pass[m] = true;
                    contrastNum++;
                    mn++;
                }else{
                    meanDist +=dist;
                    pass[m] = false;
                    mn++;
                };
                m++;
            }
        }
        meanDist/=mn;
        return pass;
    }
    private double getTwoDistance(int[] a, int[]b,int miniComp){
        double dist = 0;
        double pIBS = 0;
        int sample = 0;
        for(int i = 0; i<a.length;i++){
            if(a[i]!= -9 && b[i]!=-9){
                if(a[i]*b[i]== 1){
                    pIBS += 1;
                }else if (a[i]*b[i]== -1){
                    pIBS += 0;
                }else{
                    pIBS += 0.5;
                }
                sample++;
            }
        }
        if(sample < miniComp) return Double.NaN;
        pIBS /= sample;
        dist = 1 - pIBS;
        return dist;
    }
    private int[] getCountparrallel(int[]a,int[]b){
        // res[0]: IBD count
        int[] res = new int[2];
        for(int i = 0; i < a.length;i++){
            if(a[i]!= -9 && b[i]!=-9){
                if(a[i]*b[i]== 1){
                    res[1] += 1;
                }else if (a[i]*b[i]== -1){
                    res[0] += 1;
                }else{
                    res[0] +=1;
                }
            }
        }
        return res;
    }
    private int[] getCount(int[]a,int[]b){
        // res[0]: IBD count
        int[] res = new int[2];
        for(int i = 0; i < a.length;i++){
            if(a[i]!= -9 && b[i]!=-9){
                if(a[i]*b[i]== 1){
                    res[1] += 1;
                }else if (a[i]*b[i]== -1){
                    res[0] += 1;
                }else{
                    res[0] +=1;
                }
            }
        }
        return res;
    }
    private int[][] getGenotype(List<String[]> genotype,int start,int end){
        int[][] num = new int[end-start][genotype.get(0).length];
        for(int i= start; i< end;i++){
            num[i-start] = toNumeric(genotype.get(i));
        }
        return num;
    }
    private int[][] getGenotypeTrans(List<String[]> genotype,int start,int end){
        int[][] num = new int[genotype.get(0).length][end-start];
        for(int i= start; i< end;i++){
            int[] a = toNumeric(genotype.get(i));
            for(int j = 0; j<a.length;j++){
                num[j][i-start] = a[j];
            }
        }
        return num;
    }
    public static int[] toNumeric(String[] te){
        int[] ge = new int[te.length];
        for(int j = 0; j<te.length; j++){
            if(te[j].startsWith("0/0")){
                ge[j] = -1;
            }else if(te[j].startsWith("1/1")){
                ge[j] = 1;
            }else if (te[j].startsWith(".")){
                ge[j] = -9;
            }else{
                ge[j] = 0;
            }
        }
        return ge;
    }
}
