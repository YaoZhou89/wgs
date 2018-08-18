/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.logging.Level;
import java.util.logging.Logger;
import math.Correlation;

/** Filt SNPs based on the local LD.
 * if one SNP is obviously not in the same LD with the adjacent SNPs, it will be removed;
 *@param windowSize: how many adjacent SNPs will be considered
 * @author yaozhou
 */
public class LDFilter {
    public LDFilter(){
        
    }
    public LDFilter(String inFile, String outFile,int windowSize,double threshold){
        this.calLD(inFile, outFile, windowSize, threshold);
    }
    
    public static Integer[] calLD(String inFile, String outFile,int windowSize,double threshold1){
        Integer[] num = new Integer[2];
        try {
            StringBuilder header = new StringBuilder();
            String temp;
            String[] tem;
            Double threshold = threshold1 + 1/windowSize;
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outFile+".cor.txt");
//            StringBuilder sites = new StringBuilder();
            int currentSNP = 0, matrixSNP = 0;
            double[][] genotypeMatrix = new double[windowSize][];
            boolean redefine = true;
            Map snp = new HashMap(windowSize);
            int ldrm  = 0;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    header.append(temp+"\n");
                }else{
                    matrixSNP = currentSNP % windowSize;
                    currentSNP ++;
                    if(currentSNP%100000 == 0){
                        System.out.println("Analyzing " + currentSNP+" SNPs....");
                    }
                    tem = temp.split("\t");
                    if(redefine){
                        bw.write(header.toString());
                        genotypeMatrix = new double[windowSize][tem.length-9];
                        redefine = false;
                    }
                    
                    snp.put(matrixSNP,temp);
                  
                    for(int site = 0; site < tem.length -9;site++){
                        if(tem[site+9].startsWith("0/0")){
                            genotypeMatrix[matrixSNP][site] = 0;
                        }else if(tem[site+9].startsWith("1/1")){
                            genotypeMatrix[matrixSNP][site] = 2;
                        }else if(tem[site+9].startsWith("0/1")){
                            genotypeMatrix[matrixSNP][site] = 1;
                        }else {
                            genotypeMatrix[matrixSNP][site] = Double.NaN;
                        }
//                        bw1.write(genotypeMatrix[matrixSNP][site]+"\t");
                    }
                    if(currentSNP > windowSize - 1){
                        Correlation cor = new Correlation();
                        int index = (matrixSNP+1) % windowSize;
                        double[] cur = cor.getRow(genotypeMatrix, index);
                        double[] cr = cor.calCor(cur, genotypeMatrix);
                        double m = cor.mean(cr);
                        if(m < threshold){
                            ldrm++;
                        }else{
                            bw.write(snp.get(index).toString());
                            bw.newLine();
                        }
                    }
                }
            }
            // hand left 
            Correlation cor = new Correlation();
            
            for(int i = 0; i < currentSNP % windowSize; i++){
                double[] cur = cor.getRow(genotypeMatrix, i);
                double[] cr = cor.calCor(cur, genotypeMatrix);
                double m = cor.mean(cr);
                if(m < threshold){
                    ldrm++;
                }else{
                    bw.write(snp.get(i).toString());
                    bw.newLine();
                }
            }
            num[0] = currentSNP;
            num[1] = ldrm;
            System.out.println("Total SNPs:\t" +currentSNP+"\nSNPs removed:\t"+ldrm);
            bw.flush();
            bw.close();
//            bw1.flush();
//            bw1.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return num;
    }
    // windowSize: kb
    public static void calLD_physical(String inFile, String outFile,int windowSize,double threshold){
        try {
            StringBuilder header = new StringBuilder();
            String temp;
            String[] tem;
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outFile+".cor.txt");
//            StringBuilder sites = new StringBuilder();
            int currentSNP = 0, matrixSNP = 0;
            double[][] genotypeMatrix = new double[windowSize][];
            boolean redefine = true,change = false;
            Map snp = new HashMap();
            int ldrm  = 0;
            int pos1 = 0;
            int pos2 = 0;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    header.append(temp+"\n");
                }else{
                    
                    matrixSNP = currentSNP % windowSize;
                    currentSNP ++;
                    if(currentSNP%10000 == 0){
                        System.out.println("Analyzing " + currentSNP+" SNPs....");
                    }
                    
                    tem = temp.split("\t");
                    if(redefine){
                        bw.write(header.toString());
                        genotypeMatrix = new double[windowSize][tem.length-9];
                        redefine = false;
                    }
                    
                    if(pos1 == pos2){
                        pos1 = Integer.parseInt(tem[0]);  
                    }else{
                        pos2 = Integer.parseInt(tem[0]);
                        if((pos2 - pos1 )>windowSize) change = true;
                        if(change) pos1 = pos2;
                    }
                    
                    snp.put(matrixSNP,temp);
                  
                    for(int site = 0; site < tem.length -9;site++){
                        if(tem[site+9].startsWith("0/0")){
                            genotypeMatrix[matrixSNP][site] = 0;
                        }else if(tem[site+9].startsWith("1/1")){
                            genotypeMatrix[matrixSNP][site] = 2;
                        }else if(tem[site+9].startsWith("0/1")){
                            genotypeMatrix[matrixSNP][site] = 1;
                        }else {
                            genotypeMatrix[matrixSNP][site] = Double.NaN;
                        }
//                        bw1.write(genotypeMatrix[matrixSNP][site]+"\t");
                    }
                    if(currentSNP > windowSize - 1){
                        Correlation cor = new Correlation();
                        int index = (matrixSNP+1) % windowSize;
                        double[] cur = cor.getRow(genotypeMatrix, index);
                        double[] cr = cor.calCor(cur, genotypeMatrix);
                        double m = cor.mean(cr);
                        if(m < threshold){
                            ldrm++;
                        }else{
                            bw.write(snp.get(index).toString());
                            bw.newLine();
                        }
                    }
                }
            }
            // hand left 
            Correlation cor = new Correlation();
            for(int i = 0; i < windowSize; i++){
                double[] cur = cor.getRow(genotypeMatrix, i);
                double[] cr = cor.calCor(cur, genotypeMatrix);
                double m = cor.mean(cr);
                if(m < threshold){
                    ldrm++;
                }
                else{
                    bw.write(snp.get(matrixSNP).toString());
                    bw.newLine();
                } 
            }
            System.out.println("Total SNPs:\t" +currentSNP+"\nSNPs removed:\t"+ldrm);
            bw.flush();
            bw.close();
//            bw1.flush();
//            bw1.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
}
