/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

import io.IOUtils;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import vcf.vcf;

/**
 *
 * @author yaozhou
 */
public class Correlation {
    public Correlation(){
    }
    public Correlation(double[] a, double[][] b){
        this.calCor(a,b);
    }
    public void siteCor(String inFile,String outFile,String site){
        vcf genoe = new vcf();
        genoe.initialVCFread(inFile);
        genoe.readVCFByChr();
        List<String[]> geno = genoe.getGenotype();
        int a = geno.size();
        int b = geno.get(1).length;
        double[][] genotypeMatrix = new double[a][b];
        for(int i = 0; i < a ; i++){
            for (int j = 0; j < b; j++){
                if(geno.get(i)[j].startsWith("0/0")){
                    genotypeMatrix[i][j] = 0;
                }else if(geno.get(i)[j].startsWith("1/1")){
                    genotypeMatrix[i][j] = 2;
                }else if(geno.get(i)[j].startsWith("0/1")){
                    genotypeMatrix[i][j] = 1;
                }else {
                    genotypeMatrix[i][j] = Double.NaN;
                }
            }
        }
        double[] siteD = new double[b];
        Map<Integer,Integer> pos = genoe.getMap();
        int p = pos.get(Integer.parseInt(site));
        siteD = genotypeMatrix[p];
        double[] cor = calCor(siteD,genotypeMatrix);
        List < String[] > info = genoe.getInfo();
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        try {
            for (int i = 0,length = info.size(); i < length; i++ ){
                bw.write(info.get(i)[0]+"\t"+info.get(i)[1]+"\t"+cor[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Correlation.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
    public double[] calCor(double[] a,double[][] b){
        double cor[] =  new double[b.length];
        int br = b.length;
        int bc = b[0].length;
        if(a.length != bc) return cor;
        PearsonsCorrelation cc = new PearsonsCorrelation();
        for(int i = 0; i < br ; i++){
            double[] brv = getRow(b,i);
            double[][] c = omitNaN(a,brv);
//            for(int j = 0; j < c[0].length;j++ ){
//                System.out.println(c[0][j]+"\t"+c[1][j]);
//            }
//            System.out.println(".....");
            if(c[0].length < 20){
                cor[i] =1;
            }
            else{
                cor[i] = cc.correlation(getRow(c,0), getRow(c,1));
            }
            cor[i] = cor[i]*cor[i];// Rsquare
        }
        return(cor);
    }
    public double[] getRow(double[][] a,int i){
        int ar = a.length;
        int ac = a[0].length;
        double[] brv = new double[ac];
        for(int j = 0; j < ac; j++ ){
                brv[j] = a[i][j];
        }
        return brv;
    }
    public double mean(double[] a){
        double mean = Double.NaN;
        double sum =0;
        double count = 0;
        for (int i = 0; i< a.length ; i++){
            if(a[i]!=Double.NaN){
                sum += a[i];
                count++;
            }
        }
        if(count >0) mean = sum/count;
        return mean;
    }
    public boolean[] getNaNIndex(double[] a){
        boolean[] indexNaN = new boolean[a.length];
        for(int i =0; i < a.length;i++){
            if(Double.isNaN(a[i])){
                indexNaN[i] = true;
            }
            else indexNaN[i] = false;
        }
        return indexNaN;
    }
    
    public double[][] omitNaN(double[] a,double[] b){
        boolean[] aNaN = getNaNIndex(a);
        boolean[] bNaN = getNaNIndex(b);
        int length = 0;
        for(int i = 0; i < a.length;i++){
           if(!aNaN[i] & !bNaN[i]) length ++;
        }
        double[][] res = new double[2][length];
        length = 0;
        for(int i = 0; i < a.length;i++){
           if(!aNaN[i] & !bNaN[i]) {
               res[0][length] = a[i];
               res[1][length] = b[i];
               length++;
           }
        }
        return res;
    }
    
}
