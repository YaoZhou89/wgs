/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;

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
            cor[i] = cor[i]*cor[i];
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
