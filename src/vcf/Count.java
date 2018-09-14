/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import java.util.concurrent.Callable;

/**
 *
 * @author yaozhou
 */
public class Count implements Runnable{
    private int[][] genotype;
    private int[] res = new int[2];
    private boolean[] pass;
//    private int s1;
    public Count(int[][] genotype,boolean[] pass){
        this.genotype = genotype;
        this.pass = pass;
    }
    public void run(){
        res[0] = 0;
        res[1] = 0;
        int[] r = new int[2];
        int m = 0;
        for(int j = 0;j<genotype.length-1;j++){
            for (int i = j +1 ; i< genotype.length;i++){
                
                if(!pass[m]) {
                    m++;
                    continue;
                }
                int[] a = genotype[j];
                int[] b = genotype[i];
                r = getCount(a,b);
                res[0] += r[0];
                res[1] += r[1];
                m++;
            }  
        }
        
        
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
   
    public int[] getCount() {
        return res;
    }      
}
