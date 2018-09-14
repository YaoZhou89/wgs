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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class Frequence {
    public Frequence(){
        
    }
    public static void sumFre(String inFile, String suffix,String outFile){
        File in = new File(inFile);
        File[] ins = IOUtils.listRecursiveFiles(in);
        File[] files = IOUtils.listFilesEndsWith(ins, suffix);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String temp = "";
        try {
            for(File subFile : files){
                
                String[] te = null;
                double[] fre = new double[2];
                double maf = 0;
                BufferedReader br = IOUtils.getTextReader(subFile.toString());
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("CHROM")){
                        System.out.println(subFile.toString());
                    } else{
                        te = temp.split("\t");
                        try{
                            maf = Double.parseDouble(te[4].split(":")[1]);
                        }catch (Exception e){
                            continue;
                        }
                        fre[1] = Double.parseDouble(te[5].split(":")[1]);
                        if(maf>fre[1]) maf = fre[1];
                        bw.write(maf+"\n");
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
           
            Logger.getLogger(Frequence.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
