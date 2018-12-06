/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fasta;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class SubstractFromFasta {
    public SubstractFromFasta(String inFile,String names,String outFile){
        this.subNames(inFile, names,outFile);
    }
    public SubstractFromFasta(String inFile,String outFile){
        this.getSplited(inFile, outFile);
    }
    private void getSplited(String inFile,String outFile){
        try {
            File Genes = new File(outFile);
            if(!Genes.isDirectory()){
                Genes.mkdirs();
            }
            String temp = "";
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = null;
            int i = 0;
            boolean Write = false;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    if(bw!=null){
                        bw.flush();
                        bw.close();
                    }
                    bw = IOUtils.getTextWriter(outFile+"/"+i+".fa");
                    i++;
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
            Logger.getLogger(SubstractFromFasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private void subNames(String inFile, String names,String outFile) {
        BufferedReader br,bn ;
        BufferedWriter bw,bw1;
        if(inFile.endsWith("gz")){
            br = IOUtils.getTextGzipReader(inFile);
            bw = IOUtils.getTextGzipWriter(outFile);
            bw1 = IOUtils.getTextGzipWriter(outFile+".index");
        }
        else {
            br = IOUtils.getTextReader(inFile);
            bw = IOUtils.getTextWriter(outFile);
        }
        bn = IOUtils.getTextReader(names);
        Set set  =   new  HashSet(); 
        String temp = null,temp1 = null,temp3 = null;
        int i = 0;
       
        try {
            while ((temp = bn.readLine())!=null){
                set.add(temp);
            }
            
            while ((temp = br.readLine())!=null){
                i++;
                if(i%2 == 1){
                    temp1 = temp;
                }else{
                    if(!set.add(temp1)){
//                    System.out.println(temp.substring(temp.length()-104,temp.length()-103));
                        bw.write(temp1+"\n");
                        bw.write(temp+"\n");
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in Reading: " + inFile);
        }
    }
}
