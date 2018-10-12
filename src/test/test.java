/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author yaozhou
 */
public class test {
    public test(String inFile,String type){
        if(type.equals("null")){
            this.getTest(inFile);
        }
        if(type.equals("new")){
            this.getGenome(inFile);
        }
        if(type.equals("test")){
            String S = "NNNATGCNNNCTANN";
            String[] s = S.split("N");
            for(int i = 0; i< s.length;i++){
                System.out.println(s[i]+"\t"+s[i].length());
            }
        }
        
    }
    private void getTest(String inFile){
        try {
            BufferedReader vcf;
            if(inFile.endsWith("gz"))  vcf = IOUtils.getTextGzipReader(inFile);
            else  vcf = IOUtils.getTextReader(inFile);
            List<String> chrInfo = new ArrayList();
//            BufferedReader bc = IOUtils.getTextReader("/public-supool/home/yaozhou/data/ref/wheat/chr.txt");
            String temp=null;
//            while ((temp = bc.readLine())!=null){
//                chrInfo.add(temp);
//            }
//            bc.close();
//            Set<String> chr = new HashSet<String>();
            String[] tem =null;
//            int i = -1;
//            String chr1 = "aa";
//            BufferedWriter chra = IOUtils.getTextWriter("/public-supool/home/yaozhou/data/ref/wheat/BM.bed");
            StringBuilder sb = new StringBuilder();
           
            while ((temp = vcf.readLine())!=null){
//                tem = temp.split("\t");
                
                if(temp.startsWith(">")){
                    System.out.println(temp);
                    System.out.println(sb.toString().length());
                    sb = new StringBuilder();
                }else{
                    sb.append(temp);
                }
                
//                if(!chr1.equals(tem[0])){
//                    chr1 = tem[0];
//                    i++;
//                }
//                temp = temp.replaceFirst(tem[0],chrInfo.get(i));
//                chra.write(temp);
//                chra.newLine();
//    //                if(!temp.startsWith("#")){
////                    tem = temp.split("\t");
////                    chr.add(tem[0]);
////                }
//                i++;
            }
//            System.out.println(Integer.toString(i));
//            chra.flush();
//            chra.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    private void getGenome(String inFile){
        try {
            BufferedReader br;
            if(inFile.endsWith(".gz"))  br= IOUtils.getTextGzipReader(inFile);
            else br =IOUtils.getTextReader(inFile);
            String outFile = inFile.replace(".fa",".fasta");
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = null;
            boolean w = false;
            while((temp=br.readLine())!=null){
               if(temp.startsWith(">")){
                   if(temp.split("\t")[0].indexOf("Lachesis_group")>-1){
//                       bw.write(temp);
//                       bw.newLine();
                       w = true;
                   }else{
                       w= false;
                   }
               } 
               if(w){
                       bw.write(temp);
                       bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(test.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
