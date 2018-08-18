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
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author yaozhou
 */
public class CompareVCF {
    public CompareVCF(String in1, String in2, String suffix){
        this.compareVCF(in1,in2,suffix);
    }
    public void compareVCF(String in1, String in2, String suffix){
        File inFILE = new File(in1);
        File[] a = IOUtils.listRecursiveFiles(inFILE);
        File[] subFile = IOUtils.listFilesEndsWith(a, suffix);
        Set RefPos = new HashSet();
        String temp = "",chr="";
        BufferedReader bref = IOUtils.getTextReader(in2);
        String outFile = in1+"/CScount"+suffix+".txt";
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String[] te = null;
        int count = 0;
        int refc =0;
        try{
            while((temp = bref.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    RefPos.add(te[0]+"_"+te[1]);
                    refc++;
                }
            }
            System.out.println("Reference read!");
            bref.close();
            
            for (File f : subFile){
                int c = 0;
                System.out.println("now processing\t"+f.toString());
                BufferedReader br = IOUtils.getTextReader(f.toString());
                while((temp = br.readLine())!=null){
                    if(!temp.startsWith("#")){
                        te = temp.split("\t");
                        if(!RefPos.add(te[0]+"_"+te[1])){
                            chr = te[0];
                            c++;
                            count++;
                        }
                    }
                }
                bw.write(chr+"\t"+ c);
                bw.newLine();
                br.close();
                System.out.println("Number of CS SNPs found in this file:\t"+c);
            }
            bw.write("All\t"+count);
            bw.newLine();
            bw.write("Ref allele is: "+ refc);
            bw.newLine();
            bw.flush();
            bw.close();
            System.out.println("Number of CS SNPs counted:\t"+count);
        }catch(Exception e){
            
        }
    }
}
