/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import htsjdk.tribble.readers.TabixReader;
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
public class CompareVCF {
    public CompareVCF(String in1, String in2, String suffix){
        this.compareVCF(in1,in2,suffix);
    }
    public CompareVCF(){
    }
    
    public void compareVCF(String in1, String in2, String suffix){
        File inFILE = new File(in1);
        File[] a = IOUtils.listRecursiveFiles(inFILE);
        File[] subFile = IOUtils.listFilesEndsWith(a, suffix);
        Set RefPos = new HashSet();
        String temp = "",chr="";
        BufferedReader bref = IOUtils.getTextReader(in2);
        String outFile = in1+"/CScount"+suffix+".txt";
       
        String[] te = null;
        int count = 0;
        int refc =0;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
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
                TabixReader br = new TabixReader(f.toString());
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
    public void compare(String in1, String in2){
        try {
            TabixReader br1 = new TabixReader(in1);
            TabixReader br2 = new TabixReader(in2);
            Set S1 = new HashSet();
            String temp;
            String[] te = null;
            while((temp = br1.readLine())!=null){
                if(temp.startsWith("#")) continue;
                te = temp.split("\t");
                S1.add(te[0]+"_"+te[1]+"_"+te[3]+"_"+te[4]);
            }
            System.out.println(in1 +"\treaded!");
            BufferedWriter bwd = IOUtils.getTextWriter(in2.replace(".vcf.gz", ".unique.vcf"));
            BufferedWriter bws = IOUtils.getTextWriter(in2.replace(".vcf.gz", ".intersect.vcf"));
            while((temp = br2.readLine())!=null){
                if(temp.startsWith("#")){
                    bwd.write(temp);
                    bwd.newLine();
                    bws.write(temp);
                    bws.newLine();
                }else{
                    te = temp.split("\t");
                    if(S1.add(te[0]+"_"+te[1]+"_"+te[3]+"_"+te[4])){
                        bwd.write(temp);
                        bwd.newLine();
                    }else{
                        bws.write(temp);
                        bws.newLine();
                    }
                }
            }
            bwd.flush();
            bwd.close();
            bws.flush();
            bws.close();
        } catch (IOException ex) {
            Logger.getLogger(CompareVCF.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
