/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fasta;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author yaozhou
 */
public class Fasta {
    public Fasta(){
        
    }
    public void readFasta(String inFile,String prefix,int sampleSize){
        try {
            BufferedReader br = null;
            if(inFile.endsWith(".gz")){
                br = IOUtils.getTextGzipReader(inFile);
            }else{
                br = IOUtils.getTextReader(inFile);
            }
            String temp = "";
            int i = 0;
            int t = 0;
            boolean f = true;
            BufferedWriter bw = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    i++;
                }
                if(i > sampleSize){
                    f = true;
                    i = 0;
                }
                if(f){
                    if(t>0){
                        bw.flush();
                        bw.close();
                    }
                    t++;
                    bw = IOUtils.getTextWriter(prefix+"/" + t + ".fasta");
                }
                bw.write(temp);
                bw.newLine();
                f = false;
            }
            bw.flush();
            bw.newLine();
        } catch (IOException ex) {
            Logger.getLogger(Fasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void removeDup(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile+"/1.fasta");
//            BufferedWriter bw1 = IOUtils.getTextWriter(outFile+".dup");
            String temp = "";
            String[] te = null;
            boolean write = true;
            Map<String, String> tran = new HashMap();
            Set name = new HashSet();
            int i = 0;
            int failed = 0;
            String namepre = "";
            StringBuilder seq = new StringBuilder();
            boolean first = true;
            boolean dup = false;
            String napre = "";
            boolean put = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    i++;
//                    write =true;
//                    char s  = temp.charAt(temp.length()-1);
//                        if(Character.getNumericValue(s)>1){
//                            write = false;
//                            failed ++;
//                     } 
                    String na = null;
                    if(temp.startsWith(">SME")) {
                        na = temp.split("\\.")[1] +"."+ temp.split("\\.")[2];
//                        System.out.println(na);
                    }else if(temp.startsWith(">NTA")){
                        na =temp.split("\\.")[1];
                    }else{
                        na = temp.split("\\.")[0];
                    }
                    
                    if(!first) {
                        if(dup){
                            failed++;
                            dup = false;
                            if(tran.get(napre).length() < seq.toString().length()){
                                tran.put(napre,seq.toString());
                            }
                        }else{
                            tran.put(napre,seq.toString());
                        }
                        if(name.add(na)){
                            napre=na;
                        }else{
                            dup = true;
                            napre = na;
                        };
                    }else{
                        name.add(temp.split("\\.")[0]);
                        first = false;
                    }
                    seq = new StringBuilder();
                    seq.append(temp);
                    seq.append("\n");
                }else{
                    String t = temp.replaceAll("\\.", "");
                    String tT = t.replaceAll("\\*", "");
                    seq.append(tT);
                    seq.append("\n");
                }
                
            }
            if(dup){
                failed++;
                dup = false;
                if(tran.get(napre).length() < seq.toString().length()){
                    tran.put(napre,seq.toString());
                }
            }else{
                tran.put(napre,seq.toString());
            }
            int j = 1;
            int window = tran.size()/200/23;
            int t = 0;
            for (String k : tran.values()){
                t++;
                if(t % window == 0){
                    j++;
                    bw.flush();
                    bw.close();
                    String outFiles = outFile + "/"+ j +".fasta";
                    bw = IOUtils.getTextWriter(outFiles);
                }
                bw.write(k);
            }
            bw.flush();
            bw.close();
//            bw1.flush();
//            bw1.close();
            System.out.println("Total is: " +i);
            System.out.println("Duplicated is: " + failed);
        } catch (IOException ex) {
            Logger.getLogger(Fasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public static boolean isStartWithNumber(String str) {
        Pattern pattern = Pattern.compile("[0-9]*");
        Matcher isNum = pattern.matcher(str.charAt(0)+"");
        if (!isNum.matches()) {
            return false;
        }
        return true;
    }
    public void getGenes(String inFile,String inFile2, String outFile){
        try {
            BufferedReader br1 = IOUtils.getTextReader(inFile);
            BufferedReader br2 = IOUtils.getTextReader(inFile2);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            Set name = new HashSet();
            String temp = "";
            String[] te = null;
            while((temp = br1.readLine())!=null){
                name.add(">"+temp);
            }
            boolean write = false;
            while((temp = br2.readLine())!=null){
                if(temp.startsWith(">LOC")){
                    te = temp.split("\\.");
//                    System.out.println(temp);
                    if(!name.add(te[0])){
                        write = true;  
                    }else{
                        name.remove(te[0]);
                        write = false;
                    }
                }
                if(write){
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Fasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void splitABD(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bwA = IOUtils.getTextWriter("TaA.fa");
            BufferedWriter bwB = IOUtils.getTextWriter("TaB.fa");
            BufferedWriter bwD = IOUtils.getTextWriter("TaD.fa");
            String temp = "";
            boolean A = false, B = false, D = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">TraesCS")){
                    if(temp.contains("A0")){
                        A = true;
                        B = false;
                        D = false;
                    }else if (temp.contains("B0")){
                        A = false;
                        B = true;
                        D = false;
                    }else if(temp.contains("D0")){
                        A = false;
                        B = false;
                        D = true;
                    }else{
                        System.out.println(temp);
                        A = false;
                        B = false;
                        D = true;
                    }
                }
                if(A){
                    bwA.write(temp);
                    bwA.newLine();
                }else if(B){
                    bwB.write(temp);
                    bwB.newLine();
                }else if(D){
                    bwD.write(temp);
                    bwD.newLine();
                }
            }
            bwA.flush();
            bwB.flush();
            bwD.flush();
            bwA.close();
            bwB.close();
            bwD.close();
        } catch (IOException ex) {
            Logger.getLogger(Fasta.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
