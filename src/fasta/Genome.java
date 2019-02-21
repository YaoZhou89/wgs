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
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author yaozhou
 */
public class Genome {
    int chrNum = 0;
    Map chrInfo = new HashMap();
    Map<String, String> seq = new HashMap();
    public Genome(){
        
    }
    private void read(String inFile){
        
        try {
            BufferedReader br = null;
            if(inFile.endsWith(".gz")){
                br = IOUtils.getTextGzipReader(inFile);
            }else {
                br = IOUtils.getTextReader(inFile);
            }
            String temp = "";
            String chrName = "";
            StringBuilder fasta  = new StringBuilder();
            boolean first = true;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    System.out.println(temp);
                    if(!first){
                       seq.put(chrName,fasta.toString());
                    }
                    chrName = temp;
                    fasta = new StringBuilder();
                    first = false;
                }else{
                    fasta.append(temp);
                    fasta.append("\n");
                }
            }
            seq.put(chrName,fasta.toString());
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private void write(String outFile){
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            for(int i = 1; i< 13;i ++){
                String chr = ">Chr"+i;
                System.out.println(chr);
                bw.write(chr);
                bw.newLine();
                bw.write(seq.get(chr));
            }
            bw.write(">ChrUn");
            bw.newLine();
            bw.write(seq.get(">ChrUn"));
            bw.write(">Chl");
            bw.newLine();
            bw.write(seq.get(">Chl"));
            bw.flush();
            bw.close();
        } catch (IOException ex) {
                Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
            }
    }
    public void getSort(String inFile, String outFile){
        read(inFile);
        write(outFile);
    }
    public void readByChromosome(String inFile,int a, int b){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            StringBuilder chr = new StringBuilder();
            boolean test = false, chr1 =true;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">1")){
                    test = true;
                    System.out.println(temp);
                    chr = new StringBuilder();
                }else if (test){
                    chr.append(temp);
                }
                if(temp.startsWith(">2")){
                   break;
                }
            }
            System.out.println(chr.toString().substring(a, b));
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getChrInfo(String inFile,String outFile){
        BufferedReader br = IOUtils.getTextGzipReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String temp = "";
        int chrLength = 0;
        try {
            while((temp = br.readLine())!=null){ 
                if(temp.startsWith(">")){
                    bw.write(temp + "\t");
                    if(chrLength > 0){
                        bw.write(chrLength);
                        chrLength  = 0;
                        bw.newLine();
                    }
                }else{
                    chrLength += temp.length();
                }
            }
            bw.write(chrLength);
            bw.newLine();
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getBWAformat(String inFile, String outFile,double pattern){
        BufferedReader br = IOUtils.getTextGzipReader(inFile);
        BufferedWriter bw = IOUtils.getTextGzipWriter(outFile);
//        String[] p = inFile.split("/");
//        String path = inFile.split(p[p.length-1])[0];
//        BufferedWriter binof = IOUtils.getTextWriter(path+"/Readme.txt");
        String temp = "";
        String[] te = null;
        StringBuilder chr = null;
        boolean write = false;
        int chrNum = 0;
        boolean unStart = false;
        StringBuilder interval = null;
        Integer length = 0;
        int lineLength = 0;
        boolean first= true;
        try {
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    if(chrNum < pattern || !unStart){
                        chrNum++;
                        bw.write(">" + chrNum);
                        System.out.println(">" + chrNum);
                        bw.newLine();
//                        chr = new StringBuilder();
//                        interval = new StringBuilder();
                        length = 0;
                        unStart = true;
                    }else{
                        continue;
                    }
                    
//                        first = true;
                }else{
                    if(first) {
                        lineLength = temp.length();
                        first = false;
                    }
                    if(length < 400000001 ){
                        length += lineLength;
                        bw.write(temp);
                        if(temp.length()!=lineLength){
                            String Ns = StringUtils.repeat("N", lineLength-temp.length());
                            bw.write(Ns);
                        }
                        bw.newLine();
                    }else{
                        unStart = false;
                        if(temp.endsWith("N")||length > 500000000){
                            bw.write(temp);
                            bw.newLine();
                            chrNum++;
                            bw.write(">"+chrNum);
                            bw.newLine();
                            System.out.println(">" + chrNum);
                            length = 0;
                        }else{
                            length += lineLength;
                            bw.write(temp);
                            bw.newLine();
                        }
                    }
                    
                    
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void initialGenome(String inFile){
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inFile);
            String temp = "";
            int chrSize = 0;
            String chrName = "";
            StringBuilder fasta  = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    chrName = temp.split(">")[1];
                    if(chrSize > 0){
                        int s = 0;
                        s = chrSize;
                        chrInfo.put(chrName, s);
                        chrSize = 0;
                        seq.put(chrName,fasta.toString());
                    }
                    fasta = new StringBuilder();
                }else{
                    fasta.append(temp);
                }
                
            }
            int s = 0;
            s = chrSize;
            chrInfo.put(chrName, s);
            chrSize = 0;
            seq.put(chrName,fasta.toString());
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
       
    }
    public String getSeq(String chr){
        return seq.get(chr);
    }
}
