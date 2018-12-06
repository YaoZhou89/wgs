/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class singleRegion {
    public singleRegion(){
        
    }
    public void getSR(String inFile,String outFile, double Mdep,double rate){
        BufferedReader br = IOUtils.getTextGzipReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String[] depth = null;
        String temp = "";
        List<Integer[]> de = new ArrayList(1100);
        int posPre = 0;
        int posPos = 1;
        String chr = "";
        boolean start = true;
        boolean end = false;
        List<String> map = new ArrayList(1100);
        String pos = "";
        int SNP = 0;
        try {
            while((temp = br.readLine())!=null){
                SNP++;
                if(SNP%100000000==0){
                    System.out.println("Processing: " + SNP/1000000 +"M...");
                }
                depth = temp.split("\t");
                int length = depth.length;
                Integer[] dep = new Integer[length-2];
                if(start){
                   chr = depth[0];
                   start = false;
                   posPre = Integer.parseInt(depth[1]);
                }
                posPos = Integer.parseInt(depth[1]);
                pos = "rs"+depth[0]+"_"+depth[1];
                if(!depth[0].equals(chr)){
                    System.out.print("Chromosome:\t" + depth[0]);
//                    for (String[] a : map){
//                        bw.write(a[0]+"\t"+a[1]+"\t"+a[2]);
//                        bw.newLine();
//                    }
                    de.clear();
                    map.clear();
                    posPre = Integer.parseInt(depth[1]);
                    end = true;
                }
                chr = depth[0];
                map.add(pos);
                for (int i = 0; i < length - 2; i++){
                   dep[i] = Integer.parseInt(depth[i+2]);
                }
                de.add(dep);
                if(!end && (posPos - posPre) > 1000){
                    boolean p = checkSR(de,Mdep,rate);
                    if(p){
                        for (String a : map){
                            bw.write(a);
                            bw.newLine();
                        }
                    }
                    de.clear();
                    map.clear();
                    start = true;
                }else if(end){
                    end = false;
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private boolean checkSR(List<Integer[]> de, double Mdep,double rate){
        boolean pass = false;
        int p = 0;
        for(Integer[] a : de){
            int sum = 0;
            for (int i = 0, length = a.length; i < length;i++){
                sum += a[i];
            }
            if(sum > Mdep*0.9 && sum < Mdep * 1.1) {
                p++;
            }
        }
        if(p > de.size()*rate) pass = true;
        return pass;
    }
    public void getSS(String inFile,String outFile, double a,double b){
        BufferedReader br = IOUtils.getTextGzipReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String[] depth = null;
        String temp = "";
        String pos = "";
        long  SNP = 0;
        int sum = 0;
        int length = 0;
        try {
            while((temp = br.readLine())!=null){
                SNP++;
                if(SNP%100000000==0){
                    System.out.println("Processing: " + SNP/1000000 +"M...");
                }
                depth = temp.split("\t");
                
                pos = "rs"+depth[0]+"_"+depth[1];
                sum = 0;
                length = depth.length - 2;
                for(int i = 0 ; i < length ; i++){
                    sum += Integer.parseInt(depth[i+2]);
                }
                if(sum > a && sum < b){
                    bw.write(pos);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private double getMean(String[] de){
        int sum = 0;
        int length = de.length - 2;
        for(int i = 0 ; i < length ; i++){
            sum += Integer.parseInt(de[i+2]);
        }
        return (double) sum;
    }
    public void getBed(String inFile, String outFile){
        BufferedReader br = null;
        if(inFile.endsWith("gz")){
            br = IOUtils.getTextGzipReader(inFile);
        }else {
            br = IOUtils.getTextReader(inFile);
        }
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String temp = "";
        String[] te = null;
//        String t  = "rs89";
//        String[] tt = t.split("rs");
//        System.out.println(tt[1]);
        String chr = "";
        int begin = 0;
        int end = 0;
        boolean start = true;
        System.out.println("Generating bed file..");
        int pre = 0;
        try {
            while ((temp = br.readLine())!=null){
                te = temp.split("_");
                if(te.length<2) continue;
                chr = te[0].split("rs")[1];
                if(start) {
                    begin = Integer.parseInt(te[1]) - 1;
                    pre = begin;
                    start = false;
                }
                end = Integer.parseInt(te[1]);
                if((end - pre) > 20){
                    int len = (pre -begin );
                    if( len > 150){
                        bw.write(chr + "\t" + begin + "\t" + pre + "\t" + len + "\n");
                    }
                    start = true;
                }
                pre = end;
            }
            bw.write(chr + "\t" + begin + "\t" + Integer.parseInt(te[1]) + "\t" + (end - begin)+ "\n");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void toChr(String inFile,String outFile){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedWriter bw = null;
        String temp = "";
        String[] te  = null;
        String chr = "";
        boolean start = true;
        boolean first = true;
        String pre = "";
        try {
            while((temp = br.readLine())!=null){
                te = temp.split("_")[0].split("rs");
                if(start){
                    chr = te[1];
                    if(!first){
                        bw.flush();
                        bw.close();
                        first = false;
                    }
                    bw = IOUtils.getTextWriter(outFile+"."+chr+".txt");
                    if(!pre.equals("")){
                        bw.write(pre);
                        bw.newLine();
                    }
                    start = false;
                }
                if(chr.equals(te[1])){
                    bw.write(temp);
                    bw.newLine(); 
                }else{
                    pre = temp;
                    start = true;
                }
                
            }
            bw.flush();
            bw.close();
                    
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getLen(String inFile, String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            while ((temp = br.readLine())!=null){
                te = temp.split("\t");
                int len = Integer.parseInt(te[2]) - Integer.parseInt(te[1]);
                bw.write(te[0]+"\t"+te[1]+"\t"+te[2]+"\t"+len+"\n");
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getBedTotal(String inFile, String outFile, String region){
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            int start = Integer.parseInt(region.split(" ")[0])-1;
            int end = Integer.parseInt(region.split(" ")[1]);
            String[] te = null;
            String temp = "";
            while ((temp = br.readLine())!=null){
                te = temp.split("\t");
                if (Integer.parseInt(te[1])> start){
                    if(Integer.parseInt(te[1]) > end) break;
                    int sum = 0;
                    for(int i = 2; i < te.length; i ++){
                        sum += Integer.parseInt(te[i]);
                    }
                    bw.write(te[0]+"\t"+te[1]+"\t"+sum);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(singleRegion.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
