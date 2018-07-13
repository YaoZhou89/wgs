/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fastq;

// import format.table.RowTable;
import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
//import utils.IOUtils;
//import xuebo.analysis.annotation.FStringUtils;

/**
 *change the format of GBS data to TASSEL format, 
 * and also report the result of each sample
 * @author yaozhou
 */
public class SplitSampleGBS{
    public SplitSampleGBS(String inFolder,String sampleInfor){
        this.getSplited(inFolder,sampleInfor);
    }
    private void getSplited(String infileDirS,String sampleInfor){
//        RowTable rt=new RowTable(sampleInfor);
        HashMap hm1=new HashMap();
        HashMap hm2=new HashMap();
        List sample=new ArrayList();
        BufferedReader binfo = IOUtils.getTextReader(sampleInfor);
        String temp = null;
        String[] tem = null;
        Set b1 = new HashSet();
        Set b2 = new HashSet();
        try {
            temp = binfo.readLine();
            while ((temp=binfo.readLine())!=null){
                tem = temp.split("\t");
                StringBuilder barcode1 = new StringBuilder();
                StringBuilder barcode2 = new StringBuilder();
                barcode1.append(tem[3]);
                barcode1.append(tem[4]);
                barcode2.append(tem[4]);
                barcode2.append(tem[3]);
                b1.add(barcode1);
                b2.add(barcode2);
                hm1.put(barcode1.toString(),tem[2]);
                hm2.put(barcode2.toString(),tem[2]);
                sample.add(tem[2]);
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].getName().endsWith(".gz")){
                nameSet.add(fs[i].getName().split("_")[0]);
            }
        }
        int[] cnt=new int[sample.size()];
//        File outFileDir = new File(infileDirS+"/splited");
//        int dir = 0;
//        nameSet.parallelStream().forEach((String name) -> {
        for(String name : nameSet){
            String infile1 = new File (infileDirS, name+"_R1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_R2.fq.gz").getAbsolutePath();
            String outfile = new File (infileDirS, name+"_count.txt").getAbsolutePath();
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;     
            int dir = 0 ;
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);   
//                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                while ((temp1 = br1.readLine()) != null) {
//                    System.out.println(temp1);
                    String fc = temp1.split(":")[0].replace("@", "");
                    String lane = temp1.split(":")[1];
                    br2.readLine();
                    seq1=br1.readLine();
                    seq2 = br2.readLine();
                    int l1 = 0;
                    int l2 = 0;
                    if(seq1.indexOf("GATCC")<1|seq1.indexOf("GATCC")>10){
                        l1 = seq1.indexOf("CGG")+1;
                        if(l1<1 | l1>10) continue;
                        l2 = seq2.indexOf("GATCC")+1;
                        if(l2 <1 | l2>10) continue;
                    }else{
                       l1 = seq1.indexOf("GATCC");
                       if(l1<1 | l1>10) continue;
                       l2 = seq2.indexOf("CGG");
                       if(l2 <1 | l2>10) continue;
                    };
                    StringBuilder seqbar = new StringBuilder(seq1.substring(0,l1));
                    seqbar.append(seq2.substring(0,l2));
                    if(!b1.add(seqbar)|b2.add(seqbar)){
                        br1.close();
                        br2.close();
                        getCount(infile1,infile2,hm1,outfile,cnt,sample,fc,lane);
                        break;
                    }
                    if(b1.add(seqbar)|!b2.add(seqbar)){
                        br1.close();
                        br2.close();
                        getCount(infile1,infile2,hm2,outfile,cnt,sample,fc,lane);
                        break;
                    }
                }
                
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

        };
    }
    public String[][] RowTable(String sampleInfo){
        BufferedReader br = IOUtils.getTextReader(sampleInfo);
        String temp = null;
        String[][] table = new String[1][1];
        return table;
    }
    private void getCount(String infile1, String infile2,HashMap hm,String outfile,
            int[] cnt, List sample, String fc, String lane){
        try {
            BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
            BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            bw.write("Taxa\tCount\n");
            String name = infile1.split("/")[infile1.split("/").length-1];
            String out1 = infile1.replace(name,fc+"_"+lane+"_fastq.txt.gz");
            BufferedWriter bw1 = IOUtils.getTextGzipWriter(out1);
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;
            int num = 0;
            int j = 0;
            while((temp1 = br1.readLine())!=null){
//                System.out.println(temp1);
                temp2 = br2.readLine();
                seq1 = br1.readLine();
                seq2 = br2.readLine();
                int l1 = 0;
                int l2 = 0;
                if(seq1.indexOf("GATCC")<1|seq1.indexOf("GATCC")>10){
                    l1 = seq1.substring(4, 12).indexOf("CGG");
                    if(l1>-1){
                        l1 = l1+4;
                    }else{
                        continue;
                    }
                    l2 = seq2.substring(4, 12).indexOf("GATCC");
                    if(l2>-1) {
                        l2 = l2+4;
                    }else{
                        continue;
                    }
                    
                }else{
                    l1 = seq1.substring(4, 12).indexOf("GATCC");
                    if(l1>-1){
                        l1 = l1+4;
                    }else{
                        continue;
                    }
                    l2 = seq2.substring(4, 12).indexOf("CGG");
                    if(l2>-1) {
                        l2 = l2+4;
                    }else{
                        continue;
                    }
                };
//                if(seq1.startsWith("CAGCT")){
//                    if(seq2.startsWith("GCCGG")){
//                        j++;
//                        System.out.println(j);
//                    }
//                }
//                if(seq2.startsWith("CAGCTGATCC")){
//                    System.out.println("seq1");
//                }
                StringBuilder seqb = new StringBuilder(seq1.substring(0,l1));
                seqb.append(seq2.substring(0,l2));
                if(hm.get(seqb.toString())!=null){
                    int index=sample.indexOf(hm.get(seqb.toString()));
                    cnt[index]++;
                    bw1.write(temp1);bw1.newLine();
                    seq1.replace(seq1.substring(0,l1), seqb.toString());
                    seq2.replace(seq2.substring(0,l2), seqb.toString());
                    bw1.write(seq1);bw1.newLine();
                    bw1.write(br1.readLine());bw1.newLine();
                    bw1.write(br1.readLine());bw1.newLine();
                    bw1.write(temp2); bw1.newLine();
                    bw1.write(seq2);bw1.newLine();
                    bw1.write(br2.readLine());bw1.newLine();
                    bw1.write(br2.readLine());bw1.newLine();
                }else{
                    num++;
                    br1.readLine();br1.readLine();
                    br2.readLine();br2.readLine();
                }
            }
            for(int i=0;i<cnt.length;i++){
                bw.write(sample.get(i)+"\t"+cnt[i]);
                bw.newLine();
            }
            br1.close();br2.close();
            bw.write("Unknown\t"+num);
            bw.newLine();
            bw.flush();bw.close();
            System.out.println("unknown num is: " + num);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}

