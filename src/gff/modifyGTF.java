/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gff;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class modifyGTF {
    public modifyGTF(String inFile, String outFile){
        this.getGTF(inFile,outFile);
    }
    public modifyGTF(String inFile, String outFile,String coor){
        this.getGTFcor(inFile,outFile,coor);
    }
    private void getGTF(String inFile, String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outFile);
            Set type =  new HashSet();
            type.add("exon");
            type.add("CDS");
            type.add("stop_codon");
            type.add("start_codon");
            String temp ="";
            while((temp = br.readLine())!=null){
                String[] tem = temp.split("\t");
                if(type.add(tem[2])){
                    type.remove(tem[2]);
                    bw.write(temp);
                    bw.newLine();
                }else{
                    StringBuilder s = new StringBuilder(temp);
                    s.append("gene_biotype \"protein_coding\";");
                    bw.write(s.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private void getGTFcor(String inFile, String outFile, String coor){
        try {
            BufferedWriter bw;
            BufferedReader br, bc;
            if(inFile.endsWith(".gz")) br = IOUtils.getTextGzipReader(inFile);
            else br = IOUtils.getTextReader(inFile);
            if(outFile.endsWith(".gz")) bw = IOUtils.getTextGzipWriter(outFile);
            else bw = IOUtils.getTextWriter(outFile);
            bc = IOUtils.getTextReader(coor);
            String temp ="";
            String[][] cor = new String[45][6];
//            HashMap <String,ArrayList> co = new HashMap();
            int line = 0;
            while((temp = bc.readLine())!=null){
                String[] tem = temp.split("\t");
//                ArrayList<String> chr = new ArrayList<String>();
                for(int i =0 ; i< tem.length;i++){
                    cor[line][i] = tem[i];
//                    chr.add(temp);
                }
                line++;
            }
            while((temp = br.readLine())!=null){
                String[] tem = temp.split("\t");
                StringBuilder ntem = new StringBuilder();
                for(int i = 0; i < cor.length;i++){
                    if(tem[0].equals(cor[i][3])){
                        if(  Integer.parseInt(tem[4]) < Integer.parseInt(cor[i][5])){
                            tem[0] = cor[i][0];
                        }else{
                            tem[0] = cor[i+1][0];
                            tem[3] = Integer.toString(Integer.parseInt(tem[3])- Integer.parseInt(cor[i][5]));
                            tem[4] = Integer.toString(Integer.parseInt(tem[4])- Integer.parseInt(cor[i][5]));
                        }
                        ntem = getString(tem);
                        bw.write(ntem.toString());
                        bw.newLine();
                        break;
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    private StringBuilder getString(String[] s){
        StringBuilder ss = new StringBuilder();
        for(int i=0; i< s.length-1; i++){
            ss.append(s[i]);
            ss.append("\t");
        }
        ss.append(s[s.length-1]);
        return ss;
    }
}
