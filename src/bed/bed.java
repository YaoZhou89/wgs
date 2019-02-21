/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package bed;

import htsjdk.tribble.readers.TabixReader;
import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 *
 * @author yaozhou
 */
public class bed {
    public bed (){
        
    }
    public List<Integer[]> readBed (String inFile){
        List<Integer[]>  pos = new ArrayList();
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            
            while((temp = br.readLine())!=null){
                Integer[] p = new Integer[3];
                te = temp.split("\t");
                for (int i = 0 ; i < 3; i++){
                    p[i] = Integer.parseInt(te[i]);
                }
                pos.add(p);
            }
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
        return pos;
    }
    public List<Integer[]> sortBed(List<Integer[]> pos){
        List<Integer[]> sortPos = new ArrayList();
        Integer[]  former = new Integer[3];
        former[0] = 0;
        former[1] = 0;
        former[2] = 0;
        
        for (int i = 0; i < pos.size(); i++){
            Integer[] p = pos.get(i);
            former[0] = p[0];
            former[1] = p[1];
            if(p[1] <=  former[2]){
                
                former[2] = p[1];
            }else{
                former[2] = p[2];
               
                    Integer[]  a = new Integer[3];
                    for(int j = 0; j < 3; j++){
                        a[j] = former[j];
                    }
                    sortPos.add(a);
//                }else{
//                    Integer[]  a = new Integer[3];
//                    for(int j = 0; j < 3; j++){
//                        a[j] = p[j];
//                    }
//                    sortPos.add(a);
//                }
            }
        }
        System.out.println("bed file sorted!");
        return sortPos;
    }
    public void writeBed(String outFile, List<Integer[]> sortPos){
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        try {
            for(Integer[] p : sortPos){
                for (int i = 0; i < 2; i++){
                    bw.write(p[i]+"\t");
                }
                bw.write(p[2] + "\n");
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void sort(String inFile,String outFile){
        List<Integer[]> pos =  readBed(inFile);
        List<Integer[]> sortBed =  sortBed(pos);
        writeBed(outFile,sortBed);
    }
    public void getRegion(String inFile, String outFile){
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            for (int i = 1; i < 43; i++ ){

                    BufferedReader br  = IOUtils.getTextReader(inFile+"/chr"+i+".bed");

                    String temp = "";
                    String[] te = null;
                    int sum = 0;
                    while ((temp = br.readLine())!=null){
                        te = temp.split("\t");
                        int a = Integer.parseInt(te[1]);
                        int b = Integer.parseInt(te[2]);
                        sum += (b-a);
                    }
                    bw.write("chr" + i + "\t" + sum);
                    bw.newLine();

            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void splitByChr(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextGzipReader(inFile);
            BufferedWriter bw = null;
            String temp = "";
            String chr = "";
            String[] te = null;
            boolean first = true;
            while((temp = br.readLine())!=null){
                te = temp.split("\t");
                if(te[0].equals(chr)){
                    bw.write(temp);
                    bw.newLine();
                    chr = te[0];
                }else{
                    chr = te[0];
                    if(!first){
                      bw.flush();
                      bw.close();
                    }
                    first = false;
                    System.out.println("Analyzing " + chr);
                    bw = IOUtils.getTextGzipWriter(outFile+"/depth."+chr+".txt.gz");
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void mummer2bed(String inFile,String outFile,String type,int len){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            int sum = 0;
            boolean isSNP = true;
            int start = 0;
            int pre = 0;
            boolean f = true;
            while ((temp = br.readLine())!=null){
                if(temp.equals("")) continue;
                if(!isStartWithNumber(temp)) continue;
                te = temp.split("\t");
                if(te[1].equals("\\.") | te[2].equals("\\.")) isSNP = false;
                if(type.equals("SNP")){
                    if(!isSNP) continue;
                }else{
                    if(isSNP) continue;
                }
                int pos = Integer.parseInt(te[0]);
                while((pre+1)*len < pos ){
                    if(f){
                        bw.write(te[10] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + 0);
                        bw.newLine();
                        pre++;
                    }else{
                        bw.write(te[10] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + sum);
                        bw.newLine();
                        sum = 0;
                        f = true;
                        pre++;
                    }
                }
                f = false;
                sum++;
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void repeat2bed(String inFile,String outFile,int len){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            int sum = 0;
            boolean isSNP = true;
            int start = 0;
            int pre = 0;
            boolean f = true;
            while ((temp = br.readLine())!=null){
                if(temp.equals("")) continue;
                te = temp.split("\t");                
                int pos = Integer.parseInt(te[3]);
                while((pre+1)*len < pos ){
                    if(f){
                        bw.write(te[0] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + 0);
                        bw.newLine();
                        pre++;
                    }else{
                        if(sum>len) sum = len;
                        bw.write(te[0] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + sum*1.0/len);
                        bw.newLine();
                        sum = 0;
                        f = true;
                        pre++;
                    }
                }
                f = false;
                int a = Integer.parseInt(te[4]);
                int b = Integer.parseInt(te[3]);
                if (a > (pre+1)*len ) a = (pre+1)*len;
                sum  = sum + (a-b);
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void gene2bed(String inFile,String outFile,int len){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            int sum = 0;
            int start = 0;
            int pre = 0;
            boolean f = true;
            while ((temp = br.readLine())!=null){
                if(temp.equals("")) continue;
                te = temp.split("\t");                
                int pos = Integer.parseInt(te[3]);
                while((pre+1)*len < pos ){
                    if(f){
                        bw.write(te[0] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + 0);
                        bw.newLine();
                        pre++;
                    }else{
                        if(sum>len) sum = len;
                        bw.write(te[0] + "\t" + pre*len + "\t" + (pre+1)*len + "\t" + sum);
                        bw.newLine();
                        sum = 0;
                        f = true;
                        pre++;
                    }
                }
                f = false;
                 sum++;
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(bed.class.getName()).log(Level.SEVERE, null, ex);
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
}
