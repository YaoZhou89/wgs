/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import htsjdk.tribble.readers.TabixReader;
import io.IOUtils;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author yaozhou
 */
public class vcf {
    int snpNum = 0;
    Set chr = new HashSet();
    String chrom ="";
//    List<String> ID = new ArrayList();
    int sampleSize;
    String[] sample;
//    Map<String , String[]> genotype = new HashMap();
    List<String[]> genotype = new ArrayList();
    Map<Integer,Integer> posMap = new HashMap();
    List<Integer> pos = new ArrayList();
    List<String[]> info = new ArrayList();
//    Map<String , String[]> info = new HashMap();
    StringBuilder header = new StringBuilder();
    TabixReader br;
    BufferedWriter bw;
    Boolean first = true,end = false;
    String lineOne = null;
    String startPos,endPos;
    public vcf(){
        
    }
    public void initialVCFread(String inFile){
        try {
            br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("##")){
                    header.append(temp);
                    header.append("\n");
                }else if(temp.startsWith("#C")){
                    header.append(temp);
                    header.append("\n");
                    te = temp.split("\t");
                    sampleSize = te.length - 9;
                    sample = new String[sampleSize];
                    for(int i = 0; i < sampleSize; i++){
                        sample[i] = te[i+9];
                    }
                    break;
                }else{
                    System.out.println("ERROR: Read vcf header error!");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
    public void readVCFByChr(){
        genotype.clear();
        String temp = "";
        info.clear();
        snpNum = 0;
        int i = 0;
        try {
            String[] te = null;
            if(this.lineOne != null){
                te = lineOne.split("\t");
                String[] ge = new String[sampleSize];
                String[] information = new String[9];
                for(int j = 0; j<9;j++){
                    information[j] = te[j];
                }
                info.add(information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                genotype.add(ge);
                snpNum++;
                pos.add(Integer.parseInt(te[1]));
                posMap.put(Integer.parseInt(te[1]),i);
                i++;
            }
            while((temp = this.br.readLine())!=null){
                te = temp.split("\t");
                if(first) {
                    first = false;
                    chr.add(te[0]) ;
                    chrom =te[0];
                }
                if(chr.add(te[0])){
                    lineOne = temp;
                    break;
                }
                chrom =te[0];
                snpNum++;
                String[] ge = new String[sampleSize];
                String[] information = new String[9];
                for(int j = 0; j<9;j++){
                    information[j] = te[j];
                }
                info.add(information);
                for(int j = 0; j< sampleSize;j++){
                    ge[j] = te[j+9];
                }
                genotype.add(ge);
                pos.add(Integer.parseInt(te[1]));
                posMap.put(Integer.parseInt(te[1]), i);
                i++;
            }
            if(temp==null) this.end = true;
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void readVCFBlock(int windowSize){
        genotype.clear();
        String temp = "";
        info.clear();
        int i=0;
        try {
            String[] te = null;
            if(this.lineOne != null){
                snpNum++;
                te = lineOne.split("\t");
                pos.add(Integer.parseInt(te[1]));
                posMap.put(Integer.parseInt(te[1]), i);
                String[] ge = new String[sampleSize];
                String[] information = new String[9];
                for(int j = 0; j<9;j++){
                    information[j] = te[j];
                }
                info.add(information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                genotype.add(ge);
                startPos = te[1];
                i++;
            }
            while((temp = this.br.readLine())!=null){
                te = temp.split("\t");
                if(first) {
                    first = false;
                    chr.add(te[0]) ;
                    chrom =te[0];
                    startPos = te[1];
                }
                if(chr.add(te[0])){
                    lineOne = temp;
                    break;
                }
                pos.add(Integer.parseInt(te[1]));
                posMap.put(Integer.parseInt(te[1]),i);
                i++;
                snpNum++;
                endPos = te[1];
                if((Integer.parseInt(endPos) - Integer.parseInt(startPos)) > windowSize){
                    this.lineOne = temp;
                    break;
                }
                String[] ge = new String[this.sampleSize];
                String[] information = new String[9];
                for(int j = 0; j<9;j++){
                    information[j] = te[j];
                }
                info.add(information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                genotype.add(ge);
            }
            if(temp==null) this.end = true;
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void initialVCFwrite(String outFile){
        try {
            this.bw = IOUtils.getTextWriter(outFile);
            this.bw.write(this.header.toString());
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
    public void writeVCF(){
        try {
            for(int i = 0; i< genotype.size();i++){
                for(int j = 0;j<9;j++){
                    bw.write(info.get(i)[j]+"\t");
                }
                for(int j = 0; j< sampleSize-1;j++){
                    bw.write(genotype.get(i)[j]+"\t");
                }
                bw.write(genotype.get(i)[sampleSize-1]);
                bw.newLine();
                bw.flush();
            }
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void writeVCF(int pos){
        try{
            for(int j = 0;j<9;j++){
                bw.write(info.get(pos)[j]+"\t");
            }
            for(int j = 0; j< sampleSize-1;j++){
                bw.write(genotype.get(pos)[j]+"\t");
            }
            bw.write(genotype.get(pos)[sampleSize-1]);
            bw.newLine();
            bw.flush();
        }catch (Exception e){
            
        }        
    }
    public void close(){
        try {
            this.bw.flush();
            this.bw.close();
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private String toString(String[] a){
        StringBuilder re = new StringBuilder();
        for(int i = 0; i< a.length-1;i++){
            re.append(a[i]);
            re.append("\t");
        }
        re.append(a[a.length-1]);
        return re.toString();
    }
    public boolean checkEnd(){
        return this.end;
    }
    public String getChr(){
        return this.chrom;
    }
    
    public int getSampleSize(){
        return this.sampleSize;
    }
    public int getSNPnum(){
        return this.snpNum;
    }
    public String[] getSample(){
        return this.sample;
    }
    public String getStartPos(){
        return startPos;
    }
    public String getendPos(){
        return endPos;
    }
    public List<String[]> getGenotype(){
        return genotype;
    }
    public Map<Integer,Integer> getMap(){
        return posMap;
    }
    public List<String[]> getInfo(){
        return info;
    }
}
