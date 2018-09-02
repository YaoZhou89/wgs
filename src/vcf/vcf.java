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
    List<String> ID = new ArrayList();
    int sampleSize;
    String[] sample;
    Map<String , String[]> genotype = new HashMap();
    Map<String , String[]> info = new HashMap();
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
        ID.clear();
        info.clear();
        
        try {
            String[] te = null;
            int i = 0;
            String pos = "";
            
            if(this.lineOne != null){
                i++;
                this.snpNum++;
                te = this.lineOne.split("\t");
                pos = te[0]+"_"+te[1];
                this.ID.add(pos);
                String[] ge = new String[sampleSize];
                String[] information = new String[7];
                for(int j = 0; j<7;j++){
                    information[j] = te[j+2];
                }
                this.info.put(pos,information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                this.genotype.put(pos, ge);
            }
            while((temp = this.br.readLine())!=null){
                te = temp.split("\t");
                if(this.first) {
                    this.first = false;
                    this.chr.add(te[0]) ;
                    this.chrom =te[0];
                }
                if(this.chr.add(te[0])){
                    this.lineOne = temp;
                    break;
                }
                this.chrom =te[0];
                i++;
                this.snpNum++;
                pos = te[0]+"_"+te[1];
                this.ID.add(pos);
                String[] ge = new String[this.sampleSize];
                String[] information = new String[7];
                for(int j = 0; j<7;j++){
                    information[j] = te[j+2];
                }
                this.info.put(pos,information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                this.genotype.put(pos, ge);
            }
            if(temp==null) this.end = true;
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void readVCFBlock(int windowSize){
        genotype.clear();
        String temp = "";
        ID.clear();
        info.clear();
        try {
            String[] te = null;
            int i = 0;
            String pos = "";
            if(this.lineOne != null){
                i++;
                this.snpNum++;
                te = this.lineOne.split("\t");
                pos = te[0]+"_"+te[1];
                this.ID.add(pos);
                String[] ge = new String[sampleSize];
                String[] information = new String[7];
                for(int j = 0; j<7;j++){
                    information[j] = te[j+2];
                }
                this.info.put(pos,information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                this.genotype.put(pos, ge);
                startPos = te[1];
            }
            while((temp = this.br.readLine())!=null){
                te = temp.split("\t");
                if(this.first) {
                    this.first = false;
                    this.chr.add(te[0]) ;
                    this.chrom =te[0];
                    startPos = te[1];
                }
                if(this.chr.add(te[0])){
                    this.chrom =te[0];
                    this.lineOne = temp;
                    break;
                }
                i++;
                this.snpNum++;
                endPos = te[1];
                if((Integer.parseInt(endPos) - Integer.parseInt(startPos)) > windowSize){
                    this.lineOne = temp;
                    break;
                }
                pos = te[0]+"_"+te[1];
                this.ID.add(pos);
                String[] ge = new String[this.sampleSize];
                String[] information = new String[7];
                for(int j = 0; j<7;j++){
                    information[j] = te[j+2];
                }
                this.info.put(pos,information);
                for(int j = 0; j< this.sampleSize;j++){
                    ge[j] = te[j+9];
                }
                this.genotype.put(pos, ge);
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
            String pos = "";
            for(int i = 0; i<ID.size();i++){
                pos = ID.get(i);
                bw.write(pos.replace("_","\t"));
                bw.write("\t");
                bw.write(toString(info.get(pos)));
                bw.write("\t");
                bw.write(toString(genotype.get(pos)));
                bw.newLine();
                bw.flush();
            }
        } catch (IOException ex) {
            Logger.getLogger(vcf.class.getName()).log(Level.SEVERE, null, ex);
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
    public Map getGeno(){
        return this.genotype;
    }
    public List getPos(){
        return this.ID;
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
    
}
