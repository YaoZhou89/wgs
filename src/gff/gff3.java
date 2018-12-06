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
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class gff3 {
    public gff3(){
        
    }
    
    public void toBed(String inFile,String geneFile, String outFile,String pattern){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        BufferedReader bg = IOUtils.getTextReader(geneFile);
        String temp = "";
        String[] te = null;
        Set patternSet = new HashSet();
        try {
            while((temp = bg.readLine())!=null){
               patternSet.add(temp);
            }
        } catch (IOException ex) {
            Logger.getLogger(gff3.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    continue;
                }
                te = temp.split("\t");
                if(te.length>8){
                    if(te[2].equals(pattern) ){
                        String p = te[8].split(";")[0].split("=")[1];
                        if(!patternSet.add(p)){
                            bw.write(te[0]+"\t"+ (Integer.parseInt(te[3])-1) +"\t" + te[4] +"\t" + p);
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(gff3.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public List<String[]> initialGff3(String inFile){
        BufferedReader br = IOUtils.getTextReader(inFile);
        String temp = "";
        String[] te = null;
        Set ID = new HashSet();
        Set chr = new HashSet();
        List<String[]> cds = new ArrayList();
        String[] pos = new String[2];
        StringBuilder p = new StringBuilder();
        boolean first = true;
        try {
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    continue;
                }
                te = temp.split("\t");
                if(te[2].equals("CDS")){
                    String name =  te[8].split(".")[0].split("=")[1];
                    if(chr.add(te[0])){
                        first = true;
                        ID.add(name);
                        pos[0] = te[0];
                        p.append((Integer.parseInt(te[3])-1) + "_" + (Integer.parseInt(te[4])-1));
                    }else{
                        if(ID.add(name)){
                            if(!first) {
                                cds.add(pos);
                                first = false;
                            }
                            p.append(","+(Integer.parseInt(te[3])-1) + "_" + (Integer.parseInt(te[4])-1));
                        }else{
                            p.append(";"+(Integer.parseInt(te[3])-1) + "_" + (Integer.parseInt(te[4])-1));
                        }
                    }
                }
            }
            cds.add(pos);
        } catch (IOException ex) {
            Logger.getLogger(gff3.class.getName()).log(Level.SEVERE, null, ex);
        }
        return cds;
    }
}
