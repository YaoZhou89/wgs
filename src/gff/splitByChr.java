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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class splitByChr {
    public splitByChr(String inFile){
        this.getGtf(inFile);
    }
    public void getGtf(String inFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            String chr ="-1";
            BufferedWriter bw = null ;
            String outFile="";
            boolean first = true;
            while((temp = br.readLine())!=null){
                String[] tem = temp.split("\t");
                if(first){
                    chr = tem[0];
                }
                if(tem[0].equals(chr) && !first){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    
                    chr = tem[0];
                    if(!first){
                        bw.flush();
                        bw.close();
                    }
                    first = false;
                    String name = ".chr"+String.format("%03d", Integer.parseInt(chr))+".gtf";
                    outFile = inFile.replace(".gtf",name);
                    bw = IOUtils.getTextWriter(outFile);
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
}
