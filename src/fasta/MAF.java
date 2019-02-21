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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class MAF {
    public MAF(){
        
    }
    public void getSeq(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            boolean write = false;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">Solanum_lycopersicum")){
                    write = true;
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                if(write){
                    String t1 = temp.replaceAll("-", "");
                    bw.write(t1);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(MAF.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
