/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

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
public class Blast {
    public  Blast(){
        
    }
    public void changeBlastp(String inFile,String outFile){
        try {
            String temp = "";
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while ((temp =  br.readLine())!=null){
                String te = temp.replace("gnl|","");
                String[] t = te.split("\t");
                if(t.length != 12) continue;
                bw.write(te);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Blast.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
