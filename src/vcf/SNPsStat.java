/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class SNPsStat {
    public SNPsStat(String inFile, String outFile){
        this.statSNP(inFile,outFile);
    }

    private void statSNP(String inFile, String outFile) {
        BufferedReader br ;
        BufferedWriter bw;
        if(inFile.endsWith("gz")){
            br = IOUtils.getTextGzipReader(inFile);
            bw = IOUtils.getTextGzipWriter(outFile);
        }
        else {
            br = IOUtils.getTextReader(inFile);
            bw = IOUtils.getTextWriter(outFile);
        }
        String temp = null;
        String chr = "chr";
        int i = 0;
        String[] ts = null;
        try {
            while((temp = br.readLine())!=null){
               if(!temp.startsWith("#")) {
                   ts = temp.split("\t");
                   if (chr.equals(ts[0])){
                       i++;
                   }else{
                       bw.write(chr + "\t" + i+"\n");
                       bw.flush();
                       chr = ts[0];
                       i = 1;
                   }
               }
            }
            bw.write(chr + "\t" + i+"\n");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(SNPsStat.class.getName()).log(Level.SEVERE, null, ex);
        }
       
    }
}
