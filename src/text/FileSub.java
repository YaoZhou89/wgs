/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package text;

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
public class FileSub {
    public FileSub(){
        
    }
    public void getS(String inFile,String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            while ((temp = br.readLine())!=null){
                if(temp.contains("OSY")||temp.contains("MAY")||temp.contains("BAR")){
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(FileSub.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getSet(String inFile, String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")) continue;
                te = temp.split(" ");
                if(te.length < 2) te = temp.split("\t");
                for (int i = 0; i < te.length; i++){
                    if(te[i].startsWith("Ta")){
                        bw.write(te[i]);
                        bw.newLine();
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(FileSub.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
