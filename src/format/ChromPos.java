/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class ChromPos {
    public ChromPos(){
        
    }
    Map<String, String[]> pos = new HashMap();
    private void readPos(String inFile){
        try {
            String temp = "";
            String[] te = null;
            BufferedReader br = IOUtils.getTextReader(inFile);
            if((temp = br.readLine())!=null){
                te = temp.split("\t");
                pos.put(te[0],te);
            }
        } catch (IOException ex) {
            Logger.getLogger(ChromPos.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void changePos(String inFile, String inFile2, String outFile,int pos){
        try {
            String temp = "";
            String[] te = null;
            BufferedReader br = null;
            
            if(inFile.endsWith(".gz")){
                br = IOUtils.getTextGzipReader(inFile);
            }else {
                br = IOUtils.getTextReader(inFile);
            }
            
            readPos(inFile2);
            
            while ((temp = br.readLine())!=null){
                
            }
        } catch (IOException ex) {
            Logger.getLogger(ChromPos.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
