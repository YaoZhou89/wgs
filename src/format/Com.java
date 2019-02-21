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
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class Com {
    public Com(){
        
    }
    public void getSame(String inFile,String inFile2, String outFile){
        try {
            BufferedReader br1 = IOUtils.getTextReader(inFile);
            BufferedReader br2 = IOUtils.getTextGzipReader(inFile2);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            Set pos = new HashSet();
            while((temp = br1.readLine())!=null){
                te = temp.split("\t");
                if(te.length < 2) continue;
                pos.add(te[1]);
            }
            System.out.println("Size of file1 is "+ pos.size());
            while((temp = br2.readLine())!=null){
                if(temp.startsWith("#")) continue;
                te = temp.split("\t");
                if(!pos.add(te[1])){
                    bw.write(te[0] + "\t" + te[1]);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Com.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
