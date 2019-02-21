/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class MergeVmapI {
    public void MergeVampI(){
        
    }
    public void mergeBarley(String inFile,String barleyFile, String outFile){
        try {
            BufferedReader br = null;
            BufferedReader bb = null;
            if(inFile.endsWith(".gz")){
                br = IOUtils.getTextGzipReader(inFile);
            }else{
                br = IOUtils.getBufferedReader(inFile);
            }
            if(inFile.endsWith(".gz")){
                bb = IOUtils.getTextGzipReader(barleyFile);
            }else{
                bb = IOUtils.getBufferedReader(barleyFile);
            }
            String temp = "";
            Set barley = new HashSet();
            String[] pos = null;
            while ((temp = bb.readLine())!=null){
                pos = temp.split("\t");
                System.out.println(pos[0] );
                break;
            }
        } catch (IOException ex) {
            Logger.getLogger(MergeVmapI.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
