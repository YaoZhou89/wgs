/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package GO;

import io.IOUtils;
import java.io.BufferedReader;
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
public class GO {
    public GO(){
        
    }
    public void randomSelect(String inFile, String outFile){
        
    }
    public void toGO(String inFile,String ref, String outFile){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedReader bf = IOUtils.getTextReader(ref);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        Set GO = new HashSet();
        Map<String, List> G = new HashMap();
        String temp = "";
        String[] te = null;
        try {
            while((temp = br.readLine())!=null){
               GO.add(temp);
            }
        } catch (IOException ex) {
            Logger.getLogger(GO.class.getName()).log(Level.SEVERE, null, ex);
        }
        try {
            while((temp = bf.readLine())!=null){
                if(temp.startsWith("Gene")) continue;
                te = temp.split("\t");
                if(GO.add(te[0])){
                    GO.remove(te[0]);
                }else{
                    for(int i = 0; i < te.length;i++){
                        if(te[i].startsWith("GO")){
                            bw.write(te[0]+"\t"+te[i]+"\n");
                        }
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GO.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
