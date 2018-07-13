/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

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
public class SegeregationTest {
    public SegeregationTest(String inFile, String outFile){
        this.getST(inFile, outFile);
    }
    private void getST(String inFile, String outFile){
        try {
            BufferedReader br;
            BufferedWriter bw;
            if(inFile.endsWith(".gz")) br = IOUtils.getTextGzipReader(inFile);
            else br = IOUtils.getTextReader(inFile);
            String temp ;
            bw = IOUtils.getTextWriter(outFile);
            while ((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    String[] te = temp.split("\t");
                    int s = te.length;
                    double[][] depth = new double[2][s-9];
                    for (int i = 0; i < te.length-9; i++){
                        if(te[i+9].equals("./.")) {
                            depth[1][i] = Double.NaN;
                            depth[2][i] = Double.NaN;
                        }
                        else{
                            String[] de = te[i+9].split(":")[1].split(",");
                            depth[0][i] = Double.parseDouble(de[0]);
                            depth[1][i] = Double.parseDouble(de[1]);
                            
                        }
                    }
                    ContMatPack dpt = new ContMatPack(depth);
                    
                    double pv = dpt.pchi;
                    if(pv < 0.2){
                        pv = dpt.pvalue(depth);
                        if(pv < 0.01){
                            bw.write(temp);
                            bw.newLine();
                        }
                    }else{
                        continue;
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
