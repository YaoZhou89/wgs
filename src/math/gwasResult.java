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
public class gwasResult {
    public gwasResult(){}
    public void thinResult(String inFile, String outFile, double rate, double threshold){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String temp = "";
        String[] te = null;
        boolean w = false;
        int k = 0;
        try {
            while ((temp = br.readLine())!=null){
                k++;
                te = temp.split(" ");
                String[] info = new String[9];
                if (temp.startsWith("C")){
                    w = true;
                }else{
                    
                    int j = 0;
                    for(int i = 0; i < te.length;i++){
                        if(!te[i].equals("")){
                            info[j] = te[i];
                            j++;
                        }
                    }
                    info[1] = "rs" + info[0]+"_"+info[2];
                     try{
                        double p = Double.parseDouble(info[8]);
                        if(p < threshold){
                            w = true;
                        }else
                        if(p > threshold){
                            double r = Math.random();
                            if(r < rate) w = true;
                        }
                    }catch (Exception e){
                        w = false;
                    }
                    
                }
                if(w){
                    bw.write(info[1]+"\t"+info[0]+"\t"+info[2]+"\t"+info[8]+"\t"+info[6]+"\n");
                    w = false;
                }
//                if(k < 5){
//                    System.out.println(info[0]+"\t"+info[1]+"\t"+info[2]+"\t"+info[3]);
////                    System.out.println(info.length);
//                }
                
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(gwasResult.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void thinEMMAX (String inFile,String outFile, double rate, double threshold){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            String temp = "";
            String[] te = null;
            boolean write = false;
            while ((temp = br.readLine())!=null){
                te = temp.split("\t");
                double p = Double.parseDouble(te[2]);
                if(p < threshold) {
                    write = true;
                }else{
                    double r = Math.random();
                    if(r < rate) write = true;
                }
                if(write){
                    bw.write(temp);
                    bw.newLine();
                    write = false;
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(gwasResult.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
