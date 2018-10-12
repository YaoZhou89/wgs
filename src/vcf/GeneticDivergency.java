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
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class GeneticDivergency {
    public GeneticDivergency(){
        
    }
    public static void calPair(String inFile, String outFile){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        String temp= "";
        String[] te = null;
        boolean header = true;
        String[] names = null;
        int sampleNum = 0;
        Double[][] GD = null; 
        String a = "",b = "";
        try {
            while ((temp = br.readLine())!=null){
                if(!temp.startsWith("##")){
                    te = temp.split("\t");
                    if(temp.startsWith("#")){
                        sampleNum = te.length - 9;
                        names = new String[sampleNum];
                        GD = new Double[sampleNum][sampleNum];
                        for(int i = 0; i < sampleNum ; i++){
                            for (int j = 0; j < sampleNum; j++){
                                GD[i][j] = 0.0;
                            }
                        }
                        bw.write(te[9]);
                        for (int i = 1; i< sampleNum;i++){
                            names[i] = te[i+9];
                            bw.write("\t"+te[i+9]);
                        }
                        bw.flush();
                        bw.newLine();
                        header = false;
                    }else{
                       for(int i = 0; i < sampleNum-1;i++){
                           for (int j = i+1; j< sampleNum;j++){
                                if(!te[i+9].startsWith(".") && !te[j+9].startsWith(".")){
                                    a = te[i+9].split(":")[0];
                                    b = te[j+9].split(":")[0];
                                    GD[i][j] += 1 ;
                                    if(a.equals("0/0")){
                                        if(b.equals("0/1")) {
                                            GD[j][i] += 0.5; 
                                        }else if(b.equals("0/0")){
                                            GD[j][i] += 1; 
                                        }else{
                                            
                                        }
                                    }else if(a.equals("1/1")){
                                        if(b.equals("1/1")) {
                                            GD[j][i]+= 1; 
                                        }else if(b.equals("0/1")){
                                            GD[j][i]+= 0.5; 
                                        }else{
                                            
                                        }
                                    }else{
                                        GD[j][i]+=0.5;
                                    }
                               }
                           }
                        }
                    }
                }
            }
            for(int i = 0; i< sampleNum-1;i++){
                for (int j = i+1; j< sampleNum;j++ ){
                    GD[j][i] = GD[j][i]/GD[i][j];
                    GD[i][j] = GD[j][i];
                }
            }
            for(int i = 0; i< sampleNum;i++){
                for(int j = 0; j< sampleNum-1;j++){
                    bw.write(1-GD[i][j]+"\t");
                }
                bw.write(1-GD[i][sampleNum-1]+"\n");
                bw.flush();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            
        }
    }
    
}
