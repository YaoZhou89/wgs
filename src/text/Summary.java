/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package text;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class Summary {
    public Summary(){
        
    }
    
    public void readData(String inFile, String suffix){
        File a = new File (inFile);
        File[] as = IOUtils.listRecursiveFiles(a);
        File[] subFiles = IOUtils.listFilesEndsWith(as, suffix);
        BufferedReader br ;
        String outFile  = inFile + "/summary.txt";
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        Integer[][] sum = new Integer[4][2];
        for (int i = 0; i<4;i++){
            for (int j = 0; j<2;j++){
                sum[i][j] = 0;
            }
        }
        String temp = "";
        String[] te = null;
        try {
            for (File f : subFiles){
                br = IOUtils.getTextReader(f.toString());
                int i = 0;
                while((temp = br.readLine())!=null){
                    te = temp.split("\t");
                    sum[i][0] += Integer.parseInt(te[0]);
                    sum[i][1] += Integer.parseInt(te[1]);
                    i++;
                }
            }
            bw.write("Total SNPs"+"\tFailed SNPs\t"+"Passed SNPs\n");
            for (int i = 0; i<4;i++){
                bw.write(sum[i][0]+"\t"+sum[i][1]+"\t"+(sum[i][0]-sum[i][1]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
        } catch (IOException ex) {
            ;
        }
    }
}
