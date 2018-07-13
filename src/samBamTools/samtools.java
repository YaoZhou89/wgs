/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package samBamTools;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class samtools {
    public samtools(String inFile,String model){
        this.getMAPquality(inFile);
    }
    private void getMAPquality(String inFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = null;
            String[] tem = null;
            int QS[] = new int[201];
            BufferedWriter bw = IOUtils.getTextWriter(inFile+".stat");
            BufferedWriter bs = IOUtils.getTextWriter(inFile+".quality.bed");
            while((temp = br.readLine())!=null){
                if(!temp.startsWith("@")){
                    tem = temp.split("\t");
                    int index = Integer.parseInt(tem[4]);
                    if(index > 200) index = 200;
                    QS[index]++;
                    int end = Integer.parseInt(tem[3])+64;
                    bs.write(tem[2] + "\t" + tem[3] + "\t" + end +"\t" + tem[4]);
                    bs.newLine();
                }
            }
            for (int i =0;i<QS.length;i++){
                bw.write(String.valueOf(i)+"\t");
                bw.write(String.valueOf(QS[i]));
                bw.newLine();
            }
            bs.flush();
            bs.close();
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(samtools.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
