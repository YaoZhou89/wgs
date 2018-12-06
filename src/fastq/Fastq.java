/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fastq;

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
public class Fastq {
    public Fastq(){
        
    }
    
    public void changeName(String inFile, String outFile){
        try {
            BufferedReader br = null;
            BufferedWriter bw = IOUtils.getTextGzipWriter(outFile);
            if(inFile.endsWith(".gz")){
                br = IOUtils.getTextGzipReader(inFile);
            }else{
                br = IOUtils.getTextReader(inFile);
            }
            String temp = "";
            Integer n = 0;
            int a = 0;
            while((temp = br.readLine())!=null){
                if(a%4 == 0){
                    n++;
                    bw.write("@"+n);
                    bw.newLine();
                }else{
                    bw.write(temp);
                    bw.newLine();
                }
                a++;
            }
            System.out.println("Total sequence is: " + n);
        } catch (IOException ex) {
            Logger.getLogger(Fastq.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
