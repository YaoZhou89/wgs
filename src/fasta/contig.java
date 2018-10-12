/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fasta;

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
public class contig {
    public contig(){
        
    }
    
    public static void splitContig(String inFile, String outFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            StringBuilder b = new StringBuilder();
            String temp = "";
            boolean add = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    System.out.println(temp);
                }
                if(temp.startsWith(">26")){
                    add = true;
                    continue;
                }
                if(!add) continue;
                if(temp.startsWith(">")){
                    break;
                }
                b.append(temp);
            }
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            
            String chr = b.toString();
            String[] contig = chr.split("N");
            bw.write("chr\tcontigNum\tbegin\tend\n");
            int begin = 0;
            int end = 0;
            int c = 0;
            for(int i = 0; i < contig.length; i++){
                
                if(contig[i].length()>0){
                   c++;
                   begin = end+1;
                   end = begin+contig[i].length()-1;
                   bw.write("26\t"+c+"\t"+begin+"\t"+end+"\n");
                }else{
                    end++;
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(contig.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
