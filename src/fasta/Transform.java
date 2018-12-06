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
import org.apache.commons.lang.StringUtils;

/**
 *
 * @author yaozhou
 */
public class Transform {
    public Transform( String inFile, String outFile){
        this.toFastq(inFile, outFile);
    }
    public void toFastq(String inFile, String outFile){
        try {
            StringBuilder fas = null;
            String temp = "";
            BufferedReader br = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">")) {
                    if(fas!=null){
                       generateFastq(fas.toString(),bw);
                    }
                    fas = new StringBuilder();
                    continue;
                };
                fas.append(temp);
            }
            generateFastq(fas.toString(),bw);
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Transform.class.getName()).log(Level.SEVERE, null, ex);
        }
            
    }
    public void generateFastq(String fa,BufferedWriter bw){
        String line1 = "@EAS139:136:FC706VJ:2:2104:15343:197393 1:Y:18:ATCACG";
        String line3 = "+";
        String line4 = StringUtils.repeat("|", 125);
        int len = fa.length();
//        if(len<300) {
//            System.out.println("Warning: Input is two short.");
//        }
        String[] contig = fa.split("N");
        for (int i = 0; i < contig.length; i++){
            String f = contig[i];
            if(f.length()<125) continue;
            for(int j = 0 ; j < f.length()-125;j = j + 12){
                if((j+125) > f.length()) j = f.length()-126;
                try {
                    bw.write(line1);
                    bw.newLine();
//                    System.out.println(i);
                    bw.write(f.substring(j, j+125));
                    bw.newLine();
                    bw.write(line3);
                    bw.newLine();
                    bw.write(line4);
                    bw.newLine();
                } catch (IOException ex) {
                    Logger.getLogger(Transform.class.getName()).log(Level.SEVERE, null, ex);
                }
            }
        }
    }
    public String getReverCom(String dna){
        String out = "";
        for(int i = dna.length() - 1; i >= 0; --i) {
            char curr = dna.charAt(i);
            if(curr == 'A')
                out += 'T';
            else if(curr == 'T')
                out += 'A';
            else if(curr == 'C')
                out += 'G';
            else if(curr == 'G')
                out += 'C';
            else {
                System.out.println("ERROR: Input is not a DNA Sequence.");
                System.exit(-1);
            }
        }
        return out;
    }
    public String[] getReverse(String[] a){
        String[] b = null;
        
        return b;
    }
    public String getComplementary(String dna){
        String out = "";
        for(int i = 0; i < dna.length(); i++) {
            char curr = dna.charAt(i);
            if(curr == 'A')
                out += 'T';
            else if(curr == 'T')
                out += 'A';
            else if(curr == 'C')
                out += 'G';
            else if(curr == 'G')
                out += 'C';
            else {
                System.out.println("ERROR: Input is not a DNA Sequence.");
                System.exit(-1);
            }
        }
        return out;
    }
}
