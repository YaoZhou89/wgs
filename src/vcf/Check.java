/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

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
public class Check {
    public Check(){
        
    }
    public static void checkLastLine(String inFile,String suffix){
        File in = new File(inFile);
        String outFile = inFile+"/check.err";
        File[] all = IOUtils.listRecursiveFiles(in);
        File[] subF = IOUtils.listFilesEndsWith(all, suffix);
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            for (File f : subF){
                BufferedReader br = IOUtils.getTextReader(f.toString());
                String temp = "";
                String[] te = null;
                int lineNum = 0;
                boolean first = true;
                while((temp = br.readLine())!=null){
                    if(first && !temp.startsWith("#")){
                        te = temp.split("\t");
                        lineNum = te.length;
                        first =false;
                    }
                    if(!first){
                        te = temp.split("\t");
                        if(te.length != lineNum && te.length > 2){
                            bw.write(f.toString()+"\t"+te[0]+"\t"+te[1]);
                            bw.newLine();
                        }
                    }
                }
                br.close();
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
    public static void addST(String inFile,String outFile,String pos){
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            String[] te = null;
            int lineNum = 0;
            boolean write = false;
            while((temp = br.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    if(te[1].equals(pos)) write = true;
                }
                if(write){
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            
        }
    }
     public static void deleteST(String inFile,String outFile,String pos){
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            String[] te = null;
            int lineNum = 0;
            boolean write = true;
            while((temp = br.readLine())!=null){
                if(!temp.startsWith("#")){
                    te = temp.split("\t");
                    if(te[1].equals(pos)) write = false;
                }
                if(write){
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            
        }
    }
    
}
