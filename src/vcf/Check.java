/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import htsjdk.tribble.readers.TabixReader;
import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
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
    public void getTwo(String inFile){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            String[] te = null;
            int[] a = new int[4];
            
            int b = 0;
            int c = 0;
            int snp = 0;
            while ((temp = br.readLine())!=null){
               if(!temp.startsWith("#")) {
                    snp++;
                    for (int i = 0; i<4 ;i++){
                        a[i] = 0;
                    }
                    te = temp.split("\t");
                    for (int i = 0; i < te.length-9;i++){
                        if(te[i+9].startsWith("0/0")){
                            a[0]++;
                        }else if (te[i+9].startsWith("0/1")){
                            a[1]++;
                        }else if (te[i+9].startsWith("1/1")){
                            a[2]++;
                        }else{
                            a[3]++;
                        }
                    }
                    for (int i = 0; i< 4 ; i++){
                        if(a[1]>1) b++;
                    }
               }
            }
            System.out.println("total SNP is:\t" +snp);
            System.out.println("Number of het SNP in more than two individuals: "+b);
            System.out.println("Minor allele in homozygous individuals: "+c);

        } catch (IOException ex) {
            Logger.getLogger(Check.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void checkAll(String inFile,String outFile){
        try {
            TabixReader br = new TabixReader(inFile);
            String temp = "";
            String[] te = null;
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            int a = 0, b = 0, outlier = 0;
            while((temp = br.readLine()) != null){
                if(temp.startsWith("#")){
                    bw.write(temp);
                    bw.newLine();
                }else{
                    te = temp.split("\t");
                    bw.write(te[0]);
                    for(int i = 1, len = te.length; i < len; i++ ){
                        if(te[i].startsWith("0/0")){
                            String[] t = te[i].split(":");
                            String[] tt = t[1].split(",");
                            a = Integer.parseInt(tt[0]);
                            b = Integer.parseInt(tt[1]);
                            if( a < b){
                                outlier++;
                                if(a > 0){
                                    String ttt = te[i].replace("0/0", "0/1");
                                    bw.write("\t"+ttt);
                                }else{
                                    String ttt = te[i].replace("0/0", "1/1");
                                    bw.write("\t"+ttt);
                                }
                            }else{
                                bw.write("\t" + te[i]);
                            }
                        }else if(te[i].startsWith("0/1")){
                            bw.write("\t" + te[i]);
                        }else if(te[i].startsWith("1/1")){
                            String[] t = te[i].split(":");
                            String[] tt = t[1].split(",");
                            a = Integer.parseInt(tt[0]);
                            b = Integer.parseInt(tt[1]);
                            if( b < a){
                                outlier++;
                                if(b > 0){
                                    String ttt = te[i].replace("1/1", "0/1");
                                    bw.write("\t"+ttt);
                                }else {
                                    String ttt = te[i].replace("1/1", "0/0");
                                    bw.write("\t"+ttt);
                                }
                            }else{
                                bw.write("\t" + te[i]);
                            }
                        }else if(te[i].startsWith("./.")){
                            String[] t = te[i].split(":");
                            if(t.length > 2){
                                String[] tt = t[1].split(",");
                                a = Integer.parseInt(tt[0]);
                                b = Integer.parseInt(tt[1]);
                                if(a == 0){
                                    if( b == 0 ){
                                        bw.write("\t" + te[i]);
                                    }else{
                                        bw.write("\t" + "1/1" + ":" + a + "," + b);
                                        outlier++;
                                    }
                                }else{
                                    if(b > 0){
                                        bw.write("\t" + "0/1" + ":" + a + "," + b);
                                        outlier++;
                                    }else{
                                        bw.write("\t" + "0/0" + ":" + a + "," + b);
                                        outlier++;
                                    }
                                }
                            }else{
                                bw.write("\t" + te[i]);
                            }
                            
                        }else{
                            bw.write("\t" + te[i]);
                            
                        }
                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            System.out.println("Outiler number is: "+outlier);
        } catch (IOException ex) {
            Logger.getLogger(Check.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
