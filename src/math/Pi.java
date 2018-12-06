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
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class Pi {
    public Pi(){
        
    }
    
    public void getBedPi(String piFile, String bedFile,String outFile){
        BufferedReader pif = IOUtils.getTextReader(piFile);
        BufferedReader bedf= IOUtils.getTextReader(bedFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        List < Integer[]> bed = new ArrayList();
        String temp = "";
        String[] te = new String[3];
        
        System.out.println("Estimating the mean Pi..");
        try {
            while ((temp = bedf.readLine())!=null){
                Integer[] b = new Integer[3];
                te = temp.split("\t");
                for (int i = 0; i < 3; i++){
                    b[i] = Integer.parseInt(te[i]);
                }
                bed.add(b);
            }
        } catch (IOException ex) {
            Logger.getLogger(Pi.class.getName()).log(Level.SEVERE, null, ex);
        }
        int scanPos = 0;
        int scanPosP = 0;
//        List< Double > meanPi = new ArrayList();
        double sum = 0;
        boolean started = false;
//        System.out.println("Bed file readed, size: "+bed.size()+" Estimating...");
//        for (int i = 0 ;i < 10; i++){
//            System.out.println(bed.get(i)[0]+"\t"+bed.get(i)[1]+"\t"+bed.get(i)[2]);
//        }
        try {
            while ((temp = pif.readLine())!=null){
                if(temp.startsWith("C")) continue;
                te = temp.split("\t");
                int p = Integer.parseInt(te[1]);
//                System.out.println(p);
                scanPosP = scanPos;
                if(te[2].equals("-nan")) te[2] = "0";
                for (int i = scanPosP; i < bed.size(); i++){
//                    System.out.println(bed.get(i)[1]+"\t"+i);
                    if (p < bed.get(i)[1]) {
                        break;
                    }else if (p < bed.get(i)[2]){
//                        System.out.println(p);
                        sum += Double.parseDouble(te[2]);
                        started = true;
                        break;
                    } else {
                        int a = bed.get(i)[1];
                        int c = bed.get(i)[2];
                        bw.write(bed.get(i)[0] + "\t" + a + "\t" + c + "\t" + sum/(c-a));
                        bw.newLine();
                        started = false;
                        sum = Double.parseDouble(te[2]);
                        scanPos++;
                    }
                }
//                if(p < bed.get(scanPos)[1]){
//                    continue;
//                }else if( p < bed.get(scanPos)[2] ){
//                    sum += Double.parseDouble(te[2]);
//                } else {
//                    int a = bed.get(scanPos)[1];
//                    int c = bed.get(scanPos)[2];
//                    bw.write(bed.get(scanPos)[0]+"\t"+ a +"\t"+ c +"\t"+sum/(c-a));
//                    bw.newLine();
//                    sum = Double.parseDouble(te[2]);
//                    scanPos++;
//                }
//                if(scanPos > bed.size()) break;
            }
            if(scanPos == bed.size() -1){
                int a = bed.get(scanPos)[1];
                int c = bed.get(scanPos)[2];
                bw.write(bed.get(scanPos)[0]+"\t"+ a +"\t"+ c +"\t"+sum/(c-a));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Pi.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    public void getSitePi(String inFile,String sitePi,String outFile){
        BufferedReader br = IOUtils.getTextReader(inFile);
        BufferedWriter bw = IOUtils.getTextWriter(outFile);
        BufferedReader bs = IOUtils.getTextReader(sitePi);
        String temp = "";
        String[] te = null;
        String t1 = "";
        Set<String> sites = new HashSet();
        try {
            while ((temp = bs.readLine())!=null){
                String pos = temp.split("_")[1];
                sites.add(pos);
            }
        } catch (IOException ex) {
            Logger.getLogger(Pi.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        try {
            while((temp = br.readLine())!=null){
                te = temp.split("\t");
                if(!sites.add(te[1])){
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(Pi.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
