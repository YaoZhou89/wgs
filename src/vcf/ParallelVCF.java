/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.IOUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author yaozhou
 */
public class ParallelVCF {
    // quality:
    // depth: maxdepth 500 mindepth 125 maxSD
    // ST: none
    // LD: --windowSize 50 --threshold 0.1
    // IBD:--windowSize 50 --threshold 0.1 
    
    public ParallelVCF(String inFile, String suffix,String type,String outsuffix,
            String MQ,String FS,String MQRankSum,String ReadPosRankSum,String BSQRankSum,
            String SOR,int maxdepth, int mindepth, double maxSD, int windowSize,
            double threshold){
        
             this.getAll(inFile,suffix,type,outsuffix,MQ,FS,MQRankSum,ReadPosRankSum,
                    BSQRankSum,SOR ,maxdepth,mindepth,maxSD,windowSize,threshold);
        
    }
    public void getAll(String inFile, String suffix,String type,String outsuffix,
            String MQ,String FS,String MQRankSum,String ReadPosRankSum,String BSQRankSum,
            String SOR,int maxdepth, int mindepth, double maxSD, int windowSize,
            double threshold){
        File a = new File(inFile);
        File[] b = IOUtils.listRecursiveFiles(a);
        File[] subFile = IOUtils.listFilesEndsWith(b, suffix);
        List<String> inputs = new ArrayList<String>(subFile.length);
        for(File f : subFile){
            inputs.add(f.toString());
        }
        System.out.println("Paralle analyisizing.....");
        inputs
            .parallelStream()
            .forEach(e -> new VcfTools(e,suffix,outsuffix,MQ,FS,MQRankSum,
                    ReadPosRankSum,BSQRankSum,SOR ,maxdepth,mindepth,maxSD,windowSize,threshold));
    }
    

}
