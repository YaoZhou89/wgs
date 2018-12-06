/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package gff;

import fasta.CodonTable;
import fasta.Genome;
import java.util.List;
import java.util.Map;

/**
 *
 * @author yaozhou
 */
public class FourDegenerateSite {
    public FourDegenerateSite(){
        
    }
    public void get4D(String fa,String gff,String outFile){
        
        Map<String, Integer> codonTable = CodonTable.to4DMap();
        List<String[]> cds = new gff3().initialGff3(gff);
        Genome genome = new Genome();
        genome.initialGenome(fa);
        for(String[] a : cds){
            String ref = genome.getSeq(a[0]);
            String[] b = a[1].split(";");
            for (String p : b){
                
            }
        }
    }
    
}
