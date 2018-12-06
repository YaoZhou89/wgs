/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fasta;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 *
 * @author yaozhou
 */
public class CodonTable {
//    Map <String, Integer> codonTable = new HashMap();;
    public CodonTable(){
        
    }
    public static Map<String, Integer> to4DMap(){
        Map<String, Integer> codonTable = new HashMap();
        codonTable.put("TTT", 2);
        codonTable.put("TTC", 2);
        codonTable.put("TTA", 2);
        codonTable.put("TTG", 2);
        codonTable.put("TCT", 4);
        codonTable.put("TCC", 4);
        codonTable.put("TCA", 4);
        codonTable.put("TCG", 4);
        codonTable.put("TAT", 2);
        codonTable.put("TAC", 2);
        codonTable.put("TAA", 1);
        codonTable.put("TAG", 1);
        codonTable.put("TGT", 2);
        codonTable.put("TGC", 2);
        codonTable.put("TGA", 1);
        codonTable.put("TGG", 1);
        
        codonTable.put("CTT", 4);
        codonTable.put("CTC", 4);
        codonTable.put("CTA", 4);
        codonTable.put("CTG", 4);
        codonTable.put("CCT", 4);
        codonTable.put("CCC", 4);
        codonTable.put("CCA", 4);
        codonTable.put("CCG", 4);
        codonTable.put("CAT", 2);
        codonTable.put("CAC", 2);
        codonTable.put("CAA", 2);
        codonTable.put("CAG", 2);
        codonTable.put("CGT", 4);
        codonTable.put("CGC", 4);
        codonTable.put("CGA", 4);
        codonTable.put("CGG", 4);
        
        codonTable.put("ATT", 3);
        codonTable.put("ATC", 3);
        codonTable.put("ATA", 3);
        codonTable.put("ATG", 1);
        codonTable.put("ACT", 4);
        codonTable.put("ACC", 4);
        codonTable.put("ACA", 4);
        codonTable.put("ACG", 4);
        codonTable.put("AAT", 2);
        codonTable.put("AAC", 2);
        codonTable.put("AAA", 2);
        codonTable.put("AAG", 2);
        codonTable.put("AGT", 2);
        codonTable.put("AGC", 2);
        codonTable.put("AGA", 2);
        codonTable.put("AGG", 2);
        
        codonTable.put("GTT", 4);
        codonTable.put("GTC", 4);
        codonTable.put("GTA", 4);
        codonTable.put("GTG", 4);
        codonTable.put("GCT", 4);
        codonTable.put("GCC", 4);
        codonTable.put("GCA", 4);
        codonTable.put("GCG", 4);
        codonTable.put("GAT", 2);
        codonTable.put("GAC", 2);
        codonTable.put("GAA", 2);
        codonTable.put("GAG", 2);
        codonTable.put("GGT", 4);
        codonTable.put("GGC", 4);
        codonTable.put("GGA", 4);
        codonTable.put("GGG", 4);
        return codonTable;
    }
}
