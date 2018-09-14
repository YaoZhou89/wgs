/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package fasta;

import io.IOUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class Genome {
    public Genome(){
    }
    public void readByChromosome(String inFile,int a, int b){
        try {
            BufferedReader br = IOUtils.getTextReader(inFile);
            String temp = "";
            StringBuilder chr = new StringBuilder();
            boolean test = false, chr1 =true;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">1")){
                    test = true;
                    System.out.println(temp);
                    chr = new StringBuilder();
                }else if (test){
                    chr.append(temp);
                }
                if(temp.startsWith(">2")){
                   break;
                }
            }
            System.out.println(chr.toString().substring(a, b));
        } catch (IOException ex) {
            Logger.getLogger(Genome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
