/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import io.BuilderFromVCF;
import io.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.logging.Level;
import java.util.logging.Logger;
import net.maizegenetics.analysis.distance.*;
import net.maizegenetics.dna.snp.ExportUtils;
import net.maizegenetics.dna.snp.FilterGenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTable;
import net.maizegenetics.dna.snp.GenotypeTableBuilder;
import net.maizegenetics.dna.snp.GenotypeTableUtils;
//import net.maizegenetics.dna.snp.io.BuilderFromVCF;
import net.maizegenetics.taxa.distance.DistanceMatrix;
/**
 *
 * @author yaozhou, modified from bukowski's code
 */
public class ibdfilter {
//    public ibdfilter(){
//        
//    };
    public static void PfilterBasedOnIBD(String hm3gbsSNPsFile, String hm3SNPsFile, String outFile, int minComp, double maxIBDDist, 
        int windowSize, int numThreads, int bestContrasts) {
        Integer[] num = new Integer[2];
        int lowPhysLimit,lowAnchorLimit;
        int highPhysLimit,highAnchorLimit;
        String chr;
        BufferedWriter bw = IOUtils.getTextWriter(outFile+".info.txt");
        BufferedWriter bw1 = IOUtils.getTextWriter(outFile+".pos.txt");
        BufferedWriter bw2 = IOUtils.getTextWriter(outFile+".vcf");
        try {
            bw.write("Window\tStartSite\tTaxa\tmaxIBDDist\t IBDTaxa\tIBDContrasts\tMeanDist\n");
            GenotypeTable hm3gbs=net.maizegenetics.dna.snp.io.BuilderFromVCF.getBuilder(hm3gbsSNPsFile).build();        
            GenotypeTable hm3=net.maizegenetics.dna.snp.io.BuilderFromVCF.getBuilder(hm3SNPsFile).build();
            if(hm3 == null) { System.out.println(" inFile is null...."+hm3SNPsFile); }
            if(hm3gbs == null) { System.out.println(" anchor is null...."+hm3gbsSNPsFile); }
        // try to determine the chromosome automatically
            chr = hm3gbs.chromosomes()[0].getName();
            if(! chr.equals(hm3.chromosomes()[0].getName())){
                System.out.println("Chromosome names not the same: "+chr+", "+hm3.chromosomes()[0].getName());
                return;
            }else{
                System.out.println("Processing chromosome "+chr);
            }
        // Determine the low and high physical limits of the hm3 file
            lowPhysLimit = hm3.physicalPositions()[0];
            highPhysLimit = hm3.physicalPositions()[hm3.numberOfSites()-1];
            System.out.println("Low, High phys limit in hm3: "+ lowPhysLimit+", "+highPhysLimit);
            lowAnchorLimit = hm3gbs.physicalPositions()[0];
//            highAnchorLimit = hm3gbs.physicalPositions()[hm3gbs.numberOfSites()-1];
            String temp = "";
            String te[] = null;
            BufferedReader br = IOUtils.getTextReader(hm3SNPsFile);
            while ((temp = br.readLine())!=null){
                bw2.write(temp);
                bw2.newLine();
                bw2.flush();
                if(!temp.startsWith("#")){
                    if(temp.split("\t")[1].equals(Integer.toString(lowAnchorLimit))) break;
                }
            }
            Set PassPos = new HashSet();
        // Normalize taxa names
            List<String> hm3gbsTaxaNames = new ArrayList<String>();
            for(int i=0;i<hm3gbs.numberOfTaxa(); i++){
                hm3gbsTaxaNames.add(hm3gbs.taxaName(i).replaceAll(".h5", ""));
            }
            List<String> hm3TaxaNames = new ArrayList<String>();
            for(int i=0;i<hm3.numberOfTaxa(); i++){
                hm3TaxaNames.add(hm3.taxaName(i).replaceAll(".h5", ""));
            }
            int[] aTobTaxa=new int[hm3gbs.numberOfTaxa()];
            for (int i=0; i<hm3gbs.numberOfTaxa(); i++) {
                aTobTaxa[i]=Integer.MIN_VALUE;
                int tind = hm3TaxaNames.indexOf(hm3gbsTaxaNames.get(i));
                if(tind>-1){
                    aTobTaxa[i]=tind;
                }
                else
                {
                    System.out.println("hm3gbs taxon "+hm3gbsTaxaNames.get(i)+" has no counterpart in hm3 file");
                }
            }
        //GenotypeTableBuilder gtb=GenotypeTableBuilder.getSiteIncremental(hm3.taxa());
            GenotypeTable a;
            GenotypeTable b;
            int lastHm3SiteProcessed = 0;
            int lastHm3SiteProcPos = 0;
            for (int start=0; start < hm3gbs.numberOfSites(); start+=windowSize) {
            
                int phys_start_gbs = hm3gbs.physicalPositions()[start];
                int phys_end_gbs; 
                int lastGBSsite = start+windowSize-1;
                if(lastGBSsite >= hm3gbs.physicalPositions().length){
                    lastGBSsite = hm3gbs.physicalPositions().length-1;
                }
                phys_end_gbs = hm3gbs.physicalPositions()[lastGBSsite];
                if(phys_end_gbs < lowPhysLimit || phys_start_gbs > highPhysLimit) { continue; }
            
                int endInstance = phys_end_gbs;
                int startInstance = phys_start_gbs;
                if(endInstance > highPhysLimit){
                    endInstance = highPhysLimit;
                }
                if(startInstance<lowPhysLimit){
                    startInstance = lowPhysLimit;
                }
                System.out.println("phys_start_gbs: "+phys_start_gbs + " phys_end_gbs: "+phys_end_gbs);
                System.out.println("startInstance: "+startInstance + " endInstance: "+endInstance);
                
            //GenotypeTable a=hm3gbs;
                a=FilterGenotypeTable.getInstance(hm3gbs,start,lastGBSsite);
            //GenotypeTable b=hm3;
                b=FilterGenotypeTable.getInstance(hm3, chr, startInstance, endInstance);
                if(b == null) { continue; }
                System.out.println("Taxa in b: " + b.numberOfTaxa());
                System.out.println("Sites in b: " + b.numberOfSites());
            
                DistanceMatrix dm = IBSDistanceMatrix.getInstance(a, minComp, null);
                int numTaxa = dm.numberOfTaxa();
                int numTestSites=b.numberOfSites();
                int numTaxaChunk = 150;
            
                int [][] acounts = new int[numTaxaChunk][3]; // countIBD, taxaIBD, count
                double [][] amean = new double[numTaxaChunk][1];
                int[][] amjSame=new int[numTaxaChunk][numTestSites];
                int[][] aminorSame=new int[numTaxaChunk][numTestSites];
                int[][][] arightNuc=new int[numTaxaChunk][6][numTestSites];
                int[][][] awrongNuc=new int[numTaxaChunk][6][numTestSites];
                int[][] adiff=new int[numTaxaChunk][numTestSites];
                int[][] ahet=new int[numTaxaChunk][numTestSites];
                int[][] aminorComp=new int[numTaxaChunk][numTestSites];
                byte[][] diploids=new byte[3][numTestSites];//mj:mj,mj:mn,mn:mn
                String [][] acontrastPairs = new String[numTaxaChunk][1];
                String [][] afailedContrasts = new String[numTaxaChunk][numTestSites];
                for (int k=0; k<numTestSites; k++) {
                    diploids[0][k]=GenotypeTableUtils.getDiploidValue(b.majorAllele(k),b.majorAllele(k));
                    diploids[1][k]=GenotypeTableUtils.getDiploidValue(b.majorAllele(k),b.minorAllele(k));
                    diploids[2][k]=GenotypeTableUtils.getDiploidValue(b.minorAllele(k),b.minorAllele(k));
                }
            // Initialize string arrays holding contrast info
                for(int i=0;i<numTaxaChunk;i++){
                    acontrastPairs[i][0] = "";
                    for(int j=0;j<numTestSites;j++){
                        afailedContrasts[i][j] = "";
                    }
                }
            // If requested, collect the distribution of genetic distances - we will then pick threshold to
            // to use 50 or so best contrasts. Do not parallelize this for now...
                double maxIBDDistAdj = maxIBDDist;
                if(bestContrasts > 0){
                    double [] dstmat2sort = new double[numTaxa*(numTaxa-1)/2];
                    int counter = 0;
                    for (int i=0; i<dm.numberOfTaxa(); i++) {
                        for(int j=i+1;j<dm.numberOfTaxa(); j++){
                            if(Double.isNaN(dm.getDistance(i,j))) continue;
                            dstmat2sort[counter] = dm.getDistance(i,j);
                            counter++;
                    }
                }
                // Sort the array (in ascending order) and pick, say 50th value from the end as threshold
                // This will ensure we count only 50 best IBD pairs.
                // NOTE: numTaxa*(numTaxa-1)/2 - counter last elements (first after sorting) will be zero!
                Arrays.sort(dstmat2sort);
                // Debug print
                System.out.println("Sorted dstmat2sort matrix:");
                for(int i=0;i<dstmat2sort.length;i++){
                    System.out.print(dstmat2sort[i]+",");
                }
                System.out.println();
                maxIBDDistAdj = dstmat2sort[Math.min(dstmat2sort.length-1,bestContrasts-1+dstmat2sort.length-counter)];
                maxIBDDistAdj = Math.min(maxIBDDistAdj,maxIBDDist);
            }
            // We will parallelize over i loop
            java.util.concurrent.ExecutorService executor = Executors.newFixedThreadPool(numThreads);
            
            for (int i=0; i<dm.numberOfTaxa(); i++) {
                int ir = i % numTaxaChunk;
                pFilter pflt = new pFilter(b, dm, i, maxIBDDistAdj, aTobTaxa, numTestSites, arightNuc[ir], awrongNuc[ir], 
                        diploids, ahet[ir], aminorComp[ir], aminorSame[ir], amjSame[ir], adiff[ir], acounts[ir], amean[ir],
                        acontrastPairs[ir], afailedContrasts[ir]);
                
                executor.execute(pflt);
                
                // Re-create the thread pool 
                if(ir == numTaxaChunk-1)
                {
                    executor.shutdown();
                    try
                    {
                        executor.awaitTermination(100, TimeUnit.DAYS);
                    }
                    catch(Exception e)
                    {
                        System.out.println("executor termination problem: "+e);
                    }
                    executor = Executors.newFixedThreadPool(numThreads);
                }
            }
            
            executor.shutdown();
            try
            {
                executor.awaitTermination(100, TimeUnit.DAYS);
            }
            catch(Exception e)
            {
                System.out.println("executor termination problem: "+e);
            }
            
            // Collect cummmulatives over threads
            int countIBD=0, taxaIBD=0, count=0;
            double mean=0;
            int[] mjSame=new int[numTestSites];
            int[] minorSame=new int[numTestSites];
            int[][] rightNuc=new int[6][numTestSites];
            int[][] wrongNuc=new int[6][numTestSites];
            int[] diff=new int[numTestSites];
            int[] het=new int[numTestSites];
            int[] minorComp=new int[numTestSites];
            String contrastPairs = "";
            String [] failedContrasts = new String[numTestSites];
            for(int i=0;i<numTestSites;i++)
            {
                failedContrasts[i] = "";
            }
            for(int i=0;i<numTaxaChunk;i++)
            {
                countIBD += acounts[i][0];
                taxaIBD += acounts[i][1];
                count += acounts[i][2];
                mean += amean[i][0];
                contrastPairs += acontrastPairs[i][0];
                //int cntrslen = acontrastPairs[i][0].split("\n").length;
                //System.out.println("cntrst len for i " + i + ":" + cntrslen);
                for(int k=0;k<numTestSites;k++)
                {
                    mjSame[k] += amjSame[i][k];
                    minorSame[k] += aminorSame[i][k];
                    diff[k] += adiff[i][k];
                    het[k] += ahet[i][k];
                    minorComp[k] += aminorComp[i][k];
                    failedContrasts[k] += afailedContrasts[i][k];
                }
                for(int bs=0;bs<6;bs++)
                {
                    for(int k=0;k<numTestSites;k++)
                    {
                        rightNuc[bs][k] += arightNuc[i][bs][k];
                        wrongNuc[bs][k] += awrongNuc[i][bs][k];
                    }   
                }
            }
            
            mean/=(double)count;
//            System.out.printf("Window:%d StartSite:%d Taxa:%d maxIBDDist:%g IBDTaxa:%d IBDContrasts:%d MeanDist:%g%n",
//                    windowSize, start,hm3gbs.numberOfTaxa(),maxIBDDistAdj,taxaIBD,countIBD,mean);
                //System.out.println("Contrast pairs:");
                //System.out.print(contrastPairs);
                bw.write(windowSize+"\t"+ start+"\t"+hm3gbs.numberOfTaxa()+"\t"+maxIBDDistAdj+"\t"+taxaIBD+"\t"+countIBD+"\t"+mean+"\n");
            for (int k=0; k<numTestSites; k++) {
                double minorRate=(double)minorSame[k]/((double)diff[k]+0.5);
                int countBiased=0;
                int counttries = 0;
                for (int bs=0; bs<6; bs++) {
                    if(((double)rightNuc[bs][k]/(0.5+wrongNuc[bs][k])) >2.0) countBiased++;
                    counttries += rightNuc[bs][k] + wrongNuc[bs][k];
                }
                //if(minorComp[k]>0 && minorRate>1) {
                if(countBiased>1 || counttries == 0) {
                    // "Good" site
                    bw1.write(b.positions().get(k).getPosition()+"\t"+ het[k]+"\t"+mjSame[k]+"\t"+minorSame[k]+"\t"+diff[k]+"\t"+minorComp[k]+"\t"+minorRate+"\t");
                    for (int bs=0; bs<6; bs++) {
                        bw1.write(rightNuc[bs][k]+":"+wrongNuc[bs][k]+"\t");
                    }
                    if(counttries > 0)
                    {
                        bw1.write("R\n");
                        if(diff[k]==0 ) PassPos.add(Integer.toString(b.positions().get(k).getPosition()));
                    }else{
                        bw1.write("A\n");
                    }
                    lastHm3SiteProcPos = b.positions().get(k).getPosition();
                }else{
                    bw1.write(b.positions().get(k).getPosition()+"\t"+ het[k]+"\t"+mjSame[k]+"\t"+minorSame[k]+"\t"+diff[k]+"\t"+minorComp[k]+"\t"+minorRate+"\t");
                    for (int bs=0; bs<6; bs++) {
                        bw1.write(rightNuc[bs][k]+":"+wrongNuc[bs][k]+"\t");
                    }
                    bw1.write("W\n");
                    PassPos.add(Integer.toString(b.positions().get(k).getPosition()));
                    lastHm3SiteProcPos = b.positions().get(k).getPosition();
                }
            }
            lastHm3SiteProcessed += numTestSites;
        }
        System.out.println("Last HM3 position processed: "+ lastHm3SiteProcPos);
        while((temp=br.readLine())!=null){
            String pos = temp.split("\t")[1];
            if(!PassPos.add(pos)){
                bw2.write(temp);
                bw2.newLine();
                bw2.flush();
            }
            if(temp.split("\t")[1].equals(highPhysLimit)){
                bw2.write(temp);
                bw2.newLine();
                bw2.flush();
                break;
            }
        }
        while((temp=br.readLine())!=null){
            bw2.write(temp);
            bw2.newLine();
            bw2.flush();
        }
        bw.flush();
        bw.close();
        bw1.flush();
        bw1.close();
        bw2.flush();
        bw2.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
    
}
