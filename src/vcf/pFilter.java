/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package vcf;

import net.maizegenetics.dna.snp.GenotypeTable;
import static net.maizegenetics.dna.snp.GenotypeTable.*;
import static net.maizegenetics.dna.snp.GenotypeTableUtils.*;
import net.maizegenetics.taxa.distance.DistanceMatrix;
//import net.maizegenetics.popgen.distance.DistanceMatrix;
//import net.maizegenetics.popgen.distance.IBSDistanceMatrix;

/**
 *
 * @author bukowski
 */
public class pFilter implements Runnable {
    
    private GenotypeTable b; 
    private DistanceMatrix dm;
    private int i; // i anchor taxon of the pair
    private double maxIBDDist;
    private int [] aTobTaxa;
    private int numTestSites;
    private int [][] rightNuc;
    private int [][] wrongNuc;
    private byte [][] diploids;
    private int [] het;
    private int [] minorComp;
    private int [] minorSame;
    private int [] mjSame;
    private int [] diff;
    private int [] counts;
    private double [] mean;
    private String [] contrastPairs;
    private String [] failedContrasts;
    
    public pFilter(GenotypeTable bIn, DistanceMatrix dmIn, int iIn, double maxIBDDistIn, int [] aTobTaxaIn, int numTestSitesIn,
            int [][] rightNucIn, int [][] wrongNucIn, byte [][] diploidsIn, int [] hetIn,
            int [] minorCompIn, int [] minorSameIn, int [] mjSameIn, int [] diffIn, int [] countsIn, double [] meanIn,
            String [] contrastPairsIn, String [] failedContrastsIn)
    {
        b = bIn;
        dm = dmIn;
        i = iIn;
        maxIBDDist = maxIBDDistIn;
        aTobTaxa = aTobTaxaIn;
        numTestSites = numTestSitesIn;
        rightNuc = rightNucIn;
        wrongNuc = wrongNucIn;
        diploids = diploidsIn;
        het = hetIn;
        minorComp = minorCompIn;
        minorSame = minorSameIn;
        mjSame = mjSameIn;
        diff = diffIn;
        counts = countsIn;
        mean = meanIn;
        contrastPairs = contrastPairsIn;
        failedContrasts = failedContrastsIn;
    }
    
    @Override
    public void run()
    {
        run_pFilter();
    }
    
    private void run_pFilter()
    {
        ///Thread.currentThread().
        int ibd=0;
        int count = 0;
        int taxaIBD = 0;
        int countIBD = 0;
                for (int j=0; j<dm.numberOfTaxa(); j++) {
                    if(i==j || Double.isNaN(dm.getDistance(i,j))) continue;
                    mean[0]+=dm.getDistance(i,j);
                    count++;
                    if(dm.getDistance(i,j)<=maxIBDDist) {
                        // The following does not seem to make sense.....
                        //double[] distB=IBSDistanceMatrix.computeHetBitDistances(b,i,j);
                        //ibd++;
                        contrastPairs[0] += dm.getColumnName(i) + ":" + dm.getColumnName(j) + "\n";
                        int bI=aTobTaxa[i];
                        int bJ=aTobTaxa[j];
                        if(bI<0 || bJ<0) continue;
                        ibd++;
                        //double[] distB=IBSDistanceMatrix.computeHetBitDistances(b,bI,bJ);
                        for (int k=0; k<numTestSites; k++) {
                            byte t1b=b.genotype(bI,k);
                            byte t2b=b.genotype(bJ,k);
                            //b.genotypeAsString(bI, k)
                            if(t1b==UNKNOWN_DIPLOID_ALLELE) continue;
                            if(t2b==UNKNOWN_DIPLOID_ALLELE) continue;
                            byte[] t1ba=getDiploidValues(t1b);
                            byte[] t2ba=getDiploidValues(t2b);
                            if(t1ba[0]==t2ba[0] || t1ba[0]==t2ba[1]) {rightNuc[t1ba[0]][k]++;} else {wrongNuc[t1ba[0]][k]++;}
                            if(t1ba[1]==t2ba[0] || t1ba[1]==t2ba[1]) {rightNuc[t1ba[1]][k]++;} else {wrongNuc[t1ba[1]][k]++;}
                            // Uncomment two next lines to ignire hets
                            //if(t1b!=diploids[0][k] && t1b!=diploids[2][k]) continue;
                            //if(t2b!=diploids[0][k] && t2b!=diploids[2][k]) continue;
                            if(isHeterozygous(t1b)  || isHeterozygous(t2b)) het[k]++;
                            if(t1b==diploids[2][k] || t2b==diploids[2][k])  {
                                minorComp[k]++;
                                if(t1b==t2b) minorSame[k]++;
                            }
                            if(t1b==t2b) 
                            {
                                mjSame[k]++;
                            } 
                            else 
                            {
                                diff[k]++;
                                failedContrasts[k] += dm.getColumnName(i) + ":" + dm.getColumnName(j) + ":" + b.genotypeAsString(bI, k) + ":" + b.genotypeAsString(bJ, k) + "\n";
                            }

                        }
                        //System.out.printf("%d %d %g %g %n",i,j,dm.getDistance(i,j),distB[0]);
                    }
                }
                if(ibd>0) taxaIBD++;
                
                // int [][] acounts = new int[numTaxa][3]; // countIBD, taxaIBD, count
                counts[0] += ibd;
                counts[1] += taxaIBD;
                counts[2] += count;
    }
    
}
