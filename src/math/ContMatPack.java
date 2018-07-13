/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package math;

import net.sf.javaml.utils.*;
import jsc.combinatorics.*;
import jsc.contingencytables.ChiSquaredTest;
import jsc.contingencytables.ContingencyTable;
import jsc.distributions.*;

/**
 *
 * @author bukowski
 * 
 * Computes the p-value for a 2 byN contingency matrix.
 * If the expected values are big enough, the chisq test is done,
 * otherwise the p-value is obtained from simulation of the distribution
 * of chisq statistics (with fixed marginals).
 */
public class ContMatPack {
    
    private double[][] matrix;
    private int [] rowtotals;
    private int [] coltotals;
    private int nrows;
    private int ncols;
    private int nall; 
    private double chisq;
    private double pvalue;
    private boolean gotpvalue;
    
    public ContMatPack(double[][] data)
    {
        nrows = data.length;
        ncols = data[0].length;
        rowtotals = new int[nrows];
        coltotals = new int[ncols];
 
        // Get row totals and number of all counts
        nall = 0;
        for(int i=0;i<nrows;i++)
        {
            rowtotals[i] = 0;
            for(int j=0;j<ncols;j++)
            {
                rowtotals[i] += data[i][j];
            }
            nall += rowtotals[i];
        }
        // Make column totals
        for(int j=0;j<ncols;j++)
        {
            coltotals[j] = 0;
            for(int i=0;i<nrows;i++)
            {
                coltotals[j] += data[i][j];
            }
        }
        
        matrix = data;
        chisq = ContingencyTables.chiVal(matrix, false);
        
        /*
        ContingencyTable ctbl = new ContingencyTable(data);
        ChiSquaredTest chsq = new ChiSquaredTest(ctbl);
        
        
        if(! chsq.hasSmallExpectedFrequency())
        {
            pvalue = chsq.getSP();
            gotpvalue = true;
        }
        else
        {
            gotpvalue = false;
            // make double copy of matrix and save marginals
            matrix = new double[nrows][ncols];
            for(int i=0;i<nrows;i++)
            {
                for(int j=0;j<ncols;j++)
                {
                    matrix[i][j] = data[i][j];
                }
            }
            rowtotals = ctbl.getRowTotals();
            coltotals = ctbl.getColumnTotals();
            
            // Compute the original chisq statisics
            chisq = ContingencyTables.chiVal(matrix, false);
        }
        */
    }
    
    public double get_pvalue()
    {
        if(gotpvalue)
        {
            return pvalue;
        }
        else
        {
            // Make sute matrix is initialized!!!
            pvalue = SimulatePvalueG(matrix,1000);
            return pvalue;
        }
    }
    
    public static double pvalue(double[][] mat)
    {
        double pv = 0;
        if(ContingencyTables.cochransCriterion(mat))
        {
            pv = ContingencyTables.chiSquared(mat, false);
        }
        else
        {
            System.out.println("Swithing to simulation");
            pv = SimulatePvalueHG(mat,1000);
        }
        return pv;
    }
    
    public double ValG(double [][] data)
    {
        double Gres = 0;
        double expected=0;
        for(int i=0;i<nrows;i++)
        {
            for(int j=0;j<ncols;j++)
            {
                if(data[i][j] > 0)
                {
                    expected = (double)rowtotals[i]*coltotals[j]/nall;
                    Gres += data[i][j]*Math.log(data[i][j]/expected);
                }
            }
        }
        
        return 2*Gres;
    }
    
    public double pvalueG0(double [][] data)
    {
        double Gres = ValG(data);
        double df = (nrows-1)*(ncols-1);
        if(Gres < 0)
        {
            System.out.println("Negative Gres: "+Gres);
        }
        return ChiSquared.upperTailProb(Gres, df);
    }
    
    public double pvalueG(double [][] data)
    {
        double Gres = ValG(data);
        double df = (nrows-1)*(ncols-1);
        if(Gres < 0)
        {
            System.out.println("Negative Gres: "+Gres);
        }
        ChiSquared csq = new ChiSquared(df);
        return 1 - csq.cdf(Gres);
    }
    
    //
    // Calculate p-value from the distribution of the G statistic
    //
    public double SimulatePvalueG(double[][] mat, int numiters)
    {
        double pvalue=1;
        double mychisq;
        int [] newcolorder = new int[ncols];
        Hypergeometric hgm;
        
        int nwhite = rowtotals[0];
        int nblack = rowtotals[1];
        int nwhite_new;
        int larger_random = 0;
        
        double chisq_orig = ValG(mat);
        
        Permutations perms = new Permutations(ncols);
        for(int iter=0;iter<numiters;iter++)
        {
            newcolorder = perms.randomPermutation().toIntArray();
            nwhite = rowtotals[0];
            nblack = rowtotals[1];
            // Simulate each column except last, one by one
            for(int icol1=0;icol1<ncols-1;icol1++)
            {
                int icol = newcolorder[icol1]-1;
                
                if(nblack == 0)
                {
                    nwhite_new = coltotals[icol];
                }
                else if(nwhite==0)
                {
                    nwhite_new = 0;
                }
                else
                {
                    hgm = new Hypergeometric(coltotals[icol],nwhite+nblack,nwhite);
                    nwhite_new = (int)hgm.random();
                }
                
                //System.out.println("icol "+icol);
                mat[0][icol] = nwhite_new;
                mat[1][icol] = coltotals[icol] -mat[0][icol];
                nwhite -=mat[0][icol];
                nblack -= mat[1][icol];
            }
            // fill the last column with whatever remained:
            mat[0][newcolorder[ncols-1]-1] = nwhite;
            mat[1][newcolorder[ncols-1]-1] = nblack;
            
            // verify column total:
            if(nwhite+nblack != coltotals[newcolorder[ncols-1]-1])
            {
                System.out.println("Last column wrong!");
            }
            
            // Copute the chisq statistics
            mychisq = ValG(mat);
            //System.out.println("chisq: " + mychisq);
            if(mychisq >= chisq_orig) { larger_random++; }
        }
        pvalue = (double)larger_random/numiters;
        
        return pvalue;
    }
    
    //
    // Calculate p-value from the distribution of the G statistic
    //
    public double SimulatePvalueAsa159(int numiters, int seed)
    {
        double pvalue=1;
        double mychisq;
        int larger_random = 0;
        double chisq_orig = chisq;
        //int seed = 123456789; // should chage it somehow
        
        asa159 asaobj = new asa159(nrows,ncols,rowtotals,coltotals,seed);
        
        if(asaobj.ierror > 0) { return 100+asaobj.ierror; } // Various singular conditions, like
                                                            // Only one taxon (102), no ref allele (103)
        
        for(int iter=0;iter<numiters;iter++)
        {
            matrix = asaobj.rcont2();
            // Compute the chisq statistics
            mychisq = ContingencyTables.chiVal(matrix, false);
            //System.out.println("chisq: " + mychisq);
            if(mychisq >= chisq_orig) { larger_random++; }
        }
        pvalue = (double)larger_random/numiters;
        
        return pvalue;
    }
    //
    // Calculate the p-value from the distribution of the chi-square statistic
    //
    public static double SimulatePvalueHG(double[][] mat, int numiters)
    {
        double pvalue=1;
        double mychisq;
        int tst_number;
        int ncolumns = mat[0].length;
        int [] rtotals = new int[2];
        int [] ctotals = new int[ncolumns];
        int [] newcolorder;
        Hypergeometric hgm;
        
        // Get row and column totals
        for(int i=0;i<2;i++)
        {
            rtotals[i] = 0;
            for(int j=0;j<ncolumns;j++)
            {
                rtotals[i] += mat[i][j];
            }
        }
        // Make column totals
        for(int j=0;j<ncolumns;j++)
        {
            ctotals[j] = 0;
            for(int i=0;i<2;i++)
            {
                ctotals[j] += mat[i][j];
            }
        }
        int nwhite = rtotals[0];
        int nblack = rtotals[1];
        int nwhite_new;
        int larger_random = 0;
        
        double chisq_orig = ContingencyTables.chiVal(mat, false);
        
        Permutations perms = new Permutations(ncolumns);
        for(int iter=0;iter<numiters;iter++)
        {
            newcolorder = perms.randomPermutation().toIntArray();
            nwhite = rtotals[0];
            nblack = rtotals[1];
            // Simulate each column except last, one by one
            for(int icol1=0;icol1<ncolumns-1;icol1++)
            {
                int icol = newcolorder[icol1]-1;
                
                if(nblack == 0)
                {
                    nwhite_new = ctotals[icol];
                }
                else if(nwhite==0)
                {
                    nwhite_new = 0;
                }
                else
                {
                    hgm = new Hypergeometric(ctotals[icol],nwhite+nblack,nwhite);
                    nwhite_new = (int)hgm.random();
                }
                
                //System.out.println("icol "+icol);
                mat[0][icol] = nwhite_new;
                mat[1][icol] = ctotals[icol] -mat[0][icol];
                nwhite -=mat[0][icol];
                nblack -= mat[1][icol];
            }
            // fill the last column with whatever remained:
            mat[0][newcolorder[ncolumns-1]-1] = nwhite;
            mat[1][newcolorder[ncolumns-1]-1] = nblack;
            
            // verify column total:
            if(nwhite+nblack != ctotals[newcolorder[ncolumns-1]-1])
            {
                System.out.println("Last column wrong!");
            }
            
            //System.out.println("Matrix in iteration " + iter);
            //System.out.println(matrix[0][0] + " " + matrix[0][1]);
            //System.out.println(matrix[1][0] + " " + matrix[1double][1]);
            
            // Copute the chisq statistics
            mychisq = ContingencyTables.chiVal(mat, false);
            //System.out.println("chisq: " + mychisq + " " + chisq_orig);
            if(mychisq >= chisq_orig) { larger_random++; }
        }
        pvalue = (double)larger_random/numiters;
        
        return pvalue;
    }
    
    public double hypergeomLratio(double[][] data)
    {
        double Lratio = 0;
        
        // Construct matrix of expected frequencies
        double[][] expmat = new double[nrows][ncols];
        for(int i=0;i< nrows;i++)
        {
            for(int j=0;j<ncols;j++)
            {
                expmat[i][j] = rowtotals[i]*coltotals[j]/nall;
            }
        }
        
        double prob0 = ContingencyTables.log2MultipleHypergeometric(data);
        double prob1 = ContingencyTables.log2MultipleHypergeometric(expmat);
        Lratio = prob0 - prob1;
        Lratio = Math.pow(2, Lratio);
        //System.out.println("prob0, prob1: "+prob0 + " " + prob1);
        
        return Lratio;
    }
    
    public void close()
    {
        matrix = null;
        System.gc();
    }
    
}
