/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package wheatwgs;

import GO.GO;
import bed.bed;
import fasta.Fasta;
import fasta.Genome;
import fasta.MAF;
import fasta.SubstractFromFasta;
import fasta.Transform;
import fasta.contig;
import fastq.Fastq;
import format.Blast;
import format.Com;
import genome.GenomeDistribute;
import gff.gff3;
import gff.modifyGTF;
import gff.splitByChr;
import math.Correlation;
import math.Pi;
import math.SegeregationTest;
import math.gwasResult;
import math.singleRegion;
import script.GenerateScripts;
import test.test;
import text.FileSub;
import text.Summary;
import vcf.Check;
import vcf.CompareVCF;
import vcf.Frequence;
import vcf.GeneticDivergency;
import vcf.IBD;
import vcf.MergeVmapI;
import vcf.ParallelVCF;
import vcf.VcfTools;
import static vcf.VcfTools.merge;
import vcf.vcf;

/**
 *
 * @author yaozhou
 */
public class main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int len = args.length;
        String model = "test";
        String inFile = "",keyFileName = "",inFile2 = "";
        String outFile="",names = "";
        int ExSize = 100000;
        String RE1="GGATCC",RE2="CCGG";
        String bed=null;
        boolean gzip = false,hmp = false,byChr=false,het = false;
        boolean plink = false, blink = false,depth = false;
        String size = "0";
        //the regression for depth,like D group: a = -0.8, b = 0.28856
        //linear model: y = a + bx; x is depth; y is sd value.
        double a = 0, b=0,maxSD = 0; 
        // range of depth for analysis
        int maxdepth = 1000,mindepth = 0;
        int header = 49, subnum = 1000000;
        String type = "null",chr="";
        String anchor = "",suffix=".vcf";
        int numThreads = 24, minComp = 200,windowSize = 2000,bestContrasts = -1;
        double maxIBDDist = 0.02, threshold = 0.1;
        int start = 0, end =1;
        String MQ = "40", FS= "60", MQRankSum = "-12.5", ReadPosRankSum ="-8",BSQRankSum="0",SOR="3";
       
        for (int i = 0; i < len; i++){
            if (null != args[i])switch (args[i]) {
                case "--model":
                    model = args[i+1];
                    i++;
                    break;
                case "--size":
                    ExSize = Integer.parseInt(args[i+1]);   
                    size = args[i+1];
                    i++;
                    break;
                case "--cutter":
                    RE1 = args[i+1];
                    RE2 = args[i+2];
                    i = i+2;
                    break;
                case "--make-plink":
                    plink = true;
                    break;
                case "--make-blink":
                    blink = true;
                    break;
                case "--file":
                    inFile = args[i+1];
                    i++;
                    break;
                case "--out":
                    outFile = args[i+1];
                    i++;
                    break;
                case "--file2":
                    inFile2 = args[i+1];
                    i++;
                    break;
                case "--list":
                    inFile2 = args[i+1];
                    i++;
                    break;
                case "--gzip":
                    gzip = true;
                    break;
                case "--hmp":
                    //for model is GH
                    hmp = true;
                    break;
                case "--byChr":
                    //for model is GH
                    byChr = true;
                    break;
                case "--subname":
                    //for model is GH
                    names = args[i+1];
                    break;
                case "--pattern":
                    //for model is GH
                    names = args[i+1];
                    break;
                case "--key":
                    keyFileName = args[i+1];
                    i++;
                    break;  
                case "--bed":
                    bed = args[i+1];
                    i++;
                    break; 
                case "--depth":
                    depth = true;
                    break;
                case "--a":
                    a = Double.valueOf(args[i+1]);
                    i++;
                    break;
                case "--b":
                    b = Double.valueOf(args[i+1]);
                    i++;
                    break;
                case "--max-depth":
                    maxdepth = Integer.valueOf(args[i+1]);
                    i++;
                    break;
                case "--min-depth":
                    mindepth = Integer.valueOf(args[i+1]);
                    i++;
                    break;
                case "--maxSD":
                    maxSD = Double.valueOf(args[i+1]);
                    i++;
                    break;
                case "--r":
                    maxSD = Double.valueOf(args[i+1]);
                    i++;
                    break; 
                case "--het":
                    het = true;
                    break;
                case "--type":
                    type = args[i+1];
                    i++;
                    break;
                case "--subnum":
                    subnum = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--header":
                    header = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--MQ":
                    MQ = args[i+1];
                    i++;
                    break;
                case "--FS":
                    FS = args[i+1];
                    i++;
                    break;
                case "--MQRankSum":
                    MQRankSum = args[i+1];
                    i++;
                    break;
                case "--ReadPosRankSum":
                    ReadPosRankSum = args[i+1];
                    i++;
                    break;
                case "--BSQRankSum":
                    BSQRankSum = args[i+1];
                    i++;
                    break;
                case "--SOR":
                    SOR = args[i+1];
                    i++;
                    break;
                case "--anchorFile":
                    anchor = args[i+1];
                    i++;
                    break;
                case "--minComp":
                    minComp = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--maxIBDDist":
                    maxIBDDist = Double.parseDouble(args[i+1]) ;
                    i++;
                    break;
                case "--windowSize":
                    windowSize = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--bestContrasts":
                    bestContrasts = Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--threshold":
                    threshold = Double.parseDouble(args[i+1]);
                    i++;
                    break;
                case "--suffix":
                    suffix =args[i+1];
                    i++;
                    break;
                case "--chr":
                    chr =args[i+1];
                    i++;
                    break;
                case "--start":
                    start =Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--end":
                    end =Integer.parseInt(args[i+1]);
                    i++;
                    break;
                case "--TotalDepth":
                    maxSD = Double.parseDouble(args[i+1]);
                    i++;
                    break;
                case "--rate":
                    a = Double.parseDouble(args[i+1]);
                    i++;
                    break;
                
                default:
                    break;
            }
        }
        if(outFile.equals("")){
            String[] temp = inFile.split("/");
            outFile = temp[temp.length -1].split("\\.")[0];
        }
        if(model.equals("test")){
            new test(inFile,model);
        }
        if(model.equals("fasta")){
            if(type.equals("split")){
                new Genome().readByChromosome(inFile,start,end);
            }else if(type.equals("contig")){
                contig.splitContig(inFile, outFile);
            }else if(type.equals("getGene")){
                new SubstractFromFasta(inFile, outFile);
            }else if (type.equals("toFastq")){
                new Transform(inFile,outFile);
            }else if (type.equals("getChr")){
                new Genome().getChrInfo(inFile, outFile);
            }else if (type.equals("bwa")){
                new Genome().getBWAformat(inFile,outFile,a);
            }else if (type.equals("genomeSort")){
                new Genome().getSort(inFile, outFile);
            }else if(type.equals("changeName")){
                new Transform().changeName(inFile, outFile, names);
            }else if(type.equals("splitFasta")){
                new Fasta().readFasta(inFile, suffix, windowSize);
            }else if (type.equals("removeDup")){
                new Fasta().removeDup(inFile, outFile);
            }else if(type.equals("getGenes")){
                new Fasta().getGenes(inFile, inFile2, outFile);
            }else if (type.equals("splitABD")){
                new Fasta().splitABD(inFile, outFile);
            }
        }else if (model.equals("fastq")){
            new Fastq().changeName(inFile, outFile);
        }
        if(model.equals("conserved")){
            if(type.equals("region")){
                new singleRegion().getSR(inFile, outFile, maxSD, a);
            }else if (type.equals("site")){
                new singleRegion().getSS(inFile, outFile, a, b);
            }else if (type.equals("getBed")){
                new singleRegion().getBed(inFile, outFile);
            }else if (type.equals("getBedPi")){
                new Pi().getBedPi(inFile, bed, outFile);
            }else if (type.equals("toChr")){
                new singleRegion().toChr(inFile, outFile);
            }else if(type.equals("getModified")){
                new singleRegion().getLen(inFile, outFile);
            }else if (type.equals("getRegionSum")){
                new singleRegion().getBedTotal(inFile, outFile,names);
            }
            
        }
        if(model.equals("bed")){
            if(type.equals("sortBed")){
                new bed().sort(inFile,outFile);
            }else if(type.equals("getRegion")){
                new bed().getRegion(inFile, outFile);
            }else if (type.equals("splitBychr")){
                new bed().splitByChr(inFile, outFile);
            }else if (type.equals("mummer2bed")){
                new bed().mummer2bed(inFile, outFile, names, windowSize);
            }else if (type.equals("repeat2bed")){
                new bed().repeat2bed(inFile, outFile, windowSize);
            }else if (type.equals("gene2bed")){
                new bed().gene2bed(inFile, outFile, windowSize);
            }
           
        }
        if(model.equals("genome")){
            if(type.equals("byChr")){
                new GenomeDistribute(inFile,outFile,true);
            }
        }
        if(model.equals("VmapI")){
            if(type.equals("mergeBarley")){
                new MergeVmapI().mergeBarley(inFile, inFile2, outFile);
            }else if (type.equals("getS")){
                new FileSub().getS(inFile, outFile);
            }else if(type.equals("getSet")){
                new FileSub().getSet(inFile, outFile);
            }
            
        }
        if(model.equals("MAF")){
            if(type.equals("getSeq")){
                new MAF().getSeq(inFile,outFile);
            }
        }
        if(model.equals("GO")){
            new GO().toGO(inFile, keyFileName, outFile);
        }
        if(model.equals("GWAS")){
            if(type.equals("plink")){
                new gwasResult().thinResult(inFile, outFile, maxSD, threshold);
            }else if (type.equals("emmax")){
                new gwasResult().thinEMMAX(inFile, outFile, maxSD, threshold);
            }
        }
        if(model.equals("vcf")){
            // Segeragation test
            if(type.equals("ST")){
                new SegeregationTest(inFile,outFile);
            }else if(type.equals("IBDfilter")){
//                new IBD().ibdFilterparrallel(anchor, inFile, outFile, minComp, maxIBDDist, windowSize);
//                new IBD().ibdFilter(anchor, inFile, outFile, minComp, maxIBDDist, windowSize);
                new VcfTools(anchor, inFile, outFile, minComp, maxIBDDist, windowSize, numThreads, bestContrasts);
            }else if(type.equals("LDfilter")){
                new VcfTools(inFile,outFile,windowSize,threshold);
            }else if(type.equals("merge")){
//                VcfTools.mergeVCF(inFile,outFile);
                VcfTools.merge(inFile,inFile2,outFile);
            }else if (type.equals("quality")){
                new VcfTools(inFile,outFile,MQ,FS,MQRankSum, ReadPosRankSum,BSQRankSum,SOR);
            }else if(type.equals("compare")){
                new CompareVCF(inFile,inFile2,suffix);
            }else if (type.equals("all")){
                new ParallelVCF( inFile,  suffix, type, outFile,
             MQ, FS, MQRankSum, ReadPosRankSum, BSQRankSum,
             SOR, maxdepth,  mindepth,  maxSD,  windowSize,
             threshold);
            }else if (type.equals("check")){
                Check.checkLastLine(inFile, suffix);
            }else if(type.equals("addST")){
                Check.addST(inFile,outFile,MQ);
            }else if(type.equals("deleteST")){
                Check.deleteST(inFile,outFile,MQ);
            }else if(type.equals("subset")){
                VcfTools.getSub(inFile, outFile, ExSize, subnum);
            }else if(type.equals("DepthFilter")){
                new VcfTools(inFile,outFile,a,b,maxSD,mindepth,maxdepth);
            }else if(type.equals("rsd")){
                new VcfTools(inFile,outFile,maxSD,suffix);
            }else if (type.equals("GeneticDivergency")){
                GeneticDivergency.calPair(inFile, outFile);
            }else if(type.equals("splitByChr")){
                new VcfTools(inFile,ExSize);
//                VcfTools.splitByChr(inFile, outFile, new vcf());
            }else if(type.equals("chrs")){
                VcfTools.subChr(inFile, chr, outFile);
            }else if (type.equals("winDepth")){
                new VcfTools(inFile,outFile,windowSize);
            }else if (type.equals("splitByChrs")){
                VcfTools.splitByChrs(inFile, outFile);
            }else if (type.equals("check2vcf")){
                new CompareVCF().compare(inFile, outFile);
            }else if (type.equals("twoStatic")){
                new Check().getTwo(inFile);
            }else if (type.equals("mergeVCF")){
                VcfTools.mergeVCF(inFile,outFile);

            }else if (type.equals("getDerivedSNP")){
                new VcfTools().getDerivedSNP(inFile, inFile2, outFile);
            }else if(type.equals("checkAll")){
                new Check().checkAll(inFile, outFile);
            }else if (type.equals("changeHeader")){
                new VcfTools().changeHeader(inFile, outFile, suffix);
            }else if(type.equals("mergeAll")){
                new VcfTools().mergeAll(inFile, outFile, suffix);
            }else if (type.equals("mergeFiles")){
                new VcfTools().mergeFiles(inFile, inFile2, outFile);
            }
            else{
                new VcfTools(inFile,outFile,type,maxSD);
            }
        }
        if(model.equals("txt")){
            if(type.equals("summary")){
                new Summary().readData(inFile,suffix);
            }
        }
        if(model.equals("gff")){
            if(type.equals("modifyGTF")){
                new modifyGTF(inFile,outFile);
            }else if(type.equals("changeCoordinate")){
                new modifyGTF(inFile,outFile,inFile2);
            }else if(type.equals("splitByChr")){
                new splitByChr(inFile);
            }else if (type.equals("toBed")){
                new gff3().toBed(inFile, inFile2, outFile, names);
            }
        }
        if(model.equals("frequence")){
            Frequence.sumFre(inFile, suffix, outFile);
        }
        if(model.equals("GenerateScripts")){
            if(type.equals("mergeLineage")){
                GenerateScripts.mergeLineage(inFile, outFile, chr);
            }else if(type.equals("splitBamByChr")){
                new GenerateScripts(inFile,outFile,chr);
            }else if (type.equals("samtoolsCalling")){
                new GenerateScripts().getSamtoolsCalling(inFile, outFile,suffix,chr);
            }
            else {
                new GenerateScripts(inFile,outFile,type,suffix);
            }
            
        }
        if(model.equals("blast")){
            if(type.equals("changeBlastp")){
                new Blast().changeBlastp(inFile,outFile);
            }
        }
        if(model.equals("PI")){
            if(type.equals("sitePi")){
                new Pi().getSitePi(inFile, inFile2, outFile);
            }else if (type.equals("getBed")){
                new Pi().forPlot(inFile, outFile);
            }
        }
        if(model.equals("correlation")){
            if (type.equals("siteCor")){
                new Correlation().siteCor(inFile, outFile, size);
            }
        }
        if(model.equals("compare")){
            if(type.equals("getSame")){
                new Com().getSame(inFile, inFile2, outFile);
            }
        }
        long endTime = System.currentTimeMillis();
        int timeLast = (int) ((endTime-startTime)/1000);
        System.out.println("Process finished in: "+ timeLast + " seconds");
    }
    
}
