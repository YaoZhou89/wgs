/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package wheatwgs;

import gff.modifyGTF;
import gff.splitByChr;
import math.SegeregationTest;
import script.GenerateScripts;
import vcf.VcfTools;

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
        String inFile = "",keyFileName = "",inFile2 = null;
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
        int maxdepth = 0,mindepth = 0;
        int header = 49, subnum = 1000000;
        String type = "null";
        String anchor = "",suffix=".vcf";
        int numThreads = 4, minComp = 200,windowSize = 2000,bestContrasts = -1;
        double maxIBDDist = 0.02, threshold = 0.1;
        
        String MQ = "40", FS= "10", MQRankSum = "-8", ReadPosRankSum ="-12.5",BSQRankSum="0";
        for (int i = 0; i < len; i++){
            if (null != args[i])switch (args[i]) {
                case "--model":
                    model = args[i+1];
                    i++;
                    break;
                case "--i":
                    inFile = args[i+1];
                    i++;
                    break;
                case "--o":
                    outFile = args[i+1];
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
                case "--inFile2":
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
                default:
                    break;
            }
        }
        if(outFile.equals("")){
            String[] temp = inFile.split("/");
            outFile = temp[temp.length -1].split("\\.")[0];
        }
        if(model.equals("vcf")){
            // Segeragation test
            if(type.equals("ST")){
                new SegeregationTest(inFile,outFile);
            }else if(type.equals("IBDfilter")){
                new VcfTools(anchor, inFile, outFile, minComp, maxIBDDist, windowSize, numThreads, bestContrasts);
            }else if(type.equals("LDfilter")){
                new VcfTools(inFile,outFile,windowSize,threshold);
            }else if(type.equals("merge")){
                new VcfTools(inFile,outFile,suffix);
            }
        }
        if(model.equals("gff")){
            if(type.equals("modifyGTF")){
                new modifyGTF(inFile,outFile);
            }else if(type.equals("changeCoordinate")){
                new modifyGTF(inFile,outFile,inFile2);
            }else if(type.equals("splitByChr")){
                new splitByChr(inFile);
            }
        }
        if(model.equals("GenerateScripts")){
            new GenerateScripts(inFile,outFile,type);
        }
        long endTime = System.currentTimeMillis();
        int timeLast = (int) ((endTime-startTime)/1000);
        System.out.println("Process finished in: "+ timeLast + " seconds");
    }
    
}
