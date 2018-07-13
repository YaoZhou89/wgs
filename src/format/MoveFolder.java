/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package format;

import io.IOUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.TreeSet;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class MoveFolder {
    public MoveFolder(String Folder,String out){
        this.getMoved(Folder,out);
    }
    public MoveFolder(String Folder,String out,String type){
        if(type.equals("delete")){
            this.getDeleted(Folder);
        }
        if(type.equals("stat")){
            this.getStat(Folder, out);
        }
    }
    
    public void getMoved(String path,String outpath){ 
        BufferedWriter out = IOUtils.getTextWriter(path+"/move.sh");
        File test = new File (path);
        String SourcePath = null;
        String DesPath= null;
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, "q.gz");
        File[] reFs = getPath(subFs);   
        try {
        for (File entry : reFs){
            File[] filename = IOUtils.listRecursiveFiles(entry);
            if(filename.length>2){
                for(File name :filename){
                    String[] len = name.getAbsoluteFile().toString().split("/");
                    File des = new File(outpath+"/"+len[len.length-2]);
                    if(!des.exists()) des.mkdir();
                    if(len[len.length-1].startsWith("merg")){
                        out.write("nohup rsync -avzP\t"+name+"\t" +des+"/"+len[len.length-1] + "\t");
                        out.newLine();
                    }
                   
                }
            }else{
                for(File name :filename){
                    String[] len = name.getAbsoluteFile().toString().split("/");
                    File des = new File(outpath+"/"+len[len.length-2]);
                    if(!des.exists()) des.mkdir();
                    File desFile = new File(des+"/"+len[len.length-1]);
                    if(!desFile.exists()){
                        out.write("nohup rsync -avzP\t"+name+"\t" +des+"/"+len[len.length-1] + "\t");
                        out.newLine();
                    }
                    
                }
            }
        }
        out.flush();
        out.close();
        
        } catch (IOException ex) {
            Logger.getLogger(MoveFolder.class.getName()).log(Level.SEVERE, null, ex);
        }
        
        
    }
    private void getDeleted(String path){
        File test = new File (path);
//        File[] fs = IOUtils.listRecursiveFiles(test);
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
//        File[] subFs = IOUtils.listFilesEndsWith(fs, ".sh");
//        File[] reFs = getPath(fs);   
        
        for (File entry : test.listFiles()){
                entry.delete();
//            if(filename.length < 2) filename
        }
    }
    
    public File[] getPath(File[] folder){
        TreeSet<File> fileTree = new TreeSet();
        TreeSet<File> subfile = new TreeSet();
        for (File entry : folder){
            if(entry.isHidden()) entry.delete();
            File path = new File(entry.getParent());
            fileTree.add(path);
        }
        return fileTree.toArray(new File[fileTree.size()]);
    }
    
    private void getStat(String path,String outFile){
        try {
            File test = new File (path);
            File[] fs = IOUtils.listRecursiveFiles(test);
            File[] subFs = IOUtils.listFilesEndsWith(fs, ".gz");
            BufferedWriter bw = IOUtils.getTextWriter(outFile);
            for (File entry : subFs){
                    long size = entry.length()/1024/1024;
                    if(size<10) continue;
                    bw.write(entry.toString());
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(MoveFolder.class.getName()).log(Level.SEVERE, null, ex);
        }
        
    }
}
