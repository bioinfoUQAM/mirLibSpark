/*
 * Read promir Results
 */
package paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import miRdup.Features;
import static paper.MatureBayes.execute;

/**
 *
 * @author fyrox
 */
public class ProMir {
    
    public static void main(String[] args) {
        File infile = new File("tools"+File.separator+"ProMiR"+File.separator+"fastas"+File.separator+"ProMir.outputs");
        File outfile = new File("tools"+File.separator+"ProMiR"+File.separator+"fastas"+File.separator+"ProMir.diffs.txt");
        execute(infile, outfile);
        
        
//        String tab[]= {"all","Arthropoda","mammal","Nematoda","Pisces","Viridiplantae"};
//        
//        for (String s : tab) {
//            File infile = new File("MatureBayes."+s+".txt");
//            File outfile = new File("MatureBayes."+s+".diffs.txt");
//            execute(infile, outfile);
//        }
        
        
    }
    
    public static void execute(File infile, File outfile) {
        
        
        ArrayList<Integer> diffStarts= new ArrayList<Integer>();        
        ArrayList<Integer> diffEnds= new ArrayList<Integer>();
        
        for (int i = 0; i < 1000; i++) {
            diffStarts.add(0);
            diffEnds.add(0);
        }
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            String line="";
            while (br.ready()){    
                line=br.readLine();
                try {                    
                    String knownMirna = line.split("_")[2];                    
                    String prec = line.split("_")[3];
                    String structure = line.split("_")[4];
                    
                    int startOfKnownMiRNA = prec.indexOf(knownMirna);
                    int endOfKnownMiRNA = startOfKnownMiRNA + knownMirna.length();
                    line = br.readLine();
                    int startPosition = Integer.valueOf(line.substring(line.indexOf("[")+1,line.indexOf(",")));
                    int endPosition= Integer.valueOf(line.substring(line.indexOf(",")+1,line.indexOf("]")));                        
                    line=br.readLine();
                    line=br.readLine();
                    line=br.readLine();
                    line=br.readLine();
                    
                    // get difference from known miRNA                    
                    int diffStartFromMiRNA = Math.abs(startPosition - startOfKnownMiRNA);
                    int diffendFromMiRNA = Math.abs(endPosition - endOfKnownMiRNA);
                    
                    // get difference from known miRNA star
                    Features f = new Features(knownMirna, prec, structure, true);
                    String miRNAStar=f.getMirnaStar();
                    int startOfKnownMiRNAStar=prec.indexOf(miRNAStar);
                    int endOfKnownMiRNAStar=prec.indexOf(miRNAStar)+miRNAStar.length();

                    int diffStartFromMiRNAStar=Math.abs(startPosition-startOfKnownMiRNAStar);
                    int diffendFromMiRNAStar=Math.abs(endPosition-endOfKnownMiRNAStar);
                    boolean mirnastar=false;
                    
                    if (diffStartFromMiRNAStar<diffStartFromMiRNA){
                        mirnastar=true;
                    } else {
                        mirnastar=false;
                    }

                    //differences from known mirna
                    if (mirnastar){
                        int tmp = diffStarts.get(diffStartFromMiRNAStar);
                        tmp += 1;
                        diffStarts.set(diffStartFromMiRNAStar, tmp);

                        tmp = diffEnds.get(diffendFromMiRNAStar);
                        tmp += 1;
                        diffEnds.set(diffendFromMiRNAStar, tmp);
                        
                        System.out.println("\t"+diffStartFromMiRNAStar+"\t"+diffendFromMiRNAStar);
                    } else {
                        int tmp = diffStarts.get(diffStartFromMiRNA);
                        tmp += 1;
                        diffStarts.set(diffStartFromMiRNA, tmp);

                        tmp = diffEnds.get(diffendFromMiRNA);
                        tmp += 1;
                        diffEnds.set(diffendFromMiRNA, tmp); 
                    }
                    
                } catch (Exception e) {
                    e.printStackTrace();
                } 

            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        // print results
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(outfile));
            pw.println("Start");
            int flagStart = 0;
            for (int i = diffStarts.size() - 1; i >= 0; i--) {
                if (diffStarts.get(i) != 0) {
                    flagStart = i;
                    break;
                }                
            }
            for (int i = 0; i <= flagStart; i++) {
                pw.println(i + "\t" + diffStarts.get(i));                
            }
            
            pw.println("End");
            int flagEnd = 0;
            for (int i = diffEnds.size() - 1; i >= 0; i--) {
                if (diffEnds.get(i) != 0) {
                    flagEnd = i;
                    break;
                }                
            }
            for (int i = 0; i <= flagEnd; i++) {
                pw.println(i + "\t" + diffEnds.get(i));                
            }
            
            pw.flush();pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
