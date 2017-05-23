/*
 *  miRdup v1.0
 *  Computational prediction of the localization of microRNAs within their pre-miRNA
 *  
 *  Copyright (C) 2013  Mickael Leclercq
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Get differences from predicted file of MatureBayes
 */
package paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;

/**
 *
 * @author fyrox
 */
public class MatureBayes {
    public static void main(String[] args) {
        String tab[]= {"all","Arthropoda","mammal","Nematoda","Pisces","Viridiplantae"};
        
        for (String s : tab) {
            File infile = new File("MatureBayes."+s+".txt");
            File outfile = new File("MatureBayes."+s+".diffs.txt");
            execute(infile, outfile);
        }
        
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
                if (line.startsWith(">")){
                    try {
                        String name = line.split(" ")[0];
                        String miRNA = line.split(" ")[1].split("\t")[0];
                        line = br.readLine();
                        String precursor = line;
                        
                        int startOfKnownMiRNA = precursor.indexOf(miRNA);
                        int endOfKnownMiRNA = startOfKnownMiRNA + miRNA.length();
                        line = br.readLine();
                        int startPosition5p = Integer.valueOf(line.split("\t")[0]);
                        int endPosition5p = startPosition5p + 22;                        
                        int diffStart5p = Math.abs(startPosition5p - startOfKnownMiRNA);
                        int diffEnd5p = Math.abs(endPosition5p - endOfKnownMiRNA);
                        
                        
                        line = br.readLine();
                        int startPosition3p = Integer.valueOf(line.split("\t")[0]);                        
                        int endPosition3p = startPosition3p + 22;
                        int diffStart3p = Math.abs(startPosition3p - startOfKnownMiRNA);
                        int diffEnd3p = Math.abs(endPosition3p - endOfKnownMiRNA);
                        
                        int diffStart;
                        int diffEnd;
                        boolean fivep = true;                        
                        if (diffStart5p < diffStart3p) {
                            diffStart = diffStart5p;
                            diffEnd = diffEnd5p;
                        } else {
                            fivep = false;
                            diffStart = diffStart3p;
                            diffEnd = diffEnd3p;
                        }

//                    System.out.println(name);
//                    System.out.println("\t diffStart5p="+diffStart5p);
//                    System.out.println("\t diffEnd5p="+diffEnd5p);
//                    System.out.println("\t diffStart3p="+diffStart3p);
//                    System.out.println("\t diffEnd3p="+diffEnd3p);
//                    if (fivep){
//                        System.out.println("\t Chosen=5p");
//                    } else {
//                        System.out.println("\t Chosen=3p");
//                    }

                        //differences from known mirna                    
                        int tmp = diffStarts.get(diffStart);
                        tmp += 1;
                        diffStarts.set(diffStart, tmp);
                        
                        tmp = diffEnds.get(diffEnd);
                        tmp += 1;
                        diffEnds.set(diffEnd, tmp);
                    } catch (IOException iOException) {
                    } catch (NumberFormatException numberFormatException) {
                    }
                } else {
                    line=br.readLine();
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
