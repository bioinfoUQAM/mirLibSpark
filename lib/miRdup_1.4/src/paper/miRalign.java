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
 * Submit sequences to miRalign
 */
package paper;

import java.io.*;
import java.net.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import miRdup.Features;

/**
 *
 * @author fyrox
 */
public class miRalign {
    
    public static void main(String[] args) {
        
        String infile = "all.txt"; 
        String outputdiffs="predictor"+File.separator+"all.miRalign.diffs.txt";
        String outputPredictions="predictor"+File.separator+"all.miRalign.pred.txt";
//        String infile=args[0];  
//        String outfile=args[1];
        
        ArrayList<Integer> diffStarts= new ArrayList<Integer>();        
        ArrayList<Integer> diffEnds= new ArrayList<Integer>();
        
        for (int i = 0; i < 1000; i++) {
            diffStarts.add(0);
            diffEnds.add(0);
        }
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            PrintWriter pw = new PrintWriter(new FileWriter(outputPredictions));
            String line;
            int cpt=0;
            while(br.ready()){
                try {
                    cpt++;
                    line = br.readLine();
                    //System.out.println(line);
                    String mirnaName = line.split("\t")[0];
                    String precName = line.split("\t")[1];
                    String knownMirna = line.split("\t")[2];
                    System.out.print(mirnaName+"\t");
                    String prec = line.split("\t")[3];
                    
                    //execute miRalign and get start and end position of the predicted miRNA
                    int diffStartFromMiRNA ;
                    int diffendFromMiRNA ;
                    int diffStartFromMiRNAStar = 0;
                    int diffendFromMiRNAStar = 0;
                    String predictedMiRNA;
                    String structure = null;
                    boolean mirnastar = false;
                    try {
                        Map<String, String> data = new HashMap<String, String>();
                        data.put("sequence", prec);
                        data.put("type", "1");
                        data.put("delta", "15");
                        data.put("min_seq_sim", "30");
                        data.put("MFE", "20");
                        String results = executeMiRalign(data);
                        String code = getCode(results);

                        int startPosition = code.split("\n")[5].indexOf("*");
                        int endPosition = code.split("\n")[5].lastIndexOf("*");
                        structure = code.split("\n")[7];
                        structure=structure.substring(0, structure.lastIndexOf("(")-1);
                        predictedMiRNA=code.split("\n")[6].substring(startPosition, endPosition);
                        

                        // get difference from known miRNA
                        int startOfKnownMiRNA = prec.indexOf(knownMirna);
                        int endOfKnownMiRNA = prec.indexOf(knownMirna) + knownMirna.length();

                        diffStartFromMiRNA = Math.abs(startPosition - startOfKnownMiRNA);
                        diffendFromMiRNA = Math.abs(endPosition - endOfKnownMiRNA);

                        // get difference from known miRNA star
                        Features f = new Features(knownMirna, prec, structure, true);
                        String miRNAStar=f.getMirnaStar();
                        int startOfKnownMiRNAStar=prec.indexOf(miRNAStar);
                        int endOfKnownMiRNAStar=prec.indexOf(miRNAStar)+miRNAStar.length();

                        diffStartFromMiRNAStar=Math.abs(startPosition-startOfKnownMiRNAStar);
                        diffendFromMiRNAStar=Math.abs(endPosition-endOfKnownMiRNAStar);
                        
                        if (diffStartFromMiRNAStar<diffStartFromMiRNA){
                            mirnastar=true;
                        } else {
                            mirnastar=false;
                        }
                        
                    } catch (Exception e) {
                        predictedMiRNA="null";
                        diffStartFromMiRNA=100;
                        diffendFromMiRNA=100;
                    }
                    System.out.println(diffStartFromMiRNA+"\t"+diffendFromMiRNA);
                    pw.println(mirnaName+"\t"+precName+"\t"+knownMirna+"\t"+prec+"\t"+structure+"\t"+predictedMiRNA+"\t"+diffStartFromMiRNA+"\t"+diffendFromMiRNA+"\t"+diffStartFromMiRNAStar+"\t"+diffendFromMiRNAStar);
                    
                    
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
                                       
                    if (cpt%100==0){
                        pw.flush();
                        printResults(outputdiffs, diffStarts, diffEnds);
                    }                    
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
            pw.close();
            printResults(outputdiffs, diffStarts, diffEnds);
        } catch (Exception e) {
            e.printStackTrace();
        }
                
    }
    
    public static void printResults(String outfile,ArrayList<Integer> diffStarts,ArrayList<Integer> diffEnds){
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
    
    private static String executeMiRalign2(String konwnMirna, String prec) {
        //execute miRalign and get start and end position of the predicted miRNA
        try {
            Map<String, String> data = new HashMap<String, String>();
            data.put("sequence", prec);
            data.put("type", "1");
            data.put("delta", "15");
            data.put("min_seq_sim", "30");
            data.put("MFE", "20");
            String results = executeMiRalign(data);
            String code = getCode(results);
            
            int startPosition = code.split("\n")[5].indexOf("*");
            int endPosition = code.split("\n")[5].lastIndexOf("*");

            String miRNA=code.split("\n")[6].substring(startPosition, endPosition);
            System.out.println(miRNA);
            
            // get difference from known miRNA
            int startOfKnownMiRNA = prec.indexOf(konwnMirna);
            int endOfKnownMiRNA = prec.indexOf(konwnMirna) + konwnMirna.length();
            
            int diffStart = Math.abs(startPosition - startOfKnownMiRNA);
            int diffend = Math.abs(endPosition - endOfKnownMiRNA);
            
            return diffStart + "," + diffend;
        } catch (Exception e) {
            return 100 + "," + 100;
        }
    }
    
    
    
    public static void exemple(){
        Map<String, String> data = new HashMap<String, String>();
        data.put("sequence", "ccagccugcugaagcucagagggcucugauucagaaagaucaucggauccgucugagcuuggcuggucgg");
        data.put("type", "1");
        data.put("delta", "15");
        data.put("min_seq_sim", "30");
        data.put("MFE", "20");
        String results=executeMiRalign(data);
        String code = getCode(results);
        int mirnastart=code.split("\n")[5].indexOf("*");
        int mirnaend=code.split("\n")[5].lastIndexOf("*");
        String miRNA=code.split("\n")[6].substring(mirnastart, mirnaend);
        System.out.println(miRNA);
    }
    
    /**
     * get link to results
     * @param data
     * @return 
     */
    public static String executeMiRalign(Map<String, String> data) {
        String link="";
	try{
            URL siteUrl = new URL("http://bioinfo.au.tsinghua.edu.cn/miralign/predict.php");

            HttpURLConnection conn = (HttpURLConnection) siteUrl.openConnection();
            conn.setRequestMethod("POST");
            conn.setDoOutput(true);
            conn.setDoInput(true);

            DataOutputStream out = new DataOutputStream(conn.getOutputStream());

            Set keys = data.keySet();
            Iterator keyIter = keys.iterator();
            String content = "";
            for(int i=0; keyIter.hasNext(); i++) {
                    Object key = keyIter.next();
                    if(i!=0) {
                            content += "&";
                    }
                    content += key + "=" + URLEncoder.encode(data.get(key), "UTF-8");
            }
            out.writeBytes(content);
            out.flush();
            out.close();
            BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            String line = "";
            while((line=in.readLine())!=null) {
                if (line.contains("result")){
                    link=line.substring(line.indexOf(".")+1,line.indexOf(">"));
                }                    
            }
            in.close();   
	} catch(Exception e){
             e.printStackTrace();       
        }
        return "http://bioinfo.au.tsinghua.edu.cn/miralign"+link;
    }
    
        /**
     * Get HTML code from an url
     * @param link
     * @return 
     */
    public static String getCode (String link) {
        String sourceCode = null;
        try {
            URL url = new URL(link);
            URLConnection uc = url.openConnection();
            InputStream in = uc.getInputStream();
            int c = in.read();
            StringBuilder build = new StringBuilder();
            while (c != -1) {
                build.append((char) c);
                c = in.read();
            }
            sourceCode = build.toString();
        } catch (MalformedURLException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sourceCode;
    }


}
