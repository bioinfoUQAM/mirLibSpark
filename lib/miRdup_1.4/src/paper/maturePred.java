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
 * Submit sequences to maturepred
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
public class maturePred {
    public static boolean randomlines=true;
        
    public static void main(String[] args) {
        int num=10;
        String folder="predictor"+File.separator+"maturepred"+File.separator;
        String organisms="organisms.txt";
        String infile = folder+"all"+num+".txt"; 
        System.out.println(infile);
        String outputdiffs=folder+"all.maturepred.diffs."+num+".txt";
        String outputPredictions=folder+"all.maturepred.pred."+num+".txt";
//        String infile=args[0];  
//        String outfile=args[1];
        
        //Add organisms to HashMap
        HashMap<String,String> hmorg= new HashMap<String, String>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(organisms));
            String line="";
            while(br.ready()){
                line=br.readLine();
                hmorg.put(line.split("\t")[0], line);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        //maturepred
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
                    
                    if (randomlines){
                        int randomNum = 1 + (int)(Math.random()*100);
                        for (int i = 0; i < randomNum; i++) {
                            line = br.readLine();                            
                        }
                    }
                    
                    //System.out.println(line);
                    String mirnaName = line.split("\t")[0];
                    String org=mirnaName.split("-")[0];
                    String precName = line.split("\t")[1];
                    String knownMirna = line.split("\t")[2];
                    System.out.print(mirnaName+"\t");
                    String prec = line.split("\t")[3];
                    String structure = line.split("\t")[4];
                    //execute maturepred and get start and end position of the predicted miRNA
                    int diffStartFromMiRNA ;
                    int diffendFromMiRNA ;
                    int diffStartFromMiRNAStar = 0;
                    int diffendFromMiRNAStar = 0;
                    String predictedMiRNA;
                    boolean mirnastar = false;
                    try {
                        Map<String, String> data = new HashMap<String, String>();
                        boolean plant=false;
                        int size=22;
                        if (hmorg.get(org).contains("Viridiplantae")){
                            plant=true;
                            size=21;
                        }
                        if (plant){
                            data.put("trainmodel", "plant");
                            data.put("len", "21");
                        } else {
                            data.put("trainmodel", "animal");
                            data.put("len", "22");
                        }                        
                        data.put("ex", "blank");
                        data.put("seq", ">seq\n"+prec);
                        int startPosition=0;
                        int count=0;
                        do {                            
                            count++;
                            startPosition = submit(data)-1;
                            if (startPosition==999){
                                System.err.println("err");
                            }
                                
                        } while (startPosition==999&&count<=5);
                        
                        int endPosition = startPosition+size;
                        predictedMiRNA=prec.substring(startPosition, endPosition);
                        

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
                        e.printStackTrace();
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
    
    /**
     * get link to results
     * @param data
     * @return 
     */
    public static int submit(Map<String, String> data) {
        String mirna="";
        int mirnaStart = 1000;
	try{
            URL siteUrl = new URL("http://nclab.hit.edu.cn/maturepred/predict.php");

            HttpURLConnection conn = (HttpURLConnection) siteUrl.openConnection();
            conn.setRequestMethod("POST");
            conn.setDoOutput(true);
            conn.setDoInput(true);
            conn.setRequestProperty("Accept-Charset", "UTF-8");
            conn.setRequestProperty("Content-Type", "application/x-www-form-urlencoded;charset=UTF-8");
            
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
            content=content.replace("%0", "%0D%0");
            out.writeBytes(content);
            out.flush();
            out.close();
            BufferedReader in = new BufferedReader(new InputStreamReader(conn.getInputStream()));
            String line = "";
            while((line=in.readLine())!=null) {
                if (line.contains("id=\"mature\"")){
                    mirna=line.substring(line.indexOf("id=\"mature\" >")+13);
                    mirna=mirna.substring(0,mirna.indexOf("<"));
                } 
                if (line.contains("big_num")){
                    String tmp=line.substring(line.indexOf("class=\"big_num\">Mature")+38);
                    mirnaStart=Integer.valueOf(tmp.substring(0,tmp.indexOf("<")));
                }
            }
            in.close();   
	} catch(Exception e){
             e.printStackTrace();       
        }
        return mirnaStart;
    }
   


}
