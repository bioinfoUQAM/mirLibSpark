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
 * extract results from mirdeep 2 
 * make them compatible for mirdup
 */
package paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import miRdup.MirnaObject;
/**
 *
 * @author Mickael
 */
public class mirdeep2OutputTomiRdup {
    
    public static void main(String[] args) {
        String f = "";
        File infile = new File(f+"SRR029124.output.mirdeep2.mrd");
        File outfile = new File(f+"SRR0239124.mirdeep2.toPredict.txt");
        
        // Add precursors to Hashmap
        HashMap<String,String> hmPrec = new HashMap<String,String>();        
        try {
            BufferedReader br = new BufferedReader(new FileReader("miRNA.dat"));
            String line="";
            line = br.readLine();
            while (br.ready()){                
                if (line.startsWith("ID")){
                    String name=line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                    String seq="";
                    while (!line.startsWith("SQ")){                        
                        line=br.readLine();                        
                    }  
                    line=br.readLine(); 
                    while(!line.startsWith("//")){                        
                        seq+=line.substring(0,75).replace(" ", "").toUpperCase();
                        line=br.readLine();
                    }
                    hmPrec.put(name.toLowerCase(), seq);
                } else {
                    line = br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        // add mirnas in objects Arraylist
        HashMap<String,MirnaObject> hmmirnas=new HashMap<String, MirnaObject>();
        int cpt=0;
        try {
            BufferedReader br = new BufferedReader(new FileReader("miRNA.dat"));
            String line="";
            
            while (br.ready()){
                line=br.readLine();
                if (line.startsWith("ID")){
                    String species=line.substring(5,8);                    
                    String shortname= line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                    while (!line.startsWith("FT")){
                        line=br.readLine();
                    }
                    while (!line.startsWith("XX")){
                        if (line.contains("FT   miRNA")){
                            MirnaObject m = new MirnaObject();
                            m.setId(cpt++);
                            m.setShortName(shortname);
                            //mirna
                            String positions=line.substring(21);                                
                            int start=Integer.valueOf(positions.split("\\.")[0])-1;
                            int end=Integer.valueOf(positions.split("\\.")[2]);                                
                            String prec = hmPrec.get(m.getShortName());
                            
                            m.setPrecursorSequence(prec);
                            String mirna=prec.substring(start,end);
                            m.setMatureSequence(mirna);
                            
                            br.readLine(); // accession
                            line=br.readLine();// product
                            String product=line.substring(line.indexOf("=")+1).replace("\"", "");
                            m.setFullName(product); // set star in the same time
                            line=br.readLine();// evidence
                            boolean experimental=true;
                            if (line.contains("not")){
                                experimental=false;
                            }
                            m.setExperimental(experimental);     
                            hmmirnas.put(m.getFullName(), m);
                        }
                        line=br.readLine();
                    }
                }
            }
        } catch (Exception e){
            e.printStackTrace();
        }
          
        //read mirdeep2
        
        //rejected by mirdup
//        String tab[]={"chr22_8101","chr17_4694","chr10_177","chr2_8983","chr15_3241","chr8_14569","chr3_9885","chrX_15880","chr2_9254"};
//        HashMap<String,Integer> hmTmp = new HashMap<String, Integer>();
//        for (String s : tab) {
//            hmTmp.put(s, 0);
//        }
        
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            PrintWriter pw = new PrintWriter(new FileWriter(outfile,false));
            String line="";
            String name="";
            Boolean conserved=false;
            while (br.ready()){               
                if (line.startsWith(">")){
                    name=line.substring(1);
                    while (!line.startsWith("score for cons. seed")){
                        line=br.readLine();
                    }
                    //get conserved mirna
                    String conservedmirna = null;                    
                    line=br.readLine();         
                    if (line.startsWith("miRNA")){
                        conservedmirna=line.split("\t")[1].trim(); 
                        conserved=true;
                    } else {
                        conserved=false;
                    }
                    while (!line.startsWith("exp")){
                        line = br.readLine();
                    }
                                        
                    int mirnastart = line.indexOf("M");
                    int mirnaend= line.lastIndexOf("M")+1;
                    int start = line.indexOf("f");
                    line = br.readLine();
                    if (line.startsWith("obs")){
                        line=br.readLine();
                    }
                    String precursor = line.toUpperCase();
                    String mirna = precursor.substring(mirnastart, mirnaend);    
                    String struct = br.readLine();
                    String cons="";
                    String startchr;
                    String endchr;
                    String chr;
                    if (conserved){
                        MirnaObject mo= hmmirnas.get(conservedmirna);
                        String originalmirna = null;
                        boolean experimental=false;
                        if (mo!=null) {
                            originalmirna = mo.getMatureSequence();
                            experimental=mo.isExperimental();
                            
                        }
                        if (mirna.equals(originalmirna)/*&&hmTmp.containsKey(name)*/){
                            cons="c";
                            //cons="c."+conservedmirna;
                            
                            if (experimental){
                                cons="ce\t"+conservedmirna;
                            } else {
                                cons="c\t"+conservedmirna;
                            }
                                
                            
                            
                            //non predicted by mirdup
//                            pw.println(name+"_"+cons+"\t"+mirna+"\t"+precursor.substring(start) +"\t"+struct.substring(start,struct.lastIndexOf("\t")));
//                            pw.println(mirna+"\t"+originalmirna);
//                            pw.println(precursor.substring(start)+"\n"+mo.getPrecursorSequence());
//                            pw.println();
                        } else {
                            cons="n\t";
                        }

                    }else {
                            cons="n\t";
                        }
                    
                    pw.println(name+"_"+cons+"\t"+mirna+"\t"+precursor.substring(start) +"\t"+struct.substring(start,struct.lastIndexOf("\t")));
                    
                } else {
                    line=br.readLine();
                }
            }
            pw.flush();pw.close();
            
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        
    }
    
}
