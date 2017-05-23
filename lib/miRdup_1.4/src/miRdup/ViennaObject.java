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
 * Object holding Vienna package tools informations
 */
package miRdup;

import java.io.BufferedReader;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author mycky
 */
public class ViennaObject {
    public String mirna5p;
    public String mirna3p;
    public String bothmirnas;
    public String RNAcofoldStructure;
    public String RNAcofoldMfe;
    public String RNAcofoldStructureEnsemble;
    public String RNAcofoldMfeEnsemble;
    public String RNAcofoldFrequency;
    public String RNAcofoldDeltaG;
    public String RNAcofoldAB;
    public String RNAcofoldAA;
    public String RNAcofoldBB;
    public String RNAcofoldA;
    public String RNAcofoldB;
    
    public Double BPprobGlobal=0.0; // sum of probabilities of ALL possible pairs
    public Double BPprobTmpStructure=0.0; // sum of probabilities of possible pairs
    public HashMap<Integer,Double[]> hmBPprobTmpStructure= new HashMap<Integer,Double[]>();
    public Double BPprobFinalStructure=0.0; // sum of probabilities of the final structure
    public ArrayList<Double[]> alBPprobFinalStructure= new ArrayList<Double[]>();
    public double BPprobEnsemblizedBulges;
    public double ComplexBoltzmannProbability;
     
    public Double BPprobFinalStructureGC=0.0; //canonical
    public Double BPprobFinalStructureGU=0.0; //canonical
    public Double BPprobFinalStructureAU=0.0; //canonical
    
    //non canonical (USELESS, never happen)
    public Double BPprobFinalStructureGA=0.0;
    public Double BPprobFinalStructureGG=0.0;    
    public Double BPprobFinalStructureCC=0.0;
    public Double BPprobFinalStructureCU=0.0;
    public Double BPprobFinalStructureCA=0.0;
    public Double BPprobFinalStructureAA=0.0;    
    public Double BPprobFinalStructureUU=0.0;
    
    public Boolean error=false;
    public boolean calculateConstraint=true;
    public static DecimalFormat df = new DecimalFormat ( ) ;
    

    public void parseRNAcofold(String output) {
        df.setMaximumFractionDigits(3);
        String tab[]= output.split("\n");
        mirna5p=tab[0].split("&")[0];
        mirna3p=tab[0].split("&")[1];
        bothmirnas=mirna5p+mirna3p;
        RNAcofoldStructure=tab[1].substring(0,tab[1].lastIndexOf("("));
        RNAcofoldMfe=tab[1].substring(tab[1].lastIndexOf("(")+1, tab[1].lastIndexOf(")")).trim();
        RNAcofoldStructureEnsemble=tab[2].substring(0,tab[1].lastIndexOf("("));
        RNAcofoldMfeEnsemble=tab[2].substring(tab[2].indexOf("[")+1, tab[2].indexOf("]")).trim();
        RNAcofoldFrequency=tab[3].substring(tab[3].indexOf("ble")+3, tab[3].indexOf(",")-1).trim();
        RNAcofoldDeltaG=tab[3].substring(tab[3].indexOf("=")+1).trim();
        String tab2[]=tab[6].split("\t");
        RNAcofoldAB=tab2[0];
        RNAcofoldAA=tab2[1];
        RNAcofoldBB=tab2[2];
        RNAcofoldA=tab2[3];
        RNAcofoldB=tab2[4];
        getBasePairsProbabilities("ABdot5.ps");
        System.out.print("");
    }
    
    
    public void parseRNAcofoldShell(String output) {
        String tab[]= output.split("\n");
        RNAcofoldMfe=tab[4].substring(tab[4].indexOf("=")+1, tab[4].indexOf("kcal")-1).trim();
        RNAcofoldMfeEnsemble=tab[6].substring(tab[6].indexOf("=")+1, tab[6].indexOf("kcal")-1).trim();
        RNAcofoldFrequency=tab[7].substring(tab[7].indexOf("ble")+3, tab[7].indexOf(",")-1).trim();
        RNAcofoldDeltaG=tab[8].substring(tab[8].indexOf("=")+1).trim();
        String tab2[]=tab[10].split("\t");
        RNAcofoldAB=tab2[0];
        RNAcofoldAA=tab2[1];
        RNAcofoldBB=tab2[2];
        RNAcofoldA=tab2[3];
        RNAcofoldB=tab2[4];
        getBasePairsProbabilities("ABdot5.ps");
    }

    
    public String toStringHeaderRNAcofold(){
        String s = ""
                + "RNAcofoldMfe"+"\t"
                + "RNAcofoldMfeEnsemble"+"\t"
                + "RNAcofoldFrequency"+"\t"
                + "RNAcofoldDeltaG"+"\t"
                + "RNAcofoldAB"+"\t"
                + "RNAcofoldAA"+"\t"
                + "RNAcofoldBB"+"\t"
                + "RNAcofoldA"+"\t"
                + "RNAcofoldB"+"\t"
                + "";
        return s;
        
    }
    
    public String toStringRNAcofold(){
        String s = ""
                + RNAcofoldMfe+","
                //+ RNAcofoldMfeEnsemble+"," //useless, same as AB
                + df.format(Double.valueOf(RNAcofoldFrequency))+","
                + df.format(Double.valueOf(RNAcofoldDeltaG))+","
                + df.format(Double.valueOf(RNAcofoldAB))+","
                + df.format(Double.valueOf(RNAcofoldAA))+","
                + df.format(Double.valueOf(RNAcofoldBB))+","
                + df.format(Double.valueOf(RNAcofoldA))+","
                + df.format(Double.valueOf(RNAcofoldB))+","
                + df.format(Double.valueOf(BPprobGlobal))+","
                + df.format(Double.valueOf(BPprobTmpStructure))+","
                + df.format(Double.valueOf(BPprobFinalStructure))+","
                + df.format(Double.valueOf(BPprobEnsemblizedBulges))+","
                + df.format(Double.valueOf(ComplexBoltzmannProbability))+","
                + df.format(Double.valueOf(BPprobFinalStructureGC))+","
                + df.format(Double.valueOf(BPprobFinalStructureGU))+","
                + df.format(Double.valueOf(BPprobFinalStructureAU))+""
//                + BPprobFinalStructureGG+","
//                + BPprobFinalStructureCC+","
//                + BPprobFinalStructureCU+","
//                + BPprobFinalStructureCA+","
//                + BPprobFinalStructureAU+","
//                + BPprobFinalStructureGA+","
//                + BPprobFinalStructureUU+""
                ;
        return s;        
    }

    
    public String toStringError() {
        String s = ""
                + "0"+","
                //+ "0"+","
                + "0"+","
                + "0"+","
                + "0"+","
                + "0"+","
                + "0"+","
                + "0"+","
                + "0"+","
                
                + "0.0"+","
                + "0.0"+","
                + "0.0"+","
                + "0.0"+","
                + "0.0"+","
                
                + "0.0"+","
                + "0.0"+","
                + "0.0"+""
                
//                + "0.0"+"\t"
//                + "0.0"+"\t"                
//                + "0.0"+"\t"
//                + "0.0"+"\t"
//                + "0.0"+"\t"                
//                + "0.0"+"\t"
//                + "0.0"+"\t"                
//                + "0.0"
                
                + "";
        return s;
    }

    /**
     * if submitted sequence contains errors
     * @return 
     */
    public boolean hasError(){
        if (error==true){
            return true;            
        } else {
            return false;
        }
    }

    /**
     * Get probabilities for each pairs stored in the ps file
     * @param aBdot5ps 
     */
    private void getBasePairsProbabilities(String aBdot5ps) {

        boolean finalProbs=false;
        int nucleotide5p=0;
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(aBdot5ps));
            String line=br.readLine();
            while(!line.startsWith("%data")){
                line=br.readLine();
            }
            line=br.readLine();
            while (!line.startsWith("showpage")){
                String tab[]=line.split(" ");
                int actual5pnt=Integer.valueOf(tab[0])-1;
                int nucleotide3p=Integer.valueOf(tab[1])-1;
                double score=Double.valueOf(tab[2]);
                Double content[]= {Double.valueOf(actual5pnt),Double.valueOf(nucleotide3p),score};
                
                if (actual5pnt<nucleotide5p){
                    finalProbs=true;
                }
                nucleotide5p=actual5pnt;
//                System.out.println(line);
                assignBasePairsProbability(bothmirnas.charAt(actual5pnt), 
                        bothmirnas.charAt(nucleotide3p), score);
                
                BPprobGlobal+=score;
                if (finalProbs){
                    BPprobFinalStructure+=score;                    
                    alBPprobFinalStructure.add(content);
                }else {
                    BPprobTmpStructure+=score;
                    hmBPprobTmpStructure.put(actual5pnt,content);
                    ComplexBoltzmannProbability+=-Math.exp(-score/(8.3144621*37));
                }
                
                line=br.readLine();                
            }
            br.close();
            
            //ensemblized bulges at +4 and -4
            int numberOfUnpairedBases=0;
            double SumUnpairedBasesAroundBulges=0;
            for (int i = 0; i < alBPprobFinalStructure.size(); i++) {
                if (alBPprobFinalStructure.get(i)[0]!=i){  
                    numberOfUnpairedBases++;
                    for (int j = i-4; j <= i+4; j++) {
                        try {

                            SumUnpairedBasesAroundBulges+=hmBPprobTmpStructure.get(j)[2];
                        } catch (Exception e) {
                        }
                    }
                }
            }
            BPprobEnsemblizedBulges=SumUnpairedBasesAroundBulges/numberOfUnpairedBases;
            
            
            // ComplexBoltzmannProbability
            if (calculateConstraint) {
                String constraint = "";
                if (mirna5p.length() == mirna3p.length()) {
                    for (int i = 0; i < mirna5p.length(); i++) {
                        constraint += "(";
                    }
                    constraint += "&";
                    for (int i = 0; i < mirna3p.length(); i++) {
                        constraint += ")";
                    }
                } else if (mirna5p.length() < mirna3p.length()) {
                    for (int i = 0; i < mirna5p.length(); i++) {
                        constraint += "(";
                    }
                    constraint += "&";
                    for (int i = 0; i < mirna5p.length(); i++) {
                        constraint += ")";
                    }
                    for (int i = 0; i < mirna5p.length() - mirna5p.length(); i++) {
                        constraint += ".";
                    }
                } else {
                    for (int i = 0; i < mirna5p.length() - mirna3p.length(); i++) {
                        constraint += ".";
                    }
                    for (int i = 0; i < mirna3p.length(); i++) {
                        constraint += "(";
                    }
                    constraint += "&";
                    for (int i = 0; i < mirna3p.length(); i++) {
                        constraint += ")";
                    }
                    
                }
                
                double Q = Vienna.GetInfosDuplexRNAcofoldConstraint(mirna5p, mirna3p, constraint).ComplexBoltzmannProbability;
                ComplexBoltzmannProbability = ComplexBoltzmannProbability / Q;
            }
            
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void assignBasePairsProbability(char fiveP,char threeP, Double score) {
         //canonical
        if ((fiveP=='A'&&threeP=='U')||(fiveP=='U'&&threeP=='A')){
            BPprobFinalStructureAU+=score;
        }
        if ((fiveP=='G'&&threeP=='U')||(fiveP=='U'&&threeP=='G')){
            BPprobFinalStructureGU+=score;
        }
        if ((fiveP=='G'&&threeP=='C')||(fiveP=='C'&&threeP=='G')){
            BPprobFinalStructureGC+=score;
        }
        
         //non canonical
        if ((fiveP=='G'&&threeP=='A')||(fiveP=='A'&&threeP=='G')){
            BPprobFinalStructureGA+=score;
        }
        if ((fiveP=='C'&&threeP=='A')||(fiveP=='A'&&threeP=='C')){
            BPprobFinalStructureCA+=score;
        }
        if ((fiveP=='C'&&threeP=='U')||(fiveP=='U'&&threeP=='C')){
            BPprobFinalStructureCU+=score;
        }
        if ((fiveP=='G'&&threeP=='G')){
            BPprobFinalStructureGG+=score;
        }
        if ((fiveP=='C'&&threeP=='C')){
            BPprobFinalStructureCC+=score;
        }
        if ((fiveP=='A'&&threeP=='A')){
            BPprobFinalStructureAA+=score;
        }
        if ((fiveP=='U'&&threeP=='U')){
            BPprobFinalStructureUU+=score;
        }

    }
    
}
