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
 * Alignment object
 */
package miRdupPredictor;

import miRdup.Features;

/**
 *
 * @author fyrox
 */
public class AlignmentObject {
    
    private String id;
    private Double score;
    private String alignment;
    private String prec;
    private String struct;
    private String knownmirna;
    private String phrase;
    private String arm;
    private int start;
    private int end;
    private String mirna;

    public AlignmentObject() {
        
    }

    /**
     * @return the id
     */
    public String getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(String id) {
        this.id = id;
    }


    /**
     * @return the alignment
     */
    public String getAlignment() {
        return alignment;
    }

    /**
     * @param alignment the alignment to set
     */
    public void setAlignment(String alignment) {
        this.alignment = alignment;
        String predictedmirna=alignment.trim();
        this.setStart(this.getPrec().indexOf(predictedmirna));
        this.setEnd(this.getPrec().indexOf(predictedmirna)+predictedmirna.length());
        this.setMirna(getPrec().substring(this.getStart(), this.getEnd()));
        this.setArm(getArm(getMirna(),getPrec(),getStruct()));
    }
    
    public  String getArm(String mirna, String prec, String struct){
        int taille = mirna.length();
        int mirnaStart=prec.indexOf(mirna);
        int mirnaEnd=mirnaStart+taille;
        if (getMirnaIncludedInLoop(struct)) {
            return "5'";
        }
        String arm=null;
        try {
            if (struct.substring(mirnaStart, mirnaEnd).contains("(")) {
                arm = "5'";
            } else {
                arm = "3'";
            }
        } catch (Exception e) {
            return "5'";
        }
        return arm;
    }
    
     /**
     * Vérifie que le miRNA n'est pas complètement inclu dans la loop
     * @param prec
     * @param struct
     * @param mirna
     * @return true si le miRNA est dans la loop
     */
    public  boolean getMirnaIncludedInLoop(String precStruc){        
        try {
            int startl = prec.indexOf(mirna);
            int endl = startl + mirna.length();
            if (endl > prec.length()||startl==-1) {
                return false;
            }
            String precStrucSubs=precStruc.substring(startl, endl);
            if (precStrucSubs.contains("(")&&precStrucSubs.contains(")")) {                
                return true;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return false;
    }
    
    /**
     * @return the score
     */
    public Double getScore() {
        return score;
    }

    /**
     * @param score the score to set
     */
    public void setScore(Double score) {
        this.score = score;
    }
    
    @Override
    public String toString(){
        int spaces=getPrec().length()-getAlignment().length();
        String space="";
        for (int i = 0; i <= spaces; i++) {
            space+=" ";            
        }
        
        String s=""
                + getAlignment()
                + space
                + getScore();
        return s;
    }

    /**
     * @return the prec
     */
    public String getPrec() {
        return prec;
    }

    /**
     * @param prec the prec to set
     */
    public void setPrec(String prec) {
        this.prec = prec;
    }

    /**
     * @return the phrase
     */
    public String getPhrase() {
        return phrase;
    }

    /**
     * @param phrase the phrase to set
     */
    public void setPhrase(String phrase) {
        this.phrase = phrase;
    }

    /**
     * @return the knownmirna
     */
    public String getKnownmirna() {
        return knownmirna;
    }

    /**
     * @param knownmirna the knownmirna to set
     */
    public void setKnownmirna(String knownmirna) {
        this.knownmirna = knownmirna;
    }

    /**
     * @return the start
     */
    public int getStart() {
        return start;
    }

    /**
     * @param start the start to set
     */
    public void setStart(int start) {
        this.start = start;
    }

    /**
     * @return the end
     */
    public int getEnd() {
        return end;
    }

    /**
     * @param end the end to set
     */
    public void setEnd(int end) {
        this.end = end;
    }

    /**
     * @return the mirna
     */
    public String getMirna() {
        return mirna;
    }

    /**
     * @param mirna the mirna to set
     */
    public void setMirna(String mirna) {
        this.mirna = mirna;
    }

    /**
     * @return the struct
     */
    public String getStruct() {
        return struct;
    }

    /**
     * @param struct the struct to set
     */
    public void setStruct(String struct) {
        this.struct = struct;
    }

    /**
     * @return the arm
     */
    public String getArm() {
        return arm;
    }

    /**
     * @param arm the arm to set
     */
    public void setArm(String arm) {
        this.arm = arm;
    }

    
    
}
