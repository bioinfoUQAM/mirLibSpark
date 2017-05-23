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
 * mirna object
 */
package miRdup;

import java.text.DecimalFormat;

/**
 * 
 * @author Mickael Leclercq
 */
public final class MirnaObject {
    
    private int id;
    private String idName;
    private String matureSequence;
    private String precursorSequence;
    private String structure;
    private String structureMotif;
    private String constraint;
    private String shortName; //used for mirbase
    private String fullName; //used for mirbase
    private String arff;    
    private boolean experimental;
    private boolean star;

    
    //predictions
    private int actual;
    private int predicted;
    private double score;
    private boolean validated;
    private boolean needPrediction;
    private String predictedmiRNA;
    private String predictedmiRNAstar;
    private String error;
    static DecimalFormat dec = new DecimalFormat();
    
    public MirnaObject(){
        dec.setMaximumFractionDigits(2);
        
    }
    
    MirnaObject(String fullName, String shortName, String matureSequence, String precursorSequence){
        this.setFullName(fullName);
        this.setMatureSequence(matureSequence);
        this.setShortName(shortName);    
        this.setPrecursorSequence(precursorSequence);        
    }
    
    /**
     * @return the matureSequence
     */
    public String getMatureSequence() {
        return matureSequence;
    }

    /**
     * @param matureSequence the matureSequence to set
     */
    public void setMatureSequence(String matureSequence) {
        this.matureSequence = matureSequence;
    }

    /**
     * @return the precursor
     */
    public String getPrecursorSequence() {
        return precursorSequence;
    }

    /**
     * @param precursorSequence the precursor to set
     */
    public void setPrecursorSequence(String precursorSequence) {
        this.precursorSequence = precursorSequence;        
    }

    /**
     * @return the structure
     */
    public String getStructure() {
        return structure;
    }

    /**
     * @param structure the structure to set
     */
    public void setStructure(String structure) {
        this.structure = structure;
    }

    /**
     * @return the shortName
     */
    public String getShortName() {
        return shortName;
    }

    /**
     * @param shortName the shortName to set
     */
    public void setShortName(String shortName) {
        this.shortName = shortName;
    }

    /**
     * @return the fullName
     */
    public String getFullName() {
        return fullName;
    }

    /**
     * @param fullName the fullName to set
     */
    public void setFullName(String fullName) {
        this.fullName = fullName;
        if (fullName.endsWith("-3p")){
            this.setStar(true);
        }
    }

    
    /**
     * @return the id
     */
    public int getId() {
        return id;
    }

    /**
     * @param id the id to set
     */
    public void setId(int id) {
        this.id = id;
    }
    
    

    /**
     * @return the actual
     */
    public int getActual() {
        return actual;
    }

    /**
     * @param actual the actual to set
     */
    public void setActual(int actual) {
        this.actual = actual;
    }

    /**
     * @return the predicted
     */
    public int getPredicted() {
        return predicted;
    }

    /**
     * @param predicted the predicted to set
     */
    public void setPredicted(int predicted) {
        this.predicted = predicted;
    }

    /**
     * @return the score
     */
    public double getScore() {
        return score;
    }

    /**
     * @param score the score to set
     */
    public void setScore(double score) {
        this.score = score;
    }

    /**
     * @return the validated
     */
    public boolean isValidated() {
        return validated;
    }

    /**
     * @param validated the validated to set
     */
    public void setValidated(boolean validated) {
        this.validated = validated;
    }


    /**
     * @return the idName
     */
    public String getIdName() {
        return idName;
    }

    /**
     * @param idName the idName to set
     */
    public void setIdName(String idName) {
        this.idName = idName;
    }

    /**
     * @return the experimental
     */
    public boolean isExperimental() {
        return experimental;
    }

    /**
     * @param experimental the experimental to set
     */
    public void setExperimental(boolean experimental) {
        this.experimental = experimental;
    }

    /**
     * @return the arff
     */
    public String getArff() {
        return arff;
    }

    /**
     * @param arff the arff to set
     */
    public void setArff(String arff) {
        this.arff = arff;
    }
    
    
    @Override
    public String toString(){
        String s=""+
                this.getFullName()+"\t"+
                this.getShortName()+"\t"+
                this.getMatureSequence()+"\t"+
                this.getPrecursorSequence()+"\t"+
                this.getStructure()+"\t"+
                this.isValidated()+"\t"+
                this.getScore();
        
        return s;
    }
    
    public String toStringTXT(){
        String s=""+
                this.getFullName()+"\t"+
                this.getShortName()+"\t"+
                this.getMatureSequence()+"\t"+
                this.getPrecursorSequence()+"\t"+
                this.getStructure();
        
        return s;
    }    
    
    public String toStringPredictions(){        
        String s=""+
                this.getIdName()+"\t"+
                this.getMatureSequence()+"\t"+
                this.getPrecursorSequence()+"\t"+
                this.getStructure()+"\t"+
                this.isValidated()+"\t"+                
                dec.format(this.getScore());
        try {
            if (this.isNeedPrediction() && this.isValidated() == false) {
                s += "\t" + this.getPredictedmiRNA() + "\t"
                        + this.getPredictedmiRNAstar();
            }
        } catch (Exception e) {
            s += "\tno predictable miRNA\tno predictable miRNAstar";
        }
        return s;
    }
    
    public String toStringFullPredictions(){
        String arff[]=this.getArff().split(",");
        
        String s="\n--------\n"+
                "#ID\t"+this.getIdName()+"\n"+
                "#WE\t"+this.getArff()+"\n"+
//                "#FF\tmfe\t"+arff[0]+"\n"+
//                "#FF\tMaximumLengthWithoutBulges\t"+arff[1]+"\n"+
//                "#FF\tMaximumLengthWithoutBulgesPerc\t"+arff[2]+"\n"+
//                "#FF\tBasePairsInDuplex\t"+arff[3]+"\n"+
//                "#FF\tStartOfPerfect5MerBasePair\t"+arff[4]+"\n"+
//                "#FF\tDistanceFromTerminalLoop\t"+arff[5]+"\n"+
//                "#FF\tDistanceFromHairpinStart\t"+arff[6]+"\n"+
//                "#FF\tMirnaIncludedInLoop\t"+arff[7]+"\n"+
//                "#FF\tLengthOfOverlapInLoop\t"+arff[8]+"\n"+
//                "#FF\tAverageNumberOfPairedBasesInWindow7\t"+arff[9]+"\n"+
//                "#FF\tAverageNumberOfPairedBasesInWindow5\t"+arff[10]+"\n"+
//                "#FF\tAverageNumberOfPairedBasesInWindow3\t"+arff[11]+"\n"+
//                "#FF\tBulgeAtPosition3\t"+arff[12]+"\n"+
//                "#FF\tLenghtOfBiggestBulge\t"+arff[13]+"\n"+
//                "#FF\tLengthBiggestBulgesPerc\t"+arff[14]+"\n"+
//                "#FF\tA...\t"+arff[15]+"\n"+
//                "#FF\tC...\t"+arff[16]+"\n"+
//                "#FF\tG...\t"+arff[17]+"\n"+
//                "#FF\tU...\t"+arff[18]+"\n"+
//                "#FF\tPercBasedPairAU\t"+arff[19]+"\n"+
//                "#FF\tPercBasedPairGC\t"+arff[20]+"\n"+
//                "#FF\tPercBasedPairGU\t"+arff[21]+"\n"+
                "#PR\tPrediction\t"+this.isValidated()+"\n"+
                "#SC\tScore\t"+this.getScore()+"\n"+
                Alignment()+"\n"
                ;        
        return s;
    }
    
    public String Alignment(){
        String mirna=this.getMatureSequence();
        String struc=this.getStructure();
        String prec=this.getPrecursorSequence();
        
        int mirnaStart=prec.indexOf(mirna);
        String spaces="";
        for (int i = 0; i < mirnaStart; i++) {
            spaces+=" ";            
        }
        String aln=""+
                "#PM\t"+prec+"\n"+  // pre-miRNA
                "#SS\t"+struc+"\n"+  // secondary structure
                "#MM\t"+spaces+mirna+"\n" //mature miRNA
                ;
        try {
            if (this.isNeedPrediction() && this.isValidated() == false) {
                String predmirna = this.getPredictedmiRNA();                
                int predmirnaStart = prec.indexOf(predmirna);
                spaces = "";
                for (int i = 0; i < predmirnaStart; i++) {
                    spaces += " ";                    
                }
                aln += "#PM\t" + spaces + predmirna + "\n"; //Predicted mirna
                
                String predmirnastar = this.getPredictedmiRNAstar();
                int predmirnastarStart = prec.indexOf(predmirnastar);
                spaces = "";
                for (int i = 0; i < predmirnastarStart; i++) {
                    spaces += " ";                    
                }
                aln += "#PS\t" + spaces + predmirnastar + "\n"; //Predicted mirnaStar            
            }
        } catch (Exception e) {
            aln += "#PM\tno predictable miRNA\n"; 
            aln += "#PS\tno predictable miRNAStar\n"; 
        }
        
        return aln;
    }

    /**
     * @return the constraint
     */
    public String getConstraint() {
        int mirnaStart = precursorSequence.indexOf(matureSequence);
        int mirnaEnd = mirnaStart+matureSequence.length();
        boolean left;
        if (mirnaStart<precursorSequence.length()-mirnaEnd){
            left=true;
        } else left=false;
        
        String c="";
        for (int i = 0; i < precursorSequence.length(); i++) {
            if (i>=mirnaStart&&i<mirnaEnd){
                if (left){
                    c+="<";
                } else
                    c+=">";
            } else {
                c+=".";
            }
            
        }
        constraint=c;
        
        return constraint;
    }

    /**
     * @param constraint the constraint to set
     */
    public void setConstraint(String constraint) {
        this.constraint = constraint;
    }

    /**
     * @return the structureMotif
     */
    public String getStructureMotif() {
        return structureMotif;
    }

    /**
     * @param structureMotif the structureMotif to set
     */
    public void setStructureMotif(String structureMotif) {
        this.structureMotif = structureMotif;
    }

    /**
     * @return the star
     */
    public boolean isStar() {
        return star;
    }

    /**
     * @param star the star to set
     */
    public void setStar(boolean star) {
        this.star = star;
    }

    /**
     * @return the predictedmiRNA
     */
    public String getPredictedmiRNA() {
        return predictedmiRNA;
    }

    /**
     * @param predictedmiRNA the predictedmiRNA to set
     */
    public void setPredictedmiRNA(String predictedmiRNA) {
        this.predictedmiRNA = predictedmiRNA;
    }

    /**
     * @return the predictedmiRNAstar
     */
    public String getPredictedmiRNAstar() {
        return predictedmiRNAstar;
    }

    /**
     * @param predictedmiRNAstar the predictedmiRNAstar to set
     */
    public void setPredictedmiRNAstar(String predictedmiRNAstar) {
        this.predictedmiRNAstar = predictedmiRNAstar;
    }

    /**
     * @return the needPrediction
     */
    public boolean isNeedPrediction() {
        return needPrediction;
    }

    /**
     * @param needPrediction the needPrediction to set
     */
    public void setNeedPrediction(boolean needPrediction) {
        this.needPrediction = needPrediction;
    }

    /**
     * @return the error
     */
    public String getError() {
        return error;
    }

    /**
     * @param error the error to set
     */
    public void setError(String error) {
        this.error = error;
    }

}
