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
 * Generate all possible miRNAs from a precursor and determine wether it is 
 * well positioned or not
 */
package paper;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.TooManyListenersException;
import miRdup.Features;
import miRdup.Main;
import miRdup.Vienna;
import miRdupPredictor.AlignmentObject;
import miRdupPredictor.AlignmentObject;
import miRdupPredictor.AlignmentObject;

/**
 *
 * @author fyrox
 */
public class PredictorCompMirbase {
    public static String rnafold="/home/mycky/tools/ViennaRNA-2.0.6/Progs/";
//    public static String rnafold="/home/2011/mlecle25/tools/ViennaRNA-2.0.7/Progs/";
//    public static String rnafold="/ibrixfs1/Data/mik/tools/ViennaRNA-2.0.7/Progs/";
    public static int loopflag=0; //flag to detect looping over the same function by recursion 
    public static int maximas=0; //When looping over the same function by recursion, this value take the second (three, and so on) maxima instead of the first
    public static boolean debug=false;
    public static int threshold=10; //decisionbyproducts
    public static int nbrproducts=20; //decisionbyproducts    
    public static String struct="";
            
    public static int minlength=16; 
    public static int maxlength=30; 
    
    public static void main(String[] args) {
//        predictorBySequence();
//        args=new String[]{"predictor"+File.separator+"Arthropoda_0.splitted.txt","Arthropoda.model","predictor"+File.separator+"Arthropoda_0.splitted.results.txt"};
        predictorByFile(args);
    }

    public static void predictorBySequence(){
        String konwnMirna = "CAAAUUCGGUUCUAGAGAGGUUU";
        String prec = "CCACGUCUACCCUGUAGAUCCGAAUUUGUUUUAUACUAGCUUUAAGGACAAAUUCGGUUCUAGAGAGGUUUGUGUGG";

        String model ="Arthropoda.model";
        String outfile="test.txt";
        
        ArrayList<Integer> diffStarts= new ArrayList<Integer>();        
        ArrayList<Integer> diffEnds= new ArrayList<Integer>();
        
        for (int i = 0; i < 1000; i++) {
            diffStarts.add(0);
            diffEnds.add(0);
        }
        
        String differences = executeGenerator(konwnMirna, prec, outfile+".miRdup.predictor", model);
        //differences from known mirna
        int diffstart = Integer.valueOf(differences.split(",")[0]);
        int diffend = Integer.valueOf(differences.split(",")[1]);

        int tmp = diffStarts.get(diffstart);
        tmp += 1;
        diffStarts.set(diffstart, tmp);

        tmp = diffEnds.get(diffend);
        tmp += 1;
        diffEnds.set(diffend, tmp);
        
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
    
    
    
    public static void predictorByFile(String args[]){
//        String infile = "split\\all_0.splitted.txt";
//        String model="all.AFAdaboostRF.model";   
//        String outfile="generator\\all_0.splitted.results.txt";
        
        String infile = args[0];
        String model=args[1];   
        String outfile=args[2];
        
        ArrayList<Integer> diffStarts= new ArrayList<Integer>();        
        ArrayList<Integer> diffEnds= new ArrayList<Integer>();
        
        for (int i = 0; i < 1000; i++) {
            diffStarts.add(0);
            diffEnds.add(0);
        }
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            String line="";
            while(br.ready()){
                try {
                    loopflag=0;
                    maximas=0;
                    line = br.readLine();
                    System.out.println("\n----------");
                    System.out.println(line);
                    String konwnMirna = line.split("\t")[2];
                    String prec = line.split("\t")[3];
                    struct="";
                    String differences = executeGenerator(konwnMirna, prec, infile+".miRdup.predictor", model);                       
                    //differences from known mirna
                    int diffstart = Integer.valueOf(differences.split(",")[0]);
                    int diffend = Integer.valueOf(differences.split(",")[1]);
                    
                    int tmp = diffStarts.get(diffstart);
                    tmp += 1;
                    diffStarts.set(diffstart, tmp);
                    
                    tmp = diffEnds.get(diffend);
                    tmp += 1;
                    diffEnds.set(diffend, tmp);
                    
                    
                    
                } catch (Exception e) {
                    System.err.println("error at "+ line);
                    e.printStackTrace();
                }                 
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
        //remove tmp files
        File f = new File(infile+".miRdup.predictor");                    
        f.delete();
        f = new File(infile+".miRdup.predictor"+"."+model+".miRdup.txt");                    
        f.delete();
        f = new File(infile+".miRdup.predictor"+"."+model+".miRdup.aln.txt");                    
        f.delete();
        f = new File(infile+".miRdup.predictor.arff");                    
        f.delete();
                    
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
    
    public static String executeGenerator(String konwnMirna,String prec, String outfile, String model){
        //generate miRNAs
        generateMiRNAs(prec,outfile);
        
        //Predict miRNAs with miRdup
        predict(outfile,model);
        
        //Align true miRNAs on the precursor
        String predictionsFile=outfile+"."+model+".miRdup.txt"; // file produced by mirdup
        String alnFile=outfile+"."+model+".miRdup.aln.txt"; // output
        ArrayList<AlignmentObject> results=makeAlignment(konwnMirna,predictionsFile, alnFile);
//        return decisionByMean(results);
//        return decisionByMaximas(results);
        return decisionByWindow(results);
//        return decisionByProducts(results);
        
    }
    
    /**
     * Generate all possible miRNAs, from length 16 to 30
     * @param prec
     * @param outfile 
     */
    public static void generateMiRNAs(String prec,String outfile) {
        ArrayList<String> mirnas=new ArrayList<String>();
        String os=System.getProperty("os.name").toLowerCase();
        if (!os.startsWith("win")){
            Main.rnafoldlinux=rnafold;
        } 
        struct=Vienna.GetSecondaryStructure(prec);
        //System.out.println(prec);
        
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(outfile));
            int cpt = 0;
            for (int i = minlength; i <= maxlength; i++) {
                for (int j = 0; j <= prec.length() - i; j++) {
                    cpt++;
                    String mirna = prec.substring(j, j + i);
                    //System.out.println(mirna);
                    pw.println("mir"+cpt+"\t"+mirna+"\t"+prec+"\t"+struct);                    
                    mirnas.add(mirna);
                }                
            }
            pw.flush();pw.close();
            System.out.println("Total generated miRNAs: " + cpt);
        } catch (Exception e) {
            e.printStackTrace();
        }        
    }

    /**
     * launch miRdup
     * @param outfile
     * @param model 
     */
    private static void predict(String outfile,String model) {
        String s[]={"-c",model,"-v",outfile,"-r",rnafold};
        setOptions(s);
        miRdupExecutionEMBL();
    }
    
    private static void setOptions(String[] args) {
        Main.setOptions(args);
    }

    private static void miRdupExecutionEMBL() {
        Main.miRdupExecutionEMBL();
    }

    /**
     * Make alignment of predicted miRNAs on the pre-miRNA
     * @param infile 
     */
    private static ArrayList makeAlignment(String konwnMirna, String infile, String alnFile) {
        ArrayList<AlignmentObject> alObj = new ArrayList<AlignmentObject>();
        String struct = null;
        String prec = null;
        int tot=0;
        int trues=0;
        // read prediction file output from miRdup
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            
            String line="";
            while (br.ready()){
                line=br.readLine();
                if (line.startsWith("#ID")){
                    tot++;
                    AlignmentObject ao = new AlignmentObject();
                    ao.setId(line);
                    line=br.readLine();
                    while(line.startsWith("#FF")){
                        line=br.readLine();
                    }
                    if (line.split("\t")[2].startsWith("t")){
                        trues++;
                        ao.setScore(Double.parseDouble(br.readLine().split("\t")[2]));
                        prec=br.readLine().split("\t")[1];
                        ao.setPrec(prec);
                        ao.setKnownmirna(konwnMirna);
                        struct=br.readLine().split("\t")[1];
                        ao.setStruct(struct);
                        ao.setAlignment(br.readLine().split("\t")[1]);
                        alObj.add(ao);
                    }
                }                    
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        // get max score
        double maxScore=0.0;
        for (AlignmentObject o : alObj) {
            double score=o.getScore();
            if(score>maxScore) maxScore=score;
        }
        
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(alnFile));
            // print header
            String header = ""
                    + "\nTotal potentials miRNAs: " + tot
                    + "\nTotal predicted miRNAs: " + trues
                    + "\nBestScore= " + maxScore;
//            System.out.println(header);
//            System.out.println(prec);
//            System.out.println(struct);
            
            pw.println(header);
            pw.println(prec);
            pw.println(struct);

            // sort arraylist of miRNAs with scores
            ArrayList<String> toSort = new ArrayList<String>();
            for (AlignmentObject o : alObj) {
                String phrase = o.getScore() + "\t" + o.getId();
                o.setPhrase(phrase);
                toSort.add(phrase);
            }
            Collections.sort(toSort);
            Collections.reverse(toSort);

            // print sorted arraylist
            for (String phrase : toSort) {
                for (AlignmentObject o : alObj) {
                    if (o.getPhrase().equals(phrase)) {
                        //System.out.println(o.toString());
                        pw.println(o.toString());
                    }                    
                }
            }
            pw.flush();pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return alObj;        
    }

    /**
     * Decision by mean
     * @param results
     * @return 
     */
    private static String decisionByMean(ArrayList<AlignmentObject> results) {
        String konwnMirna = results.get(0).getKnownmirna();
        String prec=results.get(0).getPrec();
        
        ArrayList<Double> starts = new ArrayList<Double>();
        ArrayList<Double> ends = new ArrayList<Double>();
        
        // full lists with 0
        for (int i = 0; i <= prec.length(); i++) {
            starts.add(0.0);
            ends.add(0.0);
        }        
        
        int cpt=0;
        for (AlignmentObject o : results) {
            if (o.getScore()>0.99) {
                double tmp = starts.get(o.getStart());
                tmp += o.getScore();
                starts.set(o.getStart(), tmp);
                
                tmp = ends.get(o.getEnd());
                tmp += o.getScore();
                ends.set(o.getEnd(), tmp);
                cpt++;
            }
        }
        System.out.println("total="+cpt);
        
        if (debug) {
            System.out.println("Starts");
            for (int i = 0; i < starts.size(); i++) {
                System.out.println(i + "\t" + starts.get(i));                
            }
            System.out.println("\nEnds");
            for (int i = 0; i < ends.size(); i++) {
                System.out.println(i + "\t" + ends.get(i));                
            }
        }
                    
        
        //get best start and end
        double beststart = 0;
        double bestend = 0;
        try {
            beststart = Collections.max(starts);
            bestend = Collections.max(ends);
        } catch (Exception e) {
            for (Double d : starts) {
                if (d>beststart){
                    beststart=d;
                }                    
            }
            for (Double d : ends) {
                if (d>bestend){
                    bestend=d;
                }                    
            }            
        }
        
        //find corresponding positions of bestscores
        int startPosition=0;
        int endPosition=0;
        for (int i = 0; i < starts.size(); i++) {
            if (starts.get(i)==beststart){
                startPosition=i;
            }
            if (ends.get(i)==bestend){
                endPosition=i;
            }
        }        
        
        
        
        // In case where the recursion didnt remove the startPosition>endPosition problem
        // Lower the min and max with loopflag
        if (loopflag>=2){
            int loweredstart=0;
            int loweredend=0;
            for (int i = 0; i < results.size(); i++) {
                if (loweredstart!=1&&loweredend!=1) {                    
                    AlignmentObject o = results.get(i);
                    if (o.getStart() == startPosition) {
                        if (loweredstart != 1 && o.getScore() > 0.99) {
                            o.setScore(0.5);                            
                            results.set(i, o);
                            loweredstart++;
                        }
                    }
                    if (o.getEnd() == endPosition) {
                        if (loweredend != 1 && o.getScore() > 0.99) {
                            o.setScore(0.5);
                            results.set(i, o);                            
                            loweredend++;
                        }
                    }
                }
            }
        }
        
        // if start>end, we remove results from one of the arm of the hairpin
        // structure goes 5' to 3'
        if (startPosition>endPosition){
            // count 5' and 3'
            int fivePrim=0;
            int threePrim=0;
            for (AlignmentObject o : results) {                
                if (o.getArm().equals("3'")){
                    threePrim++;
                }
                if (o.getArm().equals("5'")){
                    fivePrim++;
                }
            }
            
            //remove objects from the arm having the lower count
            ArrayList<AlignmentObject> results2 = new ArrayList<AlignmentObject>();
            if (fivePrim>=threePrim){
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("5'")){
                        results2.add(o);
                    }
                }
            } else {
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("3'")){
                        results2.add(o);
                    }
                }
            }
            System.out.println("+");
            loopflag++;
            //re-execute the decision
            return decisionByMean(results2);
            
        } 
        loopflag=0;
        // Find diffstart and diffend with better accuracy
        // Sum of x times score divided by sum of scores
        
        ArrayList<Double> startScores = null;
        int whereStartBegins = 0;
        int whereStartEnds = 0;
        ArrayList<Double> endScores = null;
        int whereEndBegins = 0;
        int whereEndEnds = 0;
        
        boolean error=true; //if whereStartEnds > whereEndBegins
        while (error) {
            //search values around start
            //before max start position
            double score = -1.0;            
            startScores = new ArrayList<Double>(1000);
            //initialize
            for (int i = 0; i < starts.size(); i++) {
                startScores.add(0.0);
            }
            
            int position = startPosition;
            whereStartBegins = -1;
            whereStartEnds = -1;
            while (score != 0.0 && position != -1) {
                score = starts.get(position);
                startScores.set(position, score);
                position--;
                
            }
            whereStartBegins = position + 1;

            //after max start position
            score = -1.0;            
            position = startPosition + 1;
            
            while (score != 0 && position != starts.size()) {
                score = starts.get(position);
                startScores.set(position, score);
                position++;
                
            }
            whereStartEnds = position - 1;


            //search values around end
            //before max end position
            score = -1.0;            
            endScores = new ArrayList<Double>(1000);
            //initialize
            for (int i = 0; i < ends.size(); i++) {
                endScores.add(0.0);                
            }
            position = endPosition;
            whereEndBegins = -1;
            whereEndEnds = -1;
            
            while (score != 0 && position != -1) {
                score = ends.get(position);
                endScores.set(position, score);
                position--;
                
            }
            whereEndBegins = position + 1;
            //after max start position
            score = -1.0;            
            position = endPosition + 1;
            
            while (score != 0 && position != starts.size()) {
                score = ends.get(position);
                endScores.set(position, score);
                position++;
                
            }
            whereEndEnds = position - 1;

            // Case where whereStartEnd is superior to whereEndBegins
            if (whereStartEnds > whereEndBegins) {
                int p = whereEndBegins + ((whereStartEnds - whereEndBegins) / 2);
                starts.set(p, 0.0);
                ends.set(p, 0.0);
                error = true;
            } else {
                error=false;
            }
        }
            
        
        // do calculation for start
        int j=0;
        double SumOfScores=0.0;
        for (Double d : startScores) {
            SumOfScores+=d;
        }

        Double SumPositionsTimesScores=0.0;
        j=0;
        for (int i = whereStartBegins; i < whereStartEnds; i++) {
            SumPositionsTimesScores+=i*startScores.get(i)/SumOfScores;
            j++;
        }
        int FinalStartPosition=SumPositionsTimesScores.intValue();        
        
        //do calculation for end
        j=0;
        SumOfScores=0.0;
        for (Double d : endScores) {
            SumOfScores+=d;
        }
        
        j=0;
        SumPositionsTimesScores=0.0;
        for (int i = whereEndBegins; i < whereEndEnds; i++) {
            SumPositionsTimesScores+=i*endScores.get(i)/SumOfScores;
            j++;
        }
        int FinalEndPosition=SumPositionsTimesScores.intValue();
        
        
        // get difference from known miRNA
        int startOfKnownMiRNA=prec.indexOf(konwnMirna);
        int endOfKnownMiRNA=prec.indexOf(konwnMirna)+konwnMirna.length();
        
        int diffStart=Math.abs(FinalStartPosition-startOfKnownMiRNA);
        int diffend=Math.abs(FinalEndPosition-endOfKnownMiRNA);
        
        try {
            String consensus = prec.substring(FinalStartPosition, FinalEndPosition);

            System.out.println("Start difference from known miRNA: " + diffStart);
            System.out.println("End difference from known miRNA: " + diffend);
            
            
            System.out.println("Predicted consensus miRNA: " + consensus);

            // verify miRNA length
            if (consensus.length() < 16 || consensus.length() > 33) {
                System.out.println("Predicted miRNA too long, trying maximas method...");
                return decisionByMaximas(results);
            } else {
                return diffStart + "," + diffend;
            }
        } catch (Exception e) {
            // case where prec.substring(FinalStartPosition, FinalEndPosition) have problems
            System.out.println("Unable to process distribution, trying maximas method");
            return decisionByMaximas(results);
        }
        
    }
    
    /**
     * Decision by getting maximum values
     * @param results
     * @return 
     */
    private static String decisionByMaximas(ArrayList<AlignmentObject> results) {
        String konwnMirna = results.get(0).getKnownmirna();
        String prec=results.get(0).getPrec();
        
        ArrayList<Double> starts = new ArrayList<Double>();
        ArrayList<Double> startsSort = new ArrayList<Double>();
        ArrayList<Double> ends = new ArrayList<Double>();
        ArrayList<Double> endsSort = new ArrayList<Double>();
        
        // full lists with 0
        for (int i = 0; i <= prec.length(); i++) {
            starts.add(0.0);            
            ends.add(0.0);
            startsSort.add(0.0);            
            endsSort.add(0.0);
        }        
        
        int cpt=0;
        for (AlignmentObject o : results) {
            if (o.getScore()>0.99) {
                double tmp = starts.get(o.getStart());
                tmp += o.getScore();
                starts.set(o.getStart(), tmp);
                startsSort.set(o.getStart(), tmp);
                
                tmp = ends.get(o.getEnd());
                tmp += o.getScore();
                ends.set(o.getEnd(), tmp);
                endsSort.set(o.getEnd(), tmp);
                cpt++;
            }
        }
        System.out.println("total="+cpt);
        
        if (debug) {
            System.out.println("Starts");
            for (int i = 0; i < starts.size(); i++) {
                System.out.println(i + "\t" + starts.get(i));                
            }
            System.out.println("\nEnds");
            for (int i = 0; i < ends.size(); i++) {
                System.out.println(i + "\t" + ends.get(i));                
            }
        }
        
        // order starts and ends by score
        Collections.sort(endsSort);
        Collections.reverse(endsSort);
        Collections.sort(startsSort);
        Collections.reverse(startsSort);
        
        //get best start and end
        double beststart = startsSort.get(0+maximas);
        double bestend = endsSort.get(0+maximas);
        
        
        //find corresponding positions of bestscores
        int startPosition=0;
        int endPosition=0;
        for (int i = 0; i < starts.size(); i++) {
            if (starts.get(i)==beststart){
                startPosition=i;
            }
            if (ends.get(i)==bestend){
                endPosition=i;
            }
        }        
        
        
        
        // if start>end, we remove results from one of the arm of the hairpin
        // structure goes 5' to 3'
        if (startPosition>endPosition){
            System.out.println("Selected miRNA have its start position>to end position... "
                    + "\nRemoving results from the arm having the lowest scores");
            // count 5' and 3'
            int fivePrim=0;
            int threePrim=0;
            for (AlignmentObject o : results) {                
                if (o.getArm().equals("3'")){
                    threePrim++;
                }
                if (o.getArm().equals("5'")){
                    fivePrim++;
                }
            }
            
            //remove objects from the arm having the lower count
            ArrayList<AlignmentObject> results2 = new ArrayList<AlignmentObject>();
            if (fivePrim>=threePrim){
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("5'")){
                        results2.add(o);
                    }
                }
            } else {
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("3'")){
                        results2.add(o);
                    }
                }
            }
            
            //re-execute the decision
            loopflag++;
            if (loopflag>1){ //if we already pass through here, put maximas++
                maximas++;
            }      
            
            return decisionByMaximas(results2);
            
        } else {
            // get difference from known miRNA
            int startOfKnownMiRNA=prec.indexOf(konwnMirna);
            int endOfKnownMiRNA=prec.indexOf(konwnMirna)+konwnMirna.length();

            int diffStart=Math.abs(startPosition-startOfKnownMiRNA);
            int diffend=Math.abs(endPosition-endOfKnownMiRNA);

            System.out.println("Start difference from known miRNA: "+diffStart);
            System.out.println("End difference from known miRNA: "+diffend);

            String consensus=prec.substring(startPosition, endPosition);
            System.out.println("predicted consensus miRNA: "+consensus);
            Features f = new Features(consensus, prec, struct, true);
            String complement=f.getComplementaritySequence();
            System.out.println("predicted consensus miRNA complement: "+complement);
            
            // verify miRNA length
            if (consensus.length()<16||consensus.length()>33){
                System.out.println("Predicted miRNA too long, trying another maxima threshold...");
                maximas++;
                loopflag++;
                if (loopflag>10){
                    System.out.println("Prediction impossible, return best scored predicted miRNA");                                        
                    return results.get(1).getMirna();
                } else {
                    return decisionByMaximas(results);
                }                
            } else {
                return diffStart+","+diffend;
            }
        }
        
    }
    /**
     * Decision by sliding window
     * @param results
     * @return 
     */
    private static String decisionByWindow(ArrayList<AlignmentObject> results) {
        String knownMirna = results.get(0).getKnownmirna();
        String prec=results.get(0).getPrec();
        String structure=results.get(0).getStruct();
        ArrayList<Double> starts = new ArrayList<Double>();
        ArrayList<Double> ends = new ArrayList<Double>();
        
        // full lists with 0
        for (int i = 0; i <= prec.length(); i++) {
            starts.add(0.0);            
            ends.add(0.0);
        }        
        
        int cpt=0;
        for (AlignmentObject o : results) {
            if (o.getScore()>0.99) {
                double tmp = starts.get(o.getStart());
                tmp += o.getScore();
                starts.set(o.getStart(), tmp);
                
                tmp = ends.get(o.getEnd());
                tmp += o.getScore();
                ends.set(o.getEnd(), tmp);
                cpt++;
            }
        }
        if (debug) {
            System.out.println("total putative miRNAs=" + cpt);
        }
        
        if (debug) {
            System.out.println("Starts");
            for (int i = 0; i < starts.size(); i++) {
                System.out.println(i + "\t" + starts.get(i));                
            }
            System.out.println("\nEnds");
            for (int i = 0; i < ends.size(); i++) {
                System.out.println(i + "\t" + ends.get(i));                
            }
        }
        
        ArrayList<String> alscoresmirnas = new ArrayList<String>();
        ArrayList<Double> alscores = new ArrayList<Double>();
        for (int i = minlength; i <= maxlength; i++) {
            for (int j = 0; j < prec.length()-i; j++) {
                double ss=starts.get(j);
                double se=ends.get(j+i);
                String mirna=prec.substring(j, j+i);
                double score=ss+se;
                alscores.add(score);
                alscoresmirnas.add(score+"\t"+mirna);
            }            
        }
        
        double maxScore=0;
        for (Double d : alscores) {
            if (d>maxScore){
                maxScore=d;
            }
        }
        String consensus="";
        for (String s : alscoresmirnas) {
            if (s.contains(String.valueOf(maxScore))){
                System.out.println(s);
                consensus=s.split("\t")[1];
            }
        }
        int startPosition=prec.indexOf(consensus);
        int endPosition=startPosition+consensus.length();
        
        int startOfKnownMiRNA=prec.indexOf(knownMirna);
        int endOfKnownMiRNA=prec.indexOf(knownMirna)+knownMirna.length();

        int diffStartFromMiRNA=Math.abs(startPosition-startOfKnownMiRNA);
        int diffendFromMiRNA=Math.abs(endPosition-endOfKnownMiRNA);
        
        Features f = new Features(knownMirna, prec, structure, true);
        String miRNAStar=f.getMirnaStar();
        int startOfKnownMiRNAStar=prec.indexOf(miRNAStar);
        int endOfKnownMiRNAStar=prec.indexOf(miRNAStar)+miRNAStar.length();
        
        int diffStartFromMiRNAStar=Math.abs(startPosition-startOfKnownMiRNAStar);
        int diffendFromMiRNAStar=Math.abs(endPosition-endOfKnownMiRNAStar);
        
        System.out.println("Start difference from known miRNA: "+diffStartFromMiRNA);
        System.out.println("End difference from known miRNA: "+diffendFromMiRNA);
        System.out.println("Start difference from known miRNAstar: "+diffStartFromMiRNAStar);
        System.out.println("End difference from known miRNAstar: "+diffendFromMiRNAStar);
        
        if (diffStartFromMiRNAStar<diffStartFromMiRNA){
            System.out.println("predicted consensus miRNA is the miRNAStar: "+consensus);
            return diffStartFromMiRNAStar+","+diffendFromMiRNAStar;
        } else {
            System.out.println("predicted consensus miRNA is the miRNA: "+consensus);
            return diffStartFromMiRNA+","+diffendFromMiRNA;
        }
        
        
    }   
    
            /**
     * Decision by getting predictions consensus
     * @param results
     * @return 
     */
    private static String decisionByProducts(ArrayList<AlignmentObject> results) {
        ArrayList<String> PhrasesScores = new ArrayList<String>();   
        HashMap<String,AlignmentObject> hm = new HashMap<String, AlignmentObject>();
        
        for (AlignmentObject o : results) {
            if (o.getScore()>0.99) {
                PhrasesScores.add(o.getPhrase());
                hm.put(o.getPhrase(), o);
            }
        }
        // order starts and ends by score
        Collections.sort(PhrasesScores);
        Collections.reverse(PhrasesScores);
        
        // get best arm
        int cprim=0;
        int tprim=0;
        for (int i = 0; i < 20; i++) {
            if (hm.get(PhrasesScores.get(i)).getArm().equals("5'")){
                cprim++;
            } else {
                tprim++;
            }            
        }
        String bestArm="";
        if (cprim>tprim){
            bestArm="5'";
        } else {
            bestArm="3'";
        }
        
        if (PhrasesScores.size()==1){
            nbrproducts=1;        
        } else if (PhrasesScores.size()>=10){
            nbrproducts=10;
        } else if (PhrasesScores.size()>=5){
            nbrproducts=5;
        } else if (PhrasesScores.size()>=2){
            nbrproducts=2;
        }
        
        // get products of the alignment
        String consensus="";
        boolean miRNALengthOK=false;
        while(!miRNALengthOK){
            ArrayList<String> mirnasAln=new ArrayList<String>();
            int cpt=0;
            int maxLength=0;
            
            while (cpt!=nbrproducts){
                AlignmentObject ao=hm.get(PhrasesScores.get(cpt));
                if (ao.getArm().equals(bestArm)) {
                    mirnasAln.add(ao.getAlignment());
                    if (ao.getAlignment().length()>maxLength){
                        maxLength=ao.getAlignment().length();
                    }                      
                }   
                cpt++;
            }
            
            // getting consensus sequence of the alignment
            String bigConsensus="";
            for (int i = 0; i < maxLength; i++) {
                char c = ' ';
                for (String s : mirnasAln) {                    
                    if (c==' '){
                        try {
                            c = s.charAt(i);
                        } catch (Exception e) {
                        }
                    }
                    
                }  
                bigConsensus+=c;
            }

            //getting count in alignment colunms
            String alnCount="";
            for (int i = 0; i < maxLength; i++) {
                cpt=0;
                for (String s : mirnasAln) {
                    try {
                        if (s.charAt(i) != ' ') {
                            cpt++;
                        }
                    } catch (Exception e) {
                    }
                }    
                if (cpt>=10){
                    alnCount+="0";
                } else if (cpt>0){
                    alnCount+=cpt;
                } else {
                    alnCount+=".";
                }
            }
            
            // get positions
            ArrayList<Integer> al = new ArrayList<Integer>();
            for (int i = 0; i < alnCount.length(); i++) {
                while(alnCount.charAt(i)=='.'){
                    i++;
                }
                String s = String.valueOf(alnCount.charAt(i));
                int v=Integer.valueOf(s);
                
                if (v==0){
                    v=10;
                }
                if (v>=threshold){
                    al.add(i);                    
                }
            }
            
            int start=al.get(0);
            int end=al.get(al.size()-1);
            
            consensus=bigConsensus.substring(start,end);
            if (consensus.length()<16){
                threshold--;
                miRNALengthOK=false;
            } else {
                miRNALengthOK=true;
            }
        }
        
        String konwnMirna = results.get(0).getKnownmirna();
        String prec=results.get(0).getPrec();
        int startPosition=prec.indexOf(consensus);
        int endPosition=startPosition+consensus.length();
        
        int startOfKnownMiRNA=prec.indexOf(konwnMirna);
        int endOfKnownMiRNA=prec.indexOf(konwnMirna)+konwnMirna.length();

        int diffStart=Math.abs(startPosition-startOfKnownMiRNA);
        int diffend=Math.abs(endPosition-endOfKnownMiRNA);

        System.out.println("Start difference from known miRNA: "+diffStart);
        System.out.println("End difference from known miRNA: "+diffend);

        System.out.println("predicted consensus miRNA: "+consensus);
        return diffStart+","+diffend;
    } 
}