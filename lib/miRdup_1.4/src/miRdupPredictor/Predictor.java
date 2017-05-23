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
package miRdupPredictor;

import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import miRdup.Features;
import miRdup.Main;
import miRdup.Vienna;

/**
 *
 * @author fyrox
 */
public class Predictor {

    //public static String rnafold="/home/mycky/tools/ViennaRNA-2.0.6/Progs/";
    //public static String rnafold="/ibrixfs1/Data/mik/tools/ViennaRNA-2.0.7/Progs/";
    public static String rnafold = "";
    public static int loopflag = 0; //flag to detect looping over the same function by recursion
    public static int maximas = 0; //When looping over the same function by recursion, this value take the second (three, and so on) maxima instead of the first
    public static boolean debug = false;
    public static int threshold = 10; //decisionbyproducts
    public static int nbrproducts = 20; //decisionbyproducts

    //generate miRNAs length
    public static int minlength = 16;
    public static int maxlength = 30;

    //choose miRNAs by length
    public static int minlengthDecision = 16;
    public static int maxlengthDecision = 30;

    public static String struct = "";

    public static String predictionBySequence(String precursor, String model, String outfile) {
        System.out.println("==========\nPredict miRNA positions");
        System.out.println(precursor);
        String result = executeGenerator(precursor, model, outfile + ".generatedmirnas");
        System.out.println("Predicted miRNA in 3' and 5':" + result);
        return result;
    }

    public static void predictionByFile(String infile, String model, String outfile) {

        System.out.println("==========\nPredict miRNA positions");
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            PrintWriter pw = new PrintWriter(new FileWriter(infile + ".miRdup.predictions.txt"));
            String line = "";
            while (br.ready()) {
                try {
                    line = br.readLine();
                    loopflag = 0;
                    maximas = 0;

                    struct = "";
                    System.out.println("\n----------");

                    String name = "";
                    String prec = "";
                    if (line.split("\t").length == 2) {
                        name = line.split("\t")[0];
                        prec = line.split("\t")[1].replace("T", "U");
                        System.out.println(line + "\t");
                    } else if (line.startsWith(">")) {
                        name = line.substring(1);
                        line = br.readLine();
                        prec = line.replace("T", "U");;
                        System.out.println(name + "\t" + prec);
                    }
                    minlengthDecision = 16;
                    maxlengthDecision = 30;
                    String result = executeGenerator(prec, model, outfile + "_tmpfile");
                    System.out.println("Predicted miRNA in 3' and 5': " + result);
                    pw.println(name + "\t" + prec + "\t" + struct + "\t" + result);
                    pw.flush();
                } catch (Exception e) {
                    System.err.println("error at " + line);
                    e.printStackTrace();
                }
            }
            pw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static String executeGenerator(String prec, String model, String outfile) {
        //generate miRNAs
        generateMiRNAs(prec, outfile);

        //Predict miRNAs with miRdup
        predict(outfile, model);

        //Align true miRNAs on the precursor
        File modelFile = new File(model);
        String predictionsFile = outfile + "." + modelFile.getName() + ".miRdup.txt"; // file produced by mirdup
        String alnFile = outfile + "." + modelFile.getName() + ".miRdup.aln.txt"; // output
        ArrayList<AlignmentObject> results = makeAlignment(predictionsFile, alnFile);
//        return decisionByMean(results);
        try {
            //return decisionByMaximas(results);
            return decisionByWindow(results);
//            return decisionByProducts(results);
        } catch (Exception e) {
            e.printStackTrace();
            return "no predictable miRNA";
        }
    }

    /**
     * Generate all possible miRNAs, from length 16 to 30
     *
     * @param prec
     * @param outfile
     */
    public static void generateMiRNAs(String prec, String outfile) {
        ArrayList<String> mirnas = new ArrayList<String>();
        String os = System.getProperty("os.name").toLowerCase();
        if (!os.startsWith("win")) {
            Main.rnafoldlinux = rnafold;
        }
        struct = Vienna.GetSecondaryStructure(prec);
        System.out.println(struct);

        try {
            PrintWriter pw = new PrintWriter(new FileWriter(outfile));
            int cpt = 0;
            for (int i = minlength; i <= maxlength; i++) {
                for (int j = 0; j <= prec.length() - i; j++) {
                    cpt++;
                    String mirna = prec.substring(j, j + i);
                    //System.out.println(mirna);
                    pw.println("mir" + cpt + "\t" + mirna + "\t" + prec + "\t" + struct);
                    mirnas.add(mirna);
                }
            }
            pw.flush();
            pw.close();

            System.out.println("Total generated miRNAs: " + cpt);

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * launch miRdup
     *
     * @param outfile
     * @param model
     */
    private static void predict(String outfile, String model) {
        String s[] = {"-c", model, "-v", outfile, "-r", rnafold};
        setOptions(s);
        miRdupExecutionEMBL();
    }

    private static void setOptions(String[] args) {
        Main.setOptions(args);
        Main.predictMirnaPosition = false;
    }

    private static void miRdupExecutionEMBL() {
        Main.miRdupExecutionEMBL();
    }

    /**
     * Make alignment of predicted miRNAs on the pre-miRNA
     *
     * @param infile
     */
    private static ArrayList makeAlignment(String infile, String alnFile) {
        ArrayList<AlignmentObject> alObj = new ArrayList<AlignmentObject>();
        struct = "";
        String prec = null;
        int tot = 0;
        int trues = 0;
        // read prediction file output from miRdup
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            String line = "";
            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith("#ID")) {
                    tot++;
                    AlignmentObject ao = new AlignmentObject();
                    ao.setId(line);
                    line = br.readLine();

                    while (line.startsWith("#FF")) {
                        line = br.readLine();
                    }
                    if (line.startsWith("#WE")) {
                        line = br.readLine();
                    }
                    if (line.split("\t")[2].startsWith("t")) {
                        trues++;
                        ao.setScore(Double.parseDouble(br.readLine().split("\t")[2]));
                        prec = br.readLine().split("\t")[1];
                        ao.setPrec(prec);
                        struct = br.readLine().split("\t")[1];
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
        double maxScore = 0.0;
        for (AlignmentObject o : alObj) {
            double score = o.getScore();
            if (score > maxScore) {
                maxScore = score;
            }
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
            pw.flush();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return alObj;
    }

    /**
     * Decision by getting maximum values
     *
     * @param results
     * @return
     */
    private static String decisionByMaximas(ArrayList<AlignmentObject> results) {
        String prec = results.get(0).getPrec();

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

        int cpt = 0;
        for (AlignmentObject o : results) {
            if (o.getScore() > 0.99) {
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

        // order starts and ends by score
        Collections.sort(endsSort);
        Collections.reverse(endsSort);
        Collections.sort(startsSort);
        Collections.reverse(startsSort);

        //get best start and end
        double beststart = startsSort.get(0 + maximas);
        double bestend = endsSort.get(0 + maximas);

        //find corresponding positions of bestscores
        int startPosition = 0;
        int endPosition = 0;
        for (int i = 0; i < starts.size(); i++) {
            if (starts.get(i) == beststart) {
                startPosition = i;
            }
            if (ends.get(i) == bestend) {
                endPosition = i;
            }
        }

        // if start>end, we remove results from one of the arm of the hairpin
        // structure goes 5' to 3'
        if (startPosition > endPosition) {
            if (debug) {
                System.out.println("Selected miRNA have its start position>to end position... "
                        + "\nRemoving results from the arm having the lowest scores");
            }
            // count 5' and 3'
            int fivePrim = 0;
            int threePrim = 0;
            for (AlignmentObject o : results) {
                if (o.getArm().equals("3'")) {
                    threePrim++;
                }
                if (o.getArm().equals("5'")) {
                    fivePrim++;
                }
            }

            //remove objects from the arm having the lower count
            ArrayList<AlignmentObject> results2 = new ArrayList<AlignmentObject>();
            if (fivePrim >= threePrim) {
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("5'")) {
                        results2.add(o);
                    }
                }
            } else {
                for (AlignmentObject o : results) {
                    if (o.getArm().equals("3'")) {
                        results2.add(o);
                    }
                }
            }

            //re-execute the decision
            loopflag++;
            if (loopflag > 1) { //if we already pass through here, put maximas++
                maximas++;
            }

            return decisionByMaximas(results2);

        } else {
            String consensus = prec.substring(startPosition, endPosition);

            Features f = new Features(consensus, prec, struct, true);
            String complement = f.getComplementaritySequence();
            String miRNAstar = f.getMirnaStar();

            // verify miRNA length
            if (consensus.length() < 16 || consensus.length() > 33) {
                if (debug) {
                    System.out.println("Predicted miRNA too long, trying another maxima threshold...");
                }
                maximas++;
                loopflag++;
                if (loopflag > 10) {
                    System.out.println("Prediction impossible, return best scored predicted miRNA");
                    return results.get(0).getMirna();
                } else {
                    return decisionByMaximas(results);
                }
            } else {
                if (debug) {
                    System.out.println("predicted consensus miRNA: " + consensus);
                    //System.out.println("predicted consensus miRNA complement: "+complement);
                    System.out.println("predicted consensus miRNA star: " + miRNAstar);
                }
                return consensus;
            }
        }

    }

    /**
     * Decision by sliding window
     *
     * @param results
     * @return
     */
    private static String decisionByWindow(ArrayList<AlignmentObject> results) {
        String prec = results.get(0).getPrec();
        String structure = results.get(0).getStruct();
        ArrayList<Double> starts = new ArrayList<Double>();
        ArrayList<Double> ends = new ArrayList<Double>();

        // full lists with 0
        for (int i = 0; i <= prec.length(); i++) {
            starts.add(0.0);
            ends.add(0.0);
        }

        int cpt = 0;
        for (AlignmentObject o : results) {
            if (o.getScore() > 0.99) {
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
        for (int i = minlengthDecision; i <= maxlengthDecision; i++) {
            for (int j = 0; j < prec.length() - i; j++) {
                double ss = starts.get(j);
                double se = ends.get(j + i);
                String mirna = prec.substring(j, j + i/*+1*/);
                double score = ss + se;
                alscores.add(score);
                alscoresmirnas.add(score + "\t" + mirna);
            }
        }

        double maxScore = 0;
        for (Double d : alscores) {
            if (d > maxScore) {
                maxScore = d;
            }
        }
        String consensus = "";
        for (String s : alscoresmirnas) {
            if (s.contains(String.valueOf(maxScore))) {
                //System.out.println(s);
                consensus = s.split("\t")[1];
            }
        }
        Features f = new Features(consensus, prec, structure, true);
        String miRNAStar = f.getMirnaStar();
        String consensusArm = f.getArm();
        String starArm = "";
        if (consensusArm.equals("3'")) {
            starArm = "5'";
        } else {
            starArm = "3'";
        }
        if (miRNAStar.length() > 30) {
            minlengthDecision = minlengthDecision - 1;
            maxlengthDecision = maxlengthDecision - 1;
            return decisionByWindow(results);
        } else {
            return consensus + "(" + consensusArm + ")\t" + miRNAStar + "(" + starArm + ")";
        }
    }

    /**
     * Decision by getting predictions consensus
     *
     * @param results
     * @return
     */
    private static String decisionByProducts(ArrayList<AlignmentObject> results) {
        ArrayList<String> PhrasesScores = new ArrayList<String>();
        HashMap<String, AlignmentObject> hm = new HashMap<String, AlignmentObject>();

        for (AlignmentObject o : results) {
            if (o.getScore() > 0.99) {
                PhrasesScores.add(o.getPhrase());
                hm.put(o.getPhrase(), o);
            }
        }
        // order starts and ends by score
        Collections.sort(PhrasesScores);
        Collections.reverse(PhrasesScores);

        // get nbr products
        if (PhrasesScores.size() == 1) {
            nbrproducts = 1;
        } else if (PhrasesScores.size() >= 15) {
            nbrproducts = 15;
        } else if (PhrasesScores.size() >= 10) {
            nbrproducts = 10;
        } else if (PhrasesScores.size() >= 5) {
            nbrproducts = 5;
        } else if (PhrasesScores.size() >= 2) {
            nbrproducts = 2;
        }

        // get best arm
        int cprim = 0;
        int tprim = 0;
        for (int i = 0; i < nbrproducts; i++) {
            if (hm.get(PhrasesScores.get(i)).getArm().equals("5'")) {
                cprim++;
            } else {
                tprim++;
            }
        }
        String bestArm = "";
        if (cprim > tprim) {
            bestArm = "5'";
        } else {
            bestArm = "3'";
        }

        // get products of the alignment
        String consensus = "";
        boolean miRNALengthOK = false;
        while (!miRNALengthOK) {
            ArrayList<String> mirnasAln = new ArrayList<String>();
            int cpt = 0;
            int maxLength = 0;

            while (cpt != nbrproducts) {
                AlignmentObject ao = hm.get(PhrasesScores.get(cpt));
                if (ao.getArm().equals(bestArm)) {
                    mirnasAln.add(ao.getAlignment());
                    if (ao.getAlignment().length() > maxLength) {
                        maxLength = ao.getAlignment().length();
                    }
                }
                cpt++;
            }

            // getting consensus sequence of the alignment
            String bigConsensus = "";
            for (int i = 0; i < maxLength; i++) {
                char c = ' ';
                for (String s : mirnasAln) {
                    if (c == ' ') {
                        try {
                            c = s.charAt(i);
                        } catch (Exception e) {
                        }
                    }

                }
                bigConsensus += c;
            }

            //getting count in alignment colunms
            String alnCount = "";
            for (int i = 0; i < maxLength; i++) {
                cpt = 0;
                for (String s : mirnasAln) {
                    try {
                        if (s.charAt(i) != ' ') {
                            cpt++;
                        }
                    } catch (Exception e) {
                    }
                }
                if (cpt >= 10) {
                    alnCount += "0";
                } else if (cpt > 0) {
                    alnCount += cpt;
                } else {
                    alnCount += ".";
                }
            }

            // get positions
            ArrayList<Integer> al = new ArrayList<Integer>();
            for (int i = 0; i < alnCount.length(); i++) {
                while (alnCount.charAt(i) == '.') {
                    i++;
                }
                String s = String.valueOf(alnCount.charAt(i));
                int v = Integer.valueOf(s);

                if (v == 0) {
                    v = 10;
                }
                if (v >= threshold) {
                    al.add(i);
                }
            }

            int start = al.get(0);
            int end = al.get(al.size() - 1);

            consensus = bigConsensus.substring(start, end);
            if (consensus.length() < 16) {
                threshold--;
                miRNALengthOK = false;
            } else {
                miRNALengthOK = true;
            }
        }

        return consensus;
    }

}
