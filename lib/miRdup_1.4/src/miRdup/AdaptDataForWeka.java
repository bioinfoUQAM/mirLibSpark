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
 * Adapt data for weka
 */
package miRdup;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author Mickael Leclercq
 */
public class AdaptDataForWeka {

    public static int id = 0;

    /**
     * Create a file from an arraylist of mirnas and their precursors from
     * mirbase
     *
     * @param al
     * @param outfile
     * @param bestfeatures
     */
    static void createFileFromList(ArrayList al, File outfile, boolean bestfeatures) {
        try {
            PrintWriter pw;
            pw = new PrintWriter(new FileWriter(outfile));
            if (bestfeatures) {
                pw.println(getBestWekaFileEntries());
            } else {
                pw.println(getAllWekaFileEntries());
            }
            pw.flush();
            System.out.println("Get test positive and negative dataset features");
            ArrayList data = generateDataFromObjectList(al, bestfeatures);
            for (int i = 0; i < data.size(); i++) {
                pw.println(data.get(i));

            }

            pw.flush();
            pw.println();
            pw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param al
     * @param positiveDataset
     * @return
     */
    private static ArrayList<String> generateDataFromObjectList(ArrayList<MirnaObject> alobj, boolean bestfeatures) {
        ArrayList<String> al = new ArrayList<String>();
        int cpt = 0;
        for (MirnaObject o : alobj) {
            try {
                String mirna = o.getMatureSequence();
                String prec = o.getPrecursorSequence();
                String struc = o.getStructure();
                Features f;
                try {
                    f = new Features(mirna, prec, struc, true);
                    if (bestfeatures) {
                        al.add(f.toStringBestAttributes());
                    } else {
                        al.add(f.toStringAllAttributes());
                    }

                    mirna = mirnaNegativeData(mirna, prec, struc, alobj); //generate new mirna
                    f = new Features(mirna, prec, struc, false);

                    if (bestfeatures) {
                        al.add(f.toStringBestAttributes());
                    } else {
                        al.add(f.toStringAllAttributes());
                    }
                } catch (Exception e) {
                    System.err.println("Bad hairpin at " + o.toString());
                    f = new Features();
                    al.add(f.toStringError());
                }

                cpt++;
                if (cpt % 1000 == 0) {
                    System.out.print("*");
                }
            } catch (Exception e) {
                e.printStackTrace();

            }
        }
        return al;

    }

    /**
     * dataset must be:
     *
     * @param predictionSetTabbed
     */
    public static void createPredictionDataset(String predictionSetTabbed, boolean bestAttributes) {

        try {
//            Predicted set
            PrintWriter pw = new PrintWriter(new FileWriter(predictionSetTabbed + ".arff"));
            if (bestAttributes) {
                pw.println(getBestWekaFileEntries());
            } else {
                pw.println(getAllWekaFileEntries());
            }
            pw.flush();
            System.out.println("Storing predicted dataset features in " + predictionSetTabbed + ".arff");
            File f = new File(predictionSetTabbed);
            for (String s : generateData(f, true, bestAttributes)) {
                pw.println(s);
            }
            pw.println();
            pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Generate DataSet, positive or negative, from an infile. The boolean value
     * must set as true or false depending if we want to generate positive or
     * negative dataset respectively
     *
     * @param infile
     * @param positiveDataset
     * @param bestAttributes
     * @return Arraylist of features for each iteration in Weka format
     * (separated by ,) positive dataset must be: mirna"\t"prec"\t"struc
     * negative dataset must be: prec"\t"struc ; mirna will be randomly
     * generated
     *
     */
    public static ArrayList<String> generateData(File infile, Boolean positiveDataset, Boolean bestAttributes) {
        ArrayList<String> al = new ArrayList<String>();
        int cpt = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            File fold = new File(infile + ".folded");
            PrintWriter pw = new PrintWriter(new FileWriter(fold));
            boolean folded = false;
            String line = "";
            while (br.ready()) {
                line = br.readLine();
                String mirna = "";
                try {
                    String tab[] = line.split("\t");
                    String idName = tab[0];
                    String mirna_or = tab[1];
                    mirna = tab[1].replaceAll("-", "");
                    String prec = tab[2];//.replaceAll("[^ACGTU]", "");
                    String struc = null;
                    try {
                        struc = tab[3];
                    } catch (Exception e) {
                        struc = Vienna.GetSecondaryStructure(prec);
                        pw.println(idName + "\t" + mirna_or + "\t" + prec + "\t" + struc);
                        folded = true;
                    }
                    Features f = null;
                    try {
                        if (mirna.length() >= 16 && mirna.length() <= 31) {
                            if (positiveDataset) {
                                f = new Features(mirna, prec, struc, true);
                            } else {
                                mirna = mirnaNegativeData(mirna, prec, struc, new ArrayList<MirnaObject>()); //generate new mirna
                                f = new Features(mirna, prec, struc, false);
                            }
                            if (bestAttributes) {
                                al.add(f.toStringBestAttributes());
                            } else {
                                al.add(f.toStringAllAttributes());
                            }

                        } else {
                            System.err.println("Highly improbable miRNA (lenght not between 16 and 31 nt) " + line);
                            f = new Features();
                            al.add(f.toStringError());
                        }

                    } catch (Exception e) {
                        System.err.println("Bad hairpin at " + line);
                        f = new Features();
                        al.add(f.toStringError());
                    }

                    cpt++;
                    if (cpt % 1000 == 0) {
                        System.out.print("*");
                        pw.flush();
                    }
                } catch (Exception e) {

                }
            }
            pw.flush();
            pw.close();
            if (!folded) {
                fold.delete();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }

        return al;
    }

    /**
     * generate miRNA negative Data in moving the miRNA in another place on the
     * precursor
     *
     * @param mirna
     * @param precursor
     * @param structure
     * @return
     */
    public static String mirnaNegativeDatabak(String mirna, String precursor, String structure,
            ArrayList<MirnaObject> alobj) {
        String newmirna = "";
        Features f = new Features(mirna, precursor, structure, false);
        newmirna = f.getMirnaStar();
        return newmirna;
    }

    public static String mirnaNegativeData(String mirna, String precursor, String structure,
            ArrayList<MirnaObject> alobj) {

        boolean newmirnaIsValid = false;
        String newmirna = "";
        while (newmirnaIsValid == false) {
            Random r = new Random();
            int mirnalength = mirna.length();
            int start = r.nextInt(precursor.length() - mirnalength);
            newmirna = precursor.substring(start, start + mirnalength);
            if (!mirna.equals(newmirna)) { // in case of bad random where mirna=new mirna
                for (MirnaObject o : alobj) { // verify in all set if we don't encounter the mirna
                    if (o.getMatureSequence().equals(newmirna)) {
                        if (precursor.contains(o.getPrecursorSequence()) || o.getPrecursorSequence().contains(precursor)) {
                            newmirnaIsValid = false;
                        } else {
                            newmirnaIsValid = true;
                        }
                    } else {
                        newmirnaIsValid = true;
                    }
                }
            }
        }

        return newmirna;
    }

    /**
     * Print all features
     *
     * @param f
     */
    public static void printFeatures(Features f) {
        p(f.getLength());
        p(f.getGCperc());
        p(f.getMaximumLengthWithoutBulges());
        p(f.getMaximumLengthWithoutBulgesPerc());
        p(f.getStartLengthWithoutBulges());
        p(f.getBasePairsInDuplex());
        p(f.getPresenceOfPerfect20MerBasePair());
        p(f.getStartOfPerfect20MerBasePair());
        p(f.getPresenceOfPerfect10MerBasePair());
        p(f.getStartOfPerfect10MerBasePair());
        p(f.getPresenceOfPerfect5MerBasePair());
        p(f.getStartOfPerfect5MerBasePair());
        p(f.getPresenceOfA());
        p(f.getPresenceOfU());
        p(f.getPresenceOfG());
        p(f.getPresenceOfC());
        p(f.getDistanceFromTerminalLoop());
        p(f.getDistanceFromHairpinStart());
        p(f.getMirnaIncludedInLoop());
        p(f.getLengthOfOverlapInLoop());
        p(f.getAverageNumberOfPairedBasesInWindow7());
        p(f.getAverageNumberOfPairedBasesInWindow5());
        p(f.getAverageNumberOfPairedBasesInWindow3());
        p(f.getBulgeAtPosition1());
        p(f.getBulgeAtPositionMinus2());
        p(f.getBulgeAtPosition0());
        p(f.getBulgeAtPositionMinus1());
        p(f.getNumberOfBulges());
        p(f.getLenghtOfBiggestBulge());
        p(f.getLengthBiggestBulgesPerc());

        System.out.println("");
    }

    public static void p(String s) {
        System.out.println(s);
    }

    public static void p(int s) {
        System.out.println(s);
    }

    public static void p(double s) {
        //System.out.println(df.format(s));
    }

    public static void p(boolean b) {
        System.out.println(b);
    }

    /**
     * Return all weka entries for the beginning of arff file
     *
     * @return
     */
    public static String getAllWekaFileEntries() {
        String s = ""
                + "@relation mirnaInPrecursor" + "\n"
                + "@attribute id string" + "\n"
                + "@attribute length real" + "\n"
                + "@attribute mfe real" + "\n" //RNAcofold calculate it
                + "@attribute GCperc real" + "\n"
                //                + "@attribute GCpercNormalized real"+ "\n"
                + "@attribute MaximumLengthWithoutBulges real" + "\n"
                + "@attribute MaximumLengthWithoutBulgesPerc real" + "\n"
                + "@attribute StartLengthWithoutBulges real" + "\n"
                + "@attribute BasePairsInDuplexMirnaMirnaStar real" + "\n"
                + "@attribute PresenceOfPerfect20MerBasePair { true, false}" + "\n"
                + "@attribute StartOfPerfect20MerBasePair real" + "\n"
                + "@attribute PresenceOfPerfect10MerBasePair { true, false}" + "\n"
                + "@attribute StartOfPerfect10MerBasePair real" + "\n"
                + "@attribute PresenceOfPerfect5MerBasePair { true, false}" + "\n"
                + "@attribute StartOfPerfect5MerBasePair real" + "\n"
                + "@attribute PercOfA real" + "\n"
                + "@attribute PercOfU real" + "\n"
                + "@attribute PercOfG real" + "\n"
                + "@attribute PercOfC real" + "\n"
                + "@attribute PercOfAA real" + "\n"
                + "@attribute PercOfUA real" + "\n"
                + "@attribute PercOfGA real" + "\n"
                + "@attribute PercOfCA real" + "\n"
                + "@attribute PercOfAU real" + "\n"
                + "@attribute PercOfUU real" + "\n"
                + "@attribute PercOfGU real" + "\n"
                + "@attribute PercOfCU real" + "\n"
                + "@attribute PercOfAG real" + "\n"
                + "@attribute PercOfUG real" + "\n"
                + "@attribute PercOfGG real" + "\n"
                + "@attribute PercOfCG real" + "\n"
                + "@attribute PercOfAC real" + "\n"
                + "@attribute PercOfUC real" + "\n"
                + "@attribute PercOfGC real" + "\n"
                + "@attribute PercOfCC real" + "\n"
                + "@attribute DistanceFromTerminalLoop real" + "\n"
                + "@attribute DistanceFromHairpinStart real" + "\n"
                + "@attribute MirnaIncludedInLoop { true, false}" + "\n"
                + "@attribute LengthOfOverlapInLoop real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow7 real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow5 real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow3 real" + "\n"
                + "@attribute BulgeAtPosition0 { true, false}" + "\n"
                + "@attribute BulgeAtPositionMinus1 { true, false}" + "\n"
                + "@attribute BulgeAtPositionPlus1 { true, false}" + "\n"
                + "@attribute BulgeAtPositionMinus2 { true, false}" + "\n"
                + "@attribute BulgeAtPositionPlus2 { true, false}" + "\n"
                + "@attribute BulgeAtPositionMinus3 { true, false}" + "\n"
                + "@attribute BulgeAtPositionPlus3 { true, false}" + "\n"
                + "@attribute BulgeAtPositionMinus4 { true, false}" + "\n"
                + "@attribute BulgeAtEndPosition0 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionPlus1 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionMinus1 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionPlus2 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionMinus2 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionPlus3 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionMinus3 { true, false}" + "\n"
                + "@attribute BulgeAtEndPositionPlus4 { true, false}" + "\n"
                + "@attribute NumberOfBulges real" + "\n"
                + "@attribute LenghtOfBiggestBulge real" + "\n"
                + "@attribute LengthBiggestBulgesPerc real" + "\n"
                + "@attribute A... real" + "\n"
                + "@attribute C... real" + "\n"
                + "@attribute G... real" + "\n"
                + "@attribute U... real" + "\n"
                + "@attribute A(.. real" + "\n"
                + "@attribute C(.. real" + "\n"
                + "@attribute G(.. real" + "\n"
                + "@attribute U(.. real" + "\n"
                + "@attribute A((. real" + "\n"
                + "@attribute C((. real" + "\n"
                + "@attribute G((. real" + "\n"
                + "@attribute U((. real" + "\n"
                + "@attribute A.(( real" + "\n"
                + "@attribute C.(( real" + "\n"
                + "@attribute G.(( real" + "\n"
                + "@attribute U.(( real" + "\n"
                + "@attribute A((( real" + "\n"
                + "@attribute C((( real" + "\n"
                + "@attribute G((( real" + "\n"
                + "@attribute U((( real" + "\n"
                + "@attribute A.(. real" + "\n"
                + "@attribute C.(. real" + "\n"
                + "@attribute G.(. real" + "\n"
                + "@attribute U.(. real" + "\n"
                + "@attribute A..( real" + "\n"
                + "@attribute C..( real" + "\n"
                + "@attribute G..( real" + "\n"
                + "@attribute U..( real" + "\n"
                + "@attribute A(.( real" + "\n"
                + "@attribute C(.( real" + "\n"
                + "@attribute G(.( real" + "\n"
                + "@attribute U(.( real" + "\n"
                + "@attribute PercBasedPairAU real" + "\n"
                + "@attribute PercBasedPairGC real" + "\n"
                + "@attribute PercBasedPairGU real" + "\n"
                //                + "@attribute PercBasedPairAUnorm real"+ "\n"
                //                + "@attribute PercBasedPairGCnorm real"+ "\n"
                //                + "@attribute PercBasedPairGUnorm real"+ "\n"

                + "@attribute nucleotideAt0 { A, U, G, C, N}" + "\n"
                + "@attribute nucleotideAtMinus1 { A, U, G, C, N}" + "\n"
                //                + "@attribute nucleotideAtMinus2 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtMinus3 { A, U, G, C, N}"+ "\n"
                + "@attribute nucleotideAtPlus1 { A, U, G, C, N}" + "\n"
                //                + "@attribute nucleotideAtPlus2 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtPlus3 { A, U, G, C, N}"+ "\n"
                + "@attribute nucleotideAtEnd0 { A, U, G, C, N}" + "\n"
                + "@attribute nucleotideAtEndMinus1 { A, U, G, C, N}" + "\n"
                //                + "@attribute nucleotideAtEndMinus2 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtEndMinus3 { A, U, G, C, N}"+ "\n"
                + "@attribute nucleotideAtEndPlus1 { A, U, G, C, N}" + "\n"
                //                + "@attribute nucleotideAtEndPlus2 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtEndPlus3 { A, U, G, C, N}"+ "\n"

                //                + "@attribute RNAcofoldMfe real"+"\n"
                //                + "@attribute RNAcofoldFrequency real"+"\n"
                //                + "@attribute RNAcofoldDeltaG real"+"\n"
                //                + "@attribute RNAcofoldAB real"+"\n"
                //                + "@attribute RNAcofoldAA real"+"\n"
                //                + "@attribute RNAcofoldBB real"+"\n"
                //                + "@attribute RNAcofoldA real"+"\n"
                //                + "@attribute RNAcofoldB real"+"\n"
                //                + "@attribute BPprobGlobal real"+"\n"
                //                + "@attribute BPprobTmpStructure real"+"\n"
                //                + "@attribute BPprobFinalStructure real"+"\n"
                //                + "@attribute BPprobEnsemblizedBulges real"+"\n"
                //                + "@attribute ComplexBoltzmannProbability real"+"\n"

                //                + "@attribute BPprobFinalStructureGC real"+"\n"
                //                + "@attribute BPprobFinalStructureGU real"+"\n"
                //                + "@attribute BPprobFinalStructureAU real"+"\n"
                //                + "@attribute BPprobFinalStructureGG real"+"\n"
                //                + "@attribute BPprobFinalStructureCC real"+"\n"
                //                + "@attribute BPprobFinalStructureCU real"+"\n"
                //                + "@attribute BPprobFinalStructureCA real"+"\n"
                //                + "@attribute BPprobFinalStructureGA real"+"\n"
                //                + "@attribute BPprobFinalStructureAA real"+"\n"
                //                + "@attribute BPprobFinalStructureUU real"+"\n"

                + "@attribute Class { true, false}" + "\n"
                + "@data";
        return s;
    }

    /**
     * Return all weka entries for the beginning of arff file
     *
     * @return
     */
    public static String getBestWekaFileEntries() {
        String s = ""
                + "@relation mirnaInPrecursor" + "\n"
                + "@attribute id string" + "\n"
                //                + "@attribute length real"+ "\n"
                //                + "@attribute mfe real"+ "\n"
                //                + "@attribute GCperc real"+ "\n"
                + "@attribute MaximumLengthWithoutBulges real" + "\n"
                + "@attribute MaximumLengthWithoutBulgesPerc real" + "\n"
                //                + "@attribute StartLengthWithoutBulges real"+ "\n"
                + "@attribute BasePairsInDuplexMirnaMirnaStar real" + "\n"
                //                + "@attribute PresenceOfPerfect20MerBasePair { true, false}"+ "\n"
                //                + "@attribute StartOfPerfect20MerBasePair real"+ "\n"
                //                + "@attribute PresenceOfPerfect10MerBasePair { true, false}"+ "\n"
                //                + "@attribute StartOfPerfect10MerBasePair real"+ "\n"
                //                + "@attribute PresenceOfPerfect5MerBasePair { true, false}"+ "\n"
                + "@attribute StartOfPerfect5MerBasePair real" + "\n"
                //
                //                + "@attribute PercOfA real"+ "\n"
                //                + "@attribute PercOfU real"+ "\n"
                //                + "@attribute PercOfG real"+ "\n"
                //                + "@attribute PercOfC real"+ "\n"
                //                + "@attribute PercOfAA real"+ "\n"
                //                + "@attribute PercOfUA real"+ "\n"
                //                + "@attribute PercOfGA real"+ "\n"
                //                + "@attribute PercOfCA real"+ "\n"
                //                + "@attribute PercOfAU real"+ "\n"
                //                + "@attribute PercOfUU real"+ "\n"
                //                + "@attribute PercOfGU real"+ "\n"
                //                + "@attribute PercOfCU real"+ "\n"
                //                + "@attribute PercOfAG real"+ "\n"
                //                + "@attribute PercOfUG real"+ "\n"
                //                + "@attribute PercOfGG real"+ "\n"
                //                + "@attribute PercOfCG real"+ "\n"
                //                + "@attribute PercOfAC real"+ "\n"
                //                + "@attribute PercOfUC real"+ "\n"
                //                + "@attribute PercOfGC real"+ "\n"
                //                + "@attribute PercOfCC real"+ "\n"
                //
                + "@attribute DistanceFromTerminalLoop real" + "\n"
                + "@attribute DistanceFromHairpinStart real" + "\n"
                + "@attribute MirnaIncludedInLoop { true, false}" + "\n"
                + "@attribute LengthOfOverlapInLoop real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow7 real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow5 real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow3 real" + "\n"
                //
                //                + "@attribute BulgeAtPosition0 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionMinus1 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionPlus1 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionMinus2 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionPlus2 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionMinus3 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionPlus3 { true, false}"+ "\n"
                //                + "@attribute BulgeAtPositionMinus4 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPosition0 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionPlus1 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionMinus1 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionPlus2 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionMinus2 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionPlus3 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionMinus3 { true, false}"+ "\n"
                //                + "@attribute BulgeAtEndPositionPlus4 { true, false}"+ "\n"
                //
                //                + "@attribute NumberOfBulges real"+ "\n"
                + "@attribute LenghtOfBiggestBulge real" + "\n"
                + "@attribute LengthBiggestBulgesPerc real" + "\n"
                + "@attribute A... real" + "\n"
                + "@attribute C... real" + "\n"
                + "@attribute G... real" + "\n"
                + "@attribute U... real" + "\n"
                //                + "@attribute A(.. real"+ "\n"
                //                + "@attribute C(.. real"+ "\n"
                //                + "@attribute G(.. real"+ "\n"
                //                + "@attribute U(.. real"+ "\n"
                //                + "@attribute A((. real"+ "\n"
                //                + "@attribute C((. real"+ "\n"
                //                + "@attribute G((. real"+ "\n"
                //                + "@attribute U((. real"+ "\n"
                //                + "@attribute A.(( real"+ "\n"
                //                + "@attribute C.(( real"+ "\n"
                //                + "@attribute G.(( real"+ "\n"
                //                + "@attribute U.(( real"+ "\n"
                //                + "@attribute A((( real"+ "\n"
                //                + "@attribute C((( real"+ "\n"
                + "@attribute G((( real" + "\n"
                //                + "@attribute U((( real"+ "\n"
                //                + "@attribute A.(. real"+ "\n"
                //                + "@attribute C.(. real"+ "\n"
                //                + "@attribute G.(. real"+ "\n"
                //                + "@attribute U.(. real"+ "\n"
                //                + "@attribute A..( real"+ "\n"
                //                + "@attribute C..( real"+ "\n"
                //                + "@attribute G..( real"+ "\n"
                //                + "@attribute U..( real"+ "\n"
                //                + "@attribute A(.( real"+ "\n"
                //                + "@attribute C(.( real"+ "\n"
                //                + "@attribute G(.( real"+ "\n"
                //                + "@attribute U(.( real"+ "\n"
                + "@attribute PercBasedPairAU real" + "\n"
                + "@attribute PercBasedPairGC real" + "\n"
                + "@attribute PercBasedPairGU real" + "\n"
                //
                //                + "@attribute nucleotideAt0 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtMinus1 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtMinus2 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtMinus3 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtPlus1 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtPlus2 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtPlus3 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtEnd0 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtEndMinus1 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtEndMinus2 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtEndMinus3 { A, U, G, C, N}"+ "\n"
                //                + "@attribute nucleotideAtEndPlus1 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtEndPlus2 { A, U, G, C, N}"+ "\n"
                ////                + "@attribute nucleotideAtEndPlus3 { A, U, G, C, N}"+ "\n"

                + "@attribute Class { true, false}" + "\n"
                + "@data";
        return s;
    }

    /**
     * Return best entries for the beginning of arff file
     *
     * @return
     */
    public static String getvBestWekaFileEntries() {
        String s = ""
                + "@relation mirnaInPrecursor" + "\n"
                + "@attribute DistanceFromTerminalLoop real" + "\n"
                + "@attribute LengthOfOverlapInLoop real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow3 real" + "\n"
                + "@attribute BasePairsInDuplexMirnaMirnaStar real" + "\n"
                + "@attribute LengthBiggestBulgesPerc real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow5 real" + "\n"
                + "@attribute LenghtOfBiggestBulge real" + "\n"
                + "@attribute AverageNumberOfPairedBasesInWindow7 real" + "\n"
                + "@attribute MirnaIncludedInLoop { true, false}" + "\n"
                + "@attribute DistanceFromHairpinStart real" + "\n"
                + "@attribute StartOfPerfect5MerBasePair real" + "\n"
                + "@attribute MaximumLengthWithoutBulges real" + "\n"
                + "@attribute MaximumLengthWithoutBulgesPerc real" + "\n"
                + "@attribute mfe real" + "\n"
                + "@attribute Class { true, false}" + "\n"
                + "@data";
        return s;
    }

}
