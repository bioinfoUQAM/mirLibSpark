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
 * mirdup core
 */
package miRdup;

//import weka.filters.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;

/**
 *
 * @author Mickael Leclercq
 */
public class Main {

    public static boolean debug = false;
    public static boolean keywordSetted = false;

    //RNAfold and RNAduplex path
//    public static String rnafoldlinux="/ibrixfs1/Data/mik/tools/ViennaRNA-2.0.7/Progs/";//
//    public static String rnafoldlinux="/home/mycky/tools/ViennaRNA-2.0.6/Progs/";
    public static String rnafoldlinux = "";
    /**
     * Model name AF : all features BF : best features model type : SVM Adaboost
     * RF : random forest J48 : trees MLP : MultiLayerPerceptron RIP : Ripper
     * FURIA : Advanced ripper
     */
    // Train all attributes or best Features
    public static Boolean bestFeatures = false;
    public static String modelExtension = ".model";

    // search keyword
    public static String keyword = "all";

    // matures sequences
    private static String matures = "";

    // hairpins sequences
    private static String hairpins = "";

    // mirbase organism (used to search keyword)
    private static String organisms = "";

    // embl file
    private static String embl = "";

    // miRbase folded hairpins
    private static String structures;

    //Testing on a model, a model must be submitted
    private static String predictionsFile = "";
    private static String classifier = "";

    //arff file if already created
    private static String arffFileForTrain = "";
    private static String arffFileForValidate = "";

    //predict miRNA position
    public static boolean predictMirnaPosition = false;

    // offline mirbase if false. Files must be submitted
    private static Boolean trainFromOnlineMiRbase = true;

    // ONLY Predict of the miRNA from pre-miRNA
    private static boolean predictMiRNAFromPrec = false;
    private static String predictionInfile = "";
    private static String precursor = "";
    private static String model = "";
    private static String predictionOutfile = "";

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
//        args=tests();

        try {
            setOptions(args);
        } catch (Exception e) {
            help();
        }
        String os = System.getProperty("os.name").toLowerCase();
        if (!os.contains("win") && rnafoldlinux.isEmpty()) {
            System.out.println("You must submit the RNAfold program path (option -r). See help with -help option");
        } else {
            //check Vienna tools
            if (!Vienna.checkViennaTools()) {
                System.err.println("Vienna tools cannot be executed");
            } else if (predictMiRNAFromPrec) {
                mirdupExecutionPredictor();
            } else {
                miRdupExecutionEMBL();
            }
        }

    }

    /**
     * read options
     *
     * @param args
     */
    public static void setOptions(String[] args) {
        String cmd = " ";
        for (String s : args) {
            cmd += s + " ";
        }
        String[] options = cmd.split(" -");
        for (String s : options) {
            if (s.startsWith("help")) {
                help();
            }
            // get keyword
            if (s.startsWith("k")) {
                keyword = s.substring(1).trim();
                keywordSetted = true;
            }

            // get mirbase embl file
            if (s.startsWith("e")) {
                embl = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get mirbase mature file
            if (s.startsWith("m")) {
                matures = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get hairpins files
            if (s.startsWith("h")) {
                hairpins = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get organism file
            if (s.startsWith("o")) {
                organisms = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get structure file
            if (s.startsWith("s")) {
                structures = s.substring(1).trim();
            }

            //Validation
            if (s.startsWith("v")) {
                predictionsFile = s.substring(1).trim();
            }

            // classifier model file
            if (s.startsWith("c")) {
                classifier = s.substring(1).trim();
            }

            // get arff of training
            if (s.startsWith("a")) {
                arffFileForTrain = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get arff of test
            if (s.startsWith("b")) {
                arffFileForValidate = s.substring(1).trim();
                trainFromOnlineMiRbase = false;
            }

            // get rnafold path
            if (s.startsWith("r")) {
                rnafoldlinux = s.substring(1).trim();
            }

            //predict miRNA
            if (s.startsWith("p")) {
                predictMirnaPosition = true;
            }

            //prediction of the miRNA from precursor
            if (s.startsWith("predict")) {
                predictMiRNAFromPrec = true;
            }
            if (s.startsWith("u")) {
                precursor = s.substring(1).trim();
            }
            if (s.startsWith("d")) {
                model = s.substring(1).trim();
            }
            if (s.startsWith("f")) {
                predictionOutfile = s.substring(1).trim();
            }
            if (s.startsWith("i")) {
                predictionInfile = s.substring(1).trim();
            }
        }
    }

    /**
     * Use miRbase EMBL source EMBL permits to get experimental only
     */
    public static void miRdupExecutionEMBL() {

        //Execution depending options
        //if fastas are given
        if (!matures.isEmpty() && !hairpins.isEmpty()) {
            if (keywordSetted == false) {
                keyword = "";
            }
            miRdupExecutionFastas();
        } // Normal execution. We get sequences on mirbase using keyword ("all" if empty),
        // then we train the model on it, and if the test dataset is present we submit it to the model
        else if (!classifier.isEmpty()) {
            if (predictionsFile.isEmpty()) {
                System.out.println("You must submit a file with data you want to validate");
            } else {
                System.out.println("Validation of " + predictionsFile + " with classifier " + classifier);
                // test model
                File arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, bestFeatures);
                } else {
                    System.out.println("Using arff file " + arffFileForValidate + " to validate");
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, classifier, predictMirnaPosition);
            }
        } else if (trainFromOnlineMiRbase == true) {
            // search sequences from mirbase based on a keyword
            Mirbase m = new Mirbase();
            ArrayList altrain = m.getSequencesFromMirbaseEMBL(keyword);
            // train model
            File arff = new File(keyword + ".arff");
            AdaptDataForWeka.createFileFromList(altrain, arff, bestFeatures);
            WekaModule.trainModel(arff, keyword);

            // submit a file with predicted mirnas and their hairpins to the model
            if (!predictionsFile.isEmpty()) {
                // test model
                arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, bestFeatures);
                } else {
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, keyword + modelExtension, predictMirnaPosition);
            }
            // Case offline mirbase, and local sequences are submitted
        } else {
            // search sequences from mirbase based on a keyword
            File arff = null;
            if (arffFileForTrain.isEmpty()) {
                Mirbase m = new Mirbase();
                ArrayList altrain = m.getSequencesFromMirbaseEMBLFile(keyword, embl, structures, organisms);
                arff = new File(keyword + ".arff");
                AdaptDataForWeka.createFileFromList(altrain, arff, bestFeatures);
            } else {
                arff = new File(arffFileForTrain);
            }
            // train model
            WekaModule.trainModel(arff, keyword);

            // submit a file with predicted mirnas and their hairpins to the model
            if (!predictionsFile.isEmpty()) {
                // test model
                arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, bestFeatures);
                } else {
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, keyword + modelExtension, predictMirnaPosition);
            }
        }
    }

    /**
     * Use miRbase fastas source
     */
    public static void miRdupExecutionFastas() {
        String filename = matures.substring(0, matures.lastIndexOf("."));
        // Normal execution. We get sequences on mirbase using keyword ("all" if empty),
        // then we train the model on it, and if the test dataset is present we submit it to the model
        if (!classifier.isEmpty()) {
            if (predictionsFile.isEmpty()) {
                System.out.println("You must submit a file with data you want to validate");
            } else {
                System.out.println("Validation of " + predictionsFile + " with classifier " + classifier);
                // test model
                File arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, false);
                } else {
                    System.out.println("Using arff file " + arffFileForValidate + " to validate");
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, classifier, predictMirnaPosition);
            }
        } else if (trainFromOnlineMiRbase == true) {
            // search sequences from mirbase based on a keyword
            Mirbase m = new Mirbase();
            ArrayList altrain = m.getSequencesFromMirbase(keyword);
            // train model
            File arff = new File(keyword + ".arff");
            AdaptDataForWeka.createFileFromList(altrain, arff, false);
            WekaModule.trainModel(arff, keyword);

            // submit a file with predicted mirnas and their hairpins to the model
            if (!predictionsFile.isEmpty()) {
                // test model
                arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, false);
                } else {
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, keyword + modelExtension, predictMirnaPosition);
            }
            // Case offline mirbase, and local sequences are submitted
        } else {
            // search sequences from mirbase based on a keyword
            File arff = null;
            if (arffFileForTrain.isEmpty()) {
                Mirbase m = new Mirbase();
                ArrayList altrain = m.getSequencesFromFiles(matures, hairpins, organisms, structures, keyword, filename);
                arff = new File(keyword + "" + filename + ".arff");
                AdaptDataForWeka.createFileFromList(altrain, arff, false);
            } else {
                arff = new File(arffFileForTrain);
            }
            // train model
            WekaModule.trainModel(arff, keyword + "" + filename);
            // submit a file with predicted mirnas and their hairpins to the model
            if (!predictionsFile.isEmpty()) {
                // test model
                arff = null;
                if (arffFileForValidate.isEmpty()) {
                    arff = new File(predictionsFile + ".arff");
                    AdaptDataForWeka.createPredictionDataset(predictionsFile, false);
                } else {
                    arff = new File(arffFileForValidate);
                }
                WekaModule.testModel(arff, predictionsFile, keyword + filename + modelExtension, predictMirnaPosition);
            }
        }
    }

    public static void AttributeSelection() {
        // search sequences from mirbase based on a keyword
        System.out.println("Attribute selection");
        Mirbase m = new Mirbase();
        ArrayList altrain = m.getSequencesFromFiles(matures, hairpins, organisms, structures, keyword, "tmp");
        // train model
        File arff = new File(keyword + ".arff");
        AdaptDataForWeka.createFileFromList(altrain, arff, bestFeatures);
        System.out.println("Attribute selection for " + keyword);
        WekaModule.attributeSelection(arff, keyword + ".AttribSel.txt");

    }

    /**
     * Printing help
     */
    private static void help() {
        try {
            BufferedReader br = new BufferedReader(new FileReader("readme.txt"));
            while (br.ready()) {
                System.out.println(br.readLine());
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static void mirdupExecutionPredictor() {
        miRdupPredictor.Predictor.rnafold = rnafoldlinux;
        if (!predictionInfile.isEmpty()) {
            miRdupPredictor.Predictor.predictionByFile(predictionInfile, model, predictionOutfile);
        } else {
            String predictionMiRNA = miRdupPredictor.Predictor.predictionBySequence(
                    precursor, model, predictionOutfile);
            //System.out.println("Predicted miRNA:"+predictionMiRNA);
        }

    }

}
