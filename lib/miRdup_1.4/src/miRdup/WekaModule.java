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
 * Weka module
 */
package miRdup;

import java.awt.BorderLayout;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import weka.attributeSelection.AttributeSelection;
import weka.classifiers.Classifier;
import weka.classifiers.Evaluation;
import weka.classifiers.evaluation.ThresholdCurve;
import weka.classifiers.meta.FilteredClassifier;
import weka.filters.Filter;
import weka.core.FastVector;
import weka.core.Instances;
import weka.core.Range;
import weka.core.Utils;
import weka.core.converters.ConverterUtils.DataSource;
import weka.filters.unsupervised.attribute.Remove;
import weka.gui.visualize.PlotData2D;
import weka.gui.visualize.ThresholdVisualizePanel;

/**
 *
 * @author Mickael Leclercq
 */
public class WekaModule {

    static DecimalFormat dec = new DecimalFormat();

    public static void trainModel(File arff, String keyword) {
        dec.setMaximumFractionDigits(3);
        System.out.println("\nTraining model on file " + arff);
        try {
            // load data
            DataSource source = new DataSource(arff.toString());
            Instances data = source.getDataSet();
            if (data.classIndex() == -1) {
                data.setClassIndex(data.numAttributes() - 1);
            }

            PrintWriter pwout = new PrintWriter(new FileWriter(keyword + Main.modelExtension + "Output"));
            PrintWriter pwroc = new PrintWriter(new FileWriter(keyword + Main.modelExtension + "roc.arff"));

            //remove ID row
            Remove rm = new Remove();
            rm.setAttributeIndices("1");
            FilteredClassifier fc = new FilteredClassifier();
            fc.setFilter(rm);

//            // train model svm
//            weka.classifiers.functions.LibSVM model = new weka.classifiers.functions.LibSVM();
//            model.setOptions(weka.core.Utils.splitOptions("-S 0 -K 2 -D 3 -G 0.0 -R 0.0 -N 0.5 -M 40.0 -C 1.0 -E 0.0010 -P 0.1 -B"));
            // train model MultilayerPerceptron
//            weka.classifiers.functions.MultilayerPerceptron model = new weka.classifiers.functions.MultilayerPerceptron();
//            model.setOptions(weka.core.Utils.splitOptions("-L 0.3 -M 0.2 -N 500 -V 0 -S 0 -E 20 -H a"));
            // train model Adaboost on RIPPER
//            weka.classifiers.meta.AdaBoostM1 model = new weka.classifiers.meta.AdaBoostM1();
//            model.setOptions(weka.core.Utils.splitOptions("weka.classifiers.meta.AdaBoostM1 -P 100 -S 1 -I 10 -W weka.classifiers.rules.JRip -- -F 10 -N 2.0 -O 5 -S 1"));
            // train model Adaboost on FURIA
//            weka.classifiers.meta.AdaBoostM1 model = new weka.classifiers.meta.AdaBoostM1();
//            model.setOptions(weka.core.Utils.splitOptions("weka.classifiers.meta.AdaBoostM1 -P 100 -S 1 -I 10 -W weka.classifiers.rules.FURIA -- -F 10 -N 2.0 -O 5 -S 1 -p 0 -s 0"));
            //train model Adaboot on J48 trees
//             weka.classifiers.meta.AdaBoostM1 model = new weka.classifiers.meta.AdaBoostM1();
//             model.setOptions(
//                     weka.core.Utils.splitOptions(
//                     "-P 100 -S 1 -I 10 -W weka.classifiers.trees.J48 -- -C 0.25 -M 2"));
            //train model Adaboot on Random Forest trees
            weka.classifiers.meta.AdaBoostM1 model = new weka.classifiers.meta.AdaBoostM1();
            model.setOptions(
                    weka.core.Utils.splitOptions(
                            "-P 100 -S 1 -I 10 -W weka.classifiers.trees.RandomForest -- -I 50 -K 0 -S 1"));

            if (Main.debug) {
                System.out.print("Model options: " + model.getClass().getName().trim() + " ");
            }
            System.out.print(model.getClass() + " ");
            for (String s : model.getOptions()) {
                System.out.print(s + " ");
            }

            pwout.print("Model options: " + model.getClass().getName().trim() + " ");
            for (String s : model.getOptions()) {
                pwout.print(s + " ");
            }

            //build model
//            model.buildClassifier(data);
            fc.setClassifier(model);
            fc.buildClassifier(data);

            // cross validation 10 times on the model
            Evaluation eval = new Evaluation(data);
            //eval.crossValidateModel(model, data, 10, new Random(1));
            StringBuffer sb = new StringBuffer();
            eval.crossValidateModel(fc, data, 10, new Random(1), sb, new Range("first,last"), false);

            //System.out.println(sb);
            pwout.println(sb);
            pwout.flush();

            // output
            pwout.println("\n" + eval.toSummaryString());
            System.out.println(eval.toSummaryString());

            pwout.println(eval.toClassDetailsString());
            System.out.println(eval.toClassDetailsString());

            //calculate importants values
            String ev[] = eval.toClassDetailsString().split("\n");

            String ptmp[] = ev[3].trim().split(" ");
            String ntmp[] = ev[4].trim().split(" ");
            String avgtmp[] = ev[5].trim().split(" ");

            ArrayList<String> p = new ArrayList<String>();
            ArrayList<String> n = new ArrayList<String>();
            ArrayList<String> avg = new ArrayList<String>();

            for (String s : ptmp) {
                if (!s.trim().isEmpty()) {
                    p.add(s);
                }
            }
            for (String s : ntmp) {
                if (!s.trim().isEmpty()) {
                    n.add(s);
                }
            }
            for (String s : avgtmp) {
                if (!s.trim().isEmpty()) {
                    avg.add(s);
                }
            }

            double tp = Double.parseDouble(p.get(0));
            double fp = Double.parseDouble(p.get(1));
            double tn = Double.parseDouble(n.get(0));
            double fn = Double.parseDouble(n.get(1));
            double auc = Double.parseDouble(avg.get(7));

            pwout.println("\nTP=" + tp + "\nFP=" + fp + "\nTN=" + tn + "\nFN=" + fn);
            System.out.println("\nTP=" + tp + "\nFP=" + fp + "\nTN=" + tn + "\nFN=" + fn);

            //specificity, sensitivity, Mathew's correlation, Prediction accuracy
            double sp = ((tn) / (tn + fp));
            double se = ((tp) / (tp + fn));
            double acc = ((tp + tn) / (tp + tn + fp + fn));
            double mcc = ((tp * tn) - (fp * fn)) / Math.sqrt((tp + fp) * (tn + fn) * (tp + fn) * tn + fp);

            String output
                    = "\nse=" + dec.format(se).replace(",", ".")
                    + "\nsp=" + dec.format(sp).replace(",", ".")
                    + "\nACC=" + dec.format(acc).replace(",", ".")
                    + "\nMCC=" + dec.format(mcc).replace(",", ".")
                    + "\nAUC=" + dec.format(auc).replace(",", ".");

            pwout.println(output);
            System.out.println(output);

            pwout.println(eval.toMatrixString());
            System.out.println(eval.toMatrixString());

            pwout.flush();
            pwout.close();

            //Saving model
            System.out.println("Model saved: " + keyword + Main.modelExtension);
            weka.core.SerializationHelper.write(keyword + Main.modelExtension, fc.getClassifier() /*model*/);

            // get curve
            ThresholdCurve tc = new ThresholdCurve();
            int classIndex = 0;
            Instances result = tc.getCurve(eval.predictions(), classIndex);
            pwroc.print(result.toString());
            pwroc.flush();
            pwroc.close();

            // draw curve
            //rocCurve(eval);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void testModel(File testarff, String predictionsFile, String classifier, boolean predictMiRNA) {
        System.out.println("Testing model on " + predictionsFile + " adapted in " + testarff + ". Submitted to model " + classifier);

        try {
            //add predictions sequences to object
            ArrayList<MirnaObject> alobj = new ArrayList<MirnaObject>();
            BufferedReader br = null;
            try {
                br = new BufferedReader(new FileReader(predictionsFile + ".folded"));
            } catch (FileNotFoundException fileNotFoundException) {
                br = new BufferedReader(new FileReader(predictionsFile));
            }
            BufferedReader br2 = new BufferedReader(new FileReader(testarff));
            String line2 = br2.readLine();
            while (!line2.startsWith("@data")) {
                line2 = br2.readLine();
            }
            String line = " ";
            int cpt = 0;
            while (br.ready()) {
                line = br.readLine();
                line2 = br2.readLine();
                String[] tab = line.split("\t");
                MirnaObject m = new MirnaObject();
                m.setArff(line2);
                m.setId(cpt++);
                m.setIdName(tab[0]);
                m.setMatureSequence(tab[1]);
                m.setPrecursorSequence(tab[2]);
                m.setStructure(tab[3]);
                alobj.add(m);
            }
            br.close();
            br2.close();

            // load data
            DataSource source = new DataSource(testarff.toString());
            Instances data = source.getDataSet();
            if (data.classIndex() == -1) {
                data.setClassIndex(data.numAttributes() - 1);
            }
            //remove ID row
            data.deleteAttributeAt(0);
            //load model
            Classifier model = (Classifier) weka.core.SerializationHelper.read(classifier);

            // evaluate dataset on the model
            Evaluation eval = new Evaluation(data);

            eval.evaluateModel(model, data);

            FastVector fv = eval.predictions();

            // output
            File classifFile = new File(classifier);
            PrintWriter pw = new PrintWriter(new FileWriter(predictionsFile + "." + classifFile.getName() + ".miRdup.txt"));
            PrintWriter pwt = new PrintWriter(new FileWriter(predictionsFile + "." + classifFile.getName() + ".miRdup.tab.txt"));
            PrintWriter pwout = new PrintWriter(new FileWriter(predictionsFile + "." + classifFile.getName() + ".miRdupOutput.txt"));

            for (int i = 0; i < fv.size(); i++) {
                //System.out.println(fv.elementAt(i).toString());
                String[] tab = fv.elementAt(i).toString().split(" ");
                int actual = Integer.valueOf(tab[1].substring(0, 1));
                int predicted = Integer.valueOf(tab[2].substring(0, 1));
                double score = 0.0;
                boolean validated = false;
                if (actual == predicted) { //case validated
                    int s = tab[4].length();
                    try {
                        score = Double.valueOf(tab[4]);
                        //score = Double.valueOf(tab[4].substring(0, s - 1));
                    } catch (NumberFormatException numberFormatException) {
                        score = 0.0;
                    }

                    validated = true;
                } else {// case not validated
                    int s = tab[5].length();
                    try {
                        score = Double.valueOf(tab[5]);
                        //score = Double.valueOf(tab[5].substring(0, s - 1));
                    } catch (NumberFormatException numberFormatException) {
                        score = 0.0;
                    }
                    validated = false;
                }
                MirnaObject m = alobj.get(i);
                m.setActual(actual);
                m.setPredicted(predicted);
                m.setScore(score);
                m.setValidated(validated);
                m.setNeedPrediction(predictMiRNA);
                String predictionMiRNA = "";
                if (predictMiRNA && validated == false) {
                    predictionMiRNA = miRdupPredictor.Predictor.predictionBySequence(
                            m.getPrecursorSequence(), classifier, classifFile.getName() + ".miRdupPrediction.txt");
                    try {
                        m.setPredictedmiRNA(predictionMiRNA.split(",")[0]);
                        m.setPredictedmiRNAstar(predictionMiRNA.split(",")[1]);
                    } catch (Exception e) {
                        m.setPredictedmiRNA(predictionMiRNA);
                        m.setPredictedmiRNAstar(predictionMiRNA);
                    }
                }

                pw.println(m.toStringFullPredictions());
                pwt.println(m.toStringPredictions());
                if (i % 100 == 0) {
                    pw.flush();
                    pwt.flush();
                }
            }

            //System.out.println(eval.toSummaryString("\nSummary results of predictions\n======\n", false));
            String[] out = eval.toSummaryString("\nSummary results of predictions\n======\n", false).split("\n");
            String info = out[0] + "\n" + out[1] + "\n" + out[2] + "\n" + out[4] + "\n" + out[5] + "\n" + out[6] + "\n" + out[7] + "\n" + out[11] + "\n";
            System.out.println(info);
            //System.out.println("Predicted position of the miRNA by miRdup:"+predictionMiRNA);
            pwout.println("File " + predictionsFile + " adapted in " + testarff + " submitted to model " + classifier);
            pwout.println(info);

            pw.flush();
            pw.close();
            pwt.flush();
            pwt.close();
            pwout.flush();
            pwout.close();

            System.out.println("Results in " + predictionsFile + "." + classifFile.getName() + ".miRdup.txt");

            // draw curve
            //rocCurve(eval);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public static String testModel(File testarff, String classifier) {
        // System.out.println("Testing model on "+testarff+". Submitted to model "+classifier);
        try {

            // load data
            DataSource source = new DataSource(testarff.toString());
            Instances data = source.getDataSet();
            if (data.classIndex() == -1) {
                data.setClassIndex(data.numAttributes() - 1);
            }

            //load model
            Classifier model = (Classifier) weka.core.SerializationHelper.read(classifier);

            // evaluate dataset on the model
            Evaluation eval = new Evaluation(data);

            eval.evaluateModel(model, data);
            FastVector fv = eval.predictions();

            //calculate importants values
            String ev[] = eval.toClassDetailsString().split("\n");

            String p = ev[3].trim();
            String n = ev[4].trim();

            double tp = Double.parseDouble(p.substring(0, 6).trim());
            double fp = 0;
            try {
                fp = Double.parseDouble(p.substring(11, 16).trim());
            } catch (Exception exception) {
                fp = Double.parseDouble(p.substring(7, 16).trim());
            }
            double tn = Double.parseDouble(n.substring(0, 6).trim());
            double fn = 0;
            try {
                fn = Double.parseDouble(n.substring(11, 16).trim());
            } catch (Exception exception) {
                fn = Double.parseDouble(n.substring(7, 16).trim());
            }

            //System.out.println("\nTP="+tp+"\nFP="+fp+"\nTN="+tn+"\nFN="+fn);
            //specificity, sensitivity, Mathew's correlation, Prediction accuracy
            double sp = ((tn) / (tn + fp));
            double se = ((tp) / (tp + fn));
            double acc = ((tp + tn) / (tp + tn + fp + fn));
            double mcc = ((tp * tn) - (fp * fn)) / Math.sqrt((tp + fp) * (tn + fn) * (tp + fn) * tn + fp);
//            System.out.println("\nse="+se+"\nsp="+sp+"\nACC="+dec.format(acc).replace(",", ".")+"\nMCC="+dec.format(mcc).replace(",", "."));
//            System.out.println(eval.toMatrixString());

            String out = dec.format(acc).replace(",", ".");
            System.out.println(out);
            return out;
        } catch (Exception e) {
            e.printStackTrace();
            return "";
        }

    }

    public static void attributeSelection(File arff, String outfile) {
        // load data
        try {
            PrintWriter pw = new PrintWriter(new FileWriter(outfile));
            DataSource source = new DataSource(arff.toString());
            Instances data = source.getDataSet();
            if (data.classIndex() == -1) {
                data.setClassIndex(data.numAttributes() - 1);
            }

            AttributeSelection attrsel = new AttributeSelection();
            weka.attributeSelection.InfoGainAttributeEval eval = new weka.attributeSelection.InfoGainAttributeEval();

            weka.attributeSelection.Ranker rank = new weka.attributeSelection.Ranker();
            rank.setOptions(
                    weka.core.Utils.splitOptions(
                            "-T -1.7976931348623157E308 -N -1"));
            if (Main.debug) {
                System.out.print("Model options: " + rank.getClass().getName().trim() + " ");
            }
            for (String s : rank.getOptions()) {
                System.out.print(s + " ");
            }
            attrsel.setEvaluator(eval);
            attrsel.setSearch(rank);
            attrsel.setFolds(10);

            attrsel.SelectAttributes(data);
            //attrsel.CrossValidateAttributes();

            System.out.println(attrsel.toResultsString());
            pw.println(attrsel.toResultsString());

            //evaluation.crossValidateModel(classifier, data, 10, new Random(1));
            pw.flush();
            pw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void rocCurve(Evaluation eval) {
        try {
            // generate curve
            ThresholdCurve tc = new ThresholdCurve();
            int classIndex = 0;
            Instances result = tc.getCurve(eval.predictions(), classIndex);
            result.toString();
            // plot curve
            ThresholdVisualizePanel vmc = new ThresholdVisualizePanel();
            vmc.setROCString("(Area under ROC = "
                    + Utils.doubleToString(tc.getROCArea(result), 4) + ")");
            vmc.setName(result.relationName());
            PlotData2D tempd = new PlotData2D(result);
            tempd.setPlotName(result.relationName());
            tempd.addInstanceNumberAttribute();
            // specify which points are connected
            boolean[] cp = new boolean[result.numInstances()];
            for (int n = 1; n < cp.length; n++) {
                cp[n] = true;
            }
            tempd.setConnectPoints(cp);
            // add plot
            vmc.addPlot(tempd);

            //
            result.toString();

            // display curve
            String plotName = vmc.getName();
            final javax.swing.JFrame jf
                    = new javax.swing.JFrame("Weka Classifier Visualize: " + plotName);
            jf.setSize(500, 400);
            jf.getContentPane().setLayout(new BorderLayout());
            jf.getContentPane().add(vmc, BorderLayout.CENTER);
            jf.addWindowListener(new java.awt.event.WindowAdapter() {
                public void windowClosing(java.awt.event.WindowEvent e) {
                    jf.dispose();
                }
            });

            jf.setVisible(true);
            System.out.println("");
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

}
