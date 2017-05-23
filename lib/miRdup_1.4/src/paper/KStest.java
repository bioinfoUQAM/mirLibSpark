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
 * Kolmogorov smirnov TEST
 */
package paper;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 *
 * @author fyrox
 */
public class KStest {

    static int maxlines;
    static int features;
    static String folder = "gnuplot" + File.separator;
    static String plt = folder + "createGraphs.plt";
    static DecimalFormat dec = new DecimalFormat();
    static int[] total;

    public static void main(String[] args) {

        File gnuplot = new File(plt);
        if (gnuplot.exists()) {
            gnuplot.delete();
        }

        String infiles[] = {
            "all.arff",
            "mammal.arff",
            "Pisces.arff",
            "nematod.arff",
            "Arthropoda.arff",
            "Viridiplantae.arff"};
        analyseFiles(infiles);
        kstest(infiles);
        //executegnuplot();
    }

    private static void analyseFiles(String infiles[]) {
        int max = 0;
        int file = 0;
        for (int i = 0; i < infiles.length; i++) {
            try {
                BufferedReader br = new BufferedReader(new FileReader(infiles[i]));
                String line = br.readLine();
                int cpt = 0;
                while (br.ready()) {
                    line = br.readLine();
                    cpt++;
                }
                if (cpt > max) {
                    max = cpt;
                    file = i;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        maxlines = (max / 2);

        try {
            BufferedReader br = new BufferedReader(new FileReader(infiles[file]));
            int attributes = 0;
            String line = "";
            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith("@attribute")) {
                    attributes++;
                }
            }
            br.close();
            features = attributes;
        } catch (Exception e) {
        }
        System.out.println("Number of features = " + features);
    }

    private static void kstest(String[] infiles) {

        HashMap<String, Integer> totals = new HashMap<String, Integer>();
        Double tab[][] = new Double[maxlines][infiles.length];
        total = new int[infiles.length];
        for (int i = 0; i < features; i++) {
            String title = "";
            for (int j = 0; j < infiles.length; j++) {
                try {
                    BufferedReader br = new BufferedReader(new FileReader(infiles[j]));
                    String line = br.readLine();
                    int cpt = 0;
                    int linenumber = 0;   //to get column title
                    while (br.ready()) {
                        line = br.readLine();
                        if (line != null && !line.startsWith("@")) {
                            String classification = line.substring(line.lastIndexOf(",") + 1);
                            if (classification.equals("true")) {
                                try {
                                    tab[cpt][j] = Double.valueOf(line.split(",")[i]);
                                } catch (Exception e) {
                                    // case true/false
                                    if (line.split(",")[i].trim().equals("true")) {
                                        tab[cpt][j] = 1.0;
                                    } else if (line.split(",")[i].trim().equals("false")) {
                                        tab[cpt][j] = 0.0;
                                    } // case A U G C
                                    else if (line.split(",")[i].trim().equals("A")) {
                                        tab[cpt][j] = 1.0;
                                    } else if (line.split(",")[i].trim().equals("U")) {
                                        tab[cpt][j] = 2.0;
                                    } else if (line.split(",")[i].trim().equals("G")) {
                                        tab[cpt][j] = 3.0;
                                    } else if (line.split(",")[i].trim().equals("C")) {
                                        tab[cpt][j] = 4.0;
                                    } else {
                                        //e.printStackTrace();
                                    }
                                }
                                cpt++;
                            }
                        } else {
                            // get column title
                            if (linenumber == i) {
                                title = line.split(" ")[1];

                            }
                        }
                        linenumber++;
                    }
                    total[j] = cpt;
                    totals.put(infiles[j], cpt);
                } catch (Exception e) {
                    e.printStackTrace();
                }

            }
            if (!title.equals("id")) {
                dec.setMaximumFractionDigits(8);
                analyse(infiles, title, tab, totals);
            }
        }

    }

    /**
     *
     * @param infiles
     * @param title
     * @param rawdatatab
     */
    private static void analyse(String[] infiles, String title,
            Double[][] rawdatatab, HashMap<String, Integer> totals) {

        // System.out.println("Processing "+title);
        try {
            // analyse tab
            ArrayList<List<Map.Entry<Double, Double>>> al = new ArrayList(infiles.length);
            for (int j = 0; j < infiles.length; j++) {
                HashMap<Double, Double> hm = new HashMap<Double, Double>();
                for (int i = 0; i < rawdatatab.length; i++) {
                    if (rawdatatab[i][j] != null) {
                        double d = rawdatatab[i][j];
                        //  fill hasmap
                        if (hm.containsKey(d)) {
                            double value = hm.get(d);
                            value += 1.0;
                            hm.put(d, value);
                        } else {
                            hm.put(d, 1.0);
                        }
                    }
                }

                List<Map.Entry<Double, Double>> entries = new ArrayList<Map.Entry<Double, Double>>(hm.entrySet());

                Collections.sort(entries, new Comparator<Map.Entry<Double, Double>>() {
                    @Override
                    public int compare(final Map.Entry<Double, Double> e1, final Map.Entry<Double, Double> e2) {
                        return e1.getKey().compareTo(e2.getKey());
                    }
                });

                al.add(entries);
            }

            // get biggest hashmap
            int index = 0;
            int numberOfKeys = 0;
            for (int i = 0; i < al.size(); i++) {
                List<Map.Entry<Double, Double>> entries = al.get(i);
                if (entries.size() > numberOfKeys) {
                    numberOfKeys = entries.size();
                    index = i;
                }
            }
            List<Map.Entry<Double, Double>> biggestentries = al.get(index);
            // fill a new table
            Double tmptab[][] = new Double[1000][infiles.length + 1];
            //fill 1st column
            int cpt = 0;
            for (Map.Entry entry : biggestentries) {
                tmptab[cpt][0] = (Double) entry.getKey();
                cpt++;
            }
            //fill others columns
            for (int j = 1; j < infiles.length + 1; j++) {
                List<Map.Entry<Double, Double>> list = al.get(j - 1);
                Map<Double, Double> map = new HashMap<Double, Double>();
                for (Map.Entry<Double, Double> entry : list) {
                    map.put(entry.getKey(), entry.getValue());
                }
                for (int i = 0; i < numberOfKeys; i++) {
                    double key = tmptab[i][0];
                    if (map.get(key) == null) {
                        tmptab[i][j] = 0.0;
                    } else {
                        tmptab[i][j] = map.get(key);
                    }
                }
            }
            //get totals
            ArrayList<Double> altot = new ArrayList<Double>();
            for (int j = 1; j < infiles.length + 1; j++) {
                double total = 0.0;
                for (int i = 0; i < tmptab.length; i++) {
                    if (tmptab[i][j] != null) {
                        total += tmptab[i][j];
                    }
                }
                altot.add(total);
            }

            // transform table in percentages
            Double compiledTab[][] = new Double[numberOfKeys][infiles.length + 1];
            double maxperc = 0;
            for (int j = 0; j < infiles.length + 1; j++) {
                for (int i = 0; i < compiledTab.length; i++) {
                    if (j == 0) {
                        compiledTab[i][j] = tmptab[i][j];
                    } else {
                        double value = tmptab[i][j];
                        double perc = value * 100 / altot.get(j - 1);
                        compiledTab[i][j] = perc;
                        if (perc > maxperc) {
                            maxperc = perc;
                        }
                    }
                }
            }

            // Getting range
            String minRange = String.valueOf(compiledTab[0][0]);
            String maxRange = String.valueOf(compiledTab[compiledTab.length - 1][0]);

            if (title.toUpperCase().equals("LenghtOfBiggestBulge")) {
                System.out.println("");
            }
            //improve table
            boolean intervalled = false;
            boolean needRotation = false;
            String intervaltab[][] = null;
            if ((numberOfKeys >= 16 && !title.equals("length")
                    && !title.toUpperCase().equals("STARTOFPERFECT10MERBASEPAIR")
                    || title.toUpperCase().equals("U.(.")
                    || title.toUpperCase().equals("A.(.")
                    || title.toUpperCase().equals("C.(."))
                    && !title.equals("LenghtOfBiggestBulge")) {
                intervalled = true;
                double minvalue = compiledTab[1][0];
                double maxvalue = 0;
                for (int i = 0; i < compiledTab.length; i++) {
                    if (minvalue > compiledTab[i][0]) {
                        minvalue = compiledTab[i][0];
                    }
                    if (maxvalue < compiledTab[i][0]) {
                        maxvalue = compiledTab[i][0];
                    }
                }

                if (minvalue < -10) {
                    intervaltab = TransformIntervallesmfe(5, compiledTab, minvalue, maxvalue, infiles.length + 1);
                } //case where we have something like 1 to 15 with a lot of decimals (1.1, 2.3, ...)
                else if ((maxvalue) <= 10) {
                    intervaltab = TransformIntervalles(1, compiledTab, minvalue, maxvalue, infiles.length + 1);
                } else if ((maxvalue) < 50) {
                    intervaltab = TransformIntervalles(5, compiledTab, minvalue, maxvalue, infiles.length + 1);
                } //case 0 a plusieurs centaines
                else if (maxvalue >= 50 && maxvalue <= 100) {
                    intervaltab = TransformIntervalles(10, compiledTab, minvalue, maxvalue, infiles.length + 1);
                } else {
                    intervaltab = TransformIntervalles(10, compiledTab, minvalue, maxvalue, infiles.length + 1);
                }

                //update maxperc
                maxperc = 0;
                for (int j = 1; j < infiles.length + 1; j++) {
                    for (int i = 0; i < intervaltab.length; i++) {
                        double d = Double.parseDouble(intervaltab[i][j].replace(",", "."));
                        if (d > maxperc) {
                            maxperc = d;
                        }
                    }
                }
            }

            //regroup in one table
            String finaltab[][] = null;
            if (intervalled) {
                finaltab = new String[intervaltab.length][(infiles.length + 1)];
                for (int i = 0; i < intervaltab.length; i++) {
                    for (int j = 0; j < infiles.length + 1; j++) {
                        finaltab[i][j] = intervaltab[i][j];
                    }
                }
            } else {
                finaltab = new String[compiledTab.length][(infiles.length + 1)];
                for (int i = 0; i < compiledTab.length; i++) {
                    for (int j = 0; j < infiles.length + 1; j++) {
                        finaltab[i][j] = dec.format(compiledTab[i][j]);
                    }
                }
            }

            //do KS
            //critical D
            double n = 0;
            for (int i = 0; i < finaltab.length; i++) {
                n += Double.parseDouble(finaltab[i][1].replace(",", "."));
            }
            double critD = 1.36 / Math.sqrt(n);
            System.out.println(title + " n=" + n);
            //normalize
            Double normtab[][] = new Double[finaltab.length][(infiles.length)];
            for (int i = 0; i < finaltab.length; i++) {
                for (int j = 1; j < infiles.length + 1; j++) {
                    double value = Double.parseDouble(finaltab[i][j].replace(",", "."));
                    normtab[i][j - 1] = value / n;
                }
            }
            //cumulative
            Double cumultab[][] = new Double[normtab.length][(infiles.length)];
            for (int j = 0; j < infiles.length; j++) {
                for (int i = 0; i < normtab.length; i++) {
                    double normvalue = normtab[i][j];
                    try {
                        cumultab[i][j] = normvalue + cumultab[i - 1][j];
                    } catch (Exception e) {
                        cumultab[i][j] = normvalue + 0;
                    }

                }
            }
            //compare species between each other
            Double comparetab[][] = new Double[cumultab.length][15];
            int m = 0;
            for (int i = 0; i < infiles.length - 1; i++) {
                for (int j = i + 1; j < infiles.length; j++) {
                    for (int k = 0; k < cumultab.length; k++) {
                        double comp = Math.abs(cumultab[k][i] - cumultab[k][j]);
                        comparetab[k][m] = comp;
                    }
                    m++;
                }
            }
            dec.setMaximumFractionDigits(3);
            //check if a value is > to critD
            boolean ks = false;
            for (int i = 0; i < infiles.length; i++) {
                for (int j = 0; j < comparetab.length; j++) {
                    double value = comparetab[j][i];
                    if (value > critD) {
                        //System.out.println(title+" "+finaltab[j][0]+" "+j+" "+dec.format(value)+">"+dec.format(critD));
                        ks = true;
                    }
                }
            }
            if (ks) {
                System.out.println(title);
            }
            //System.out.println("");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static String[][] TransformIntervalles(int interval, Double[][] finaltab, double minvalue, double maxvalue, int columns) {
        int firstminimum = Math.round((float) minvalue);
        int lastmaximum = Math.round((float) Math.ceil(maxvalue));

        int numberOfLines = (Math.max(Math.abs(firstminimum), Math.abs(lastmaximum)) / interval);
        String intervaltab[][] = new String[numberOfLines][columns];

        //prepare 1st column
        ArrayList<String> al = new ArrayList<String>();
        for (int i = firstminimum; i < lastmaximum; i = i + interval) {
            al.add(i + ";" + (i + interval));
        }

        //fill 1st column
        for (int i = 0; i < intervaltab.length; i++) {
            intervaltab[i][0] = "]" + al.get(i) + "]";
        }

        //fill other columns
        for (int k = 0; k < intervaltab.length; k++) {
            String intervals = al.get(k);
            double intervala = Integer.valueOf(intervals.split(";")[0]);
            double intervalb = Integer.valueOf(intervals.split(";")[1]);
            Double line[] = new Double[columns];
            for (int j = 0; j < columns; j++) {
                line[j] = 0.0;
            }
            for (int i = 0; i < finaltab.length; i++) {
                if (finaltab[i][0] > intervala && finaltab[i][0] <= intervalb) {
                    for (int j = 1; j < columns; j++) {
                        line[j] += finaltab[i][j];
                    }
                }
            }
            for (int j = 1; j < columns; j++) {
                intervaltab[k][j] = String.valueOf(dec.format(line[j]));
            }
        }

        return intervaltab;
    }

    private static String[][] TransformIntervallesmfe(int interval, Double[][] finaltab, double minvalue, double maxvalue, int columns) {
        int firstminimum = Math.round((float) minvalue);
        int lastmaximum = Math.round((float) Math.ceil(maxvalue));

        int numberOfLines = (Math.max(Math.abs(firstminimum), Math.abs(lastmaximum)) / interval);
        String intervaltab[][] = new String[numberOfLines][columns];

        //prepare 1st column
        ArrayList<String> al = new ArrayList<String>();
        for (int i = firstminimum; i < lastmaximum; i = i + interval) {
            al.add(i + ";" + (i + interval));
        }

        //fill 1st column
        for (int i = 0; i < intervaltab.length; i++) {
            intervaltab[i][0] = "]" + al.get(i) + "]";
        }

        //fill other columns
        for (int k = 0; k < intervaltab.length; k++) {
            String intervals = al.get(k);
            double intervala = Integer.valueOf(intervals.split(";")[0]);
            double intervalb = Integer.valueOf(intervals.split(";")[1]);
            Double line[] = new Double[columns];
            for (int j = 0; j < columns; j++) {
                line[j] = 0.0;
            }
            for (int i = 0; i < finaltab.length; i++) {
                if (finaltab[i][0] > intervala && finaltab[i][0] <= intervalb) {
                    for (int j = 1; j < columns; j++) {
                        line[j] += finaltab[i][j];
                    }
                }
            }
            for (int j = 1; j < columns; j++) {
                intervaltab[k][j] = String.valueOf(dec.format(line[j]));
            }
        }
        return intervaltab;
    }

}
