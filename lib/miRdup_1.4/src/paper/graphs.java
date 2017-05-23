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
 * Generate GnuPlot scripts
 * 
 */
package paper;

import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import miRdup.Tools;

/**
 *
 * @author Mickael Leclercq
 */
public class graphs {
    
    static int maxlines;
    static int features;
    static String folder="gnuplot"+File.separator;
    static String plt=folder+"createGraphs.plt";
    static DecimalFormat dec = new DecimalFormat();
    static int[] total;

    
    public static void main(String[] args) {
        dec.setMaximumFractionDigits(2);
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
        generateCharts(infiles);
        //executegnuplot();
    }
    
    private static void analyseFiles(String infiles[]) {
        int max=0;
        int file=0;
        for (int i = 0; i < infiles.length; i++) {           
            try {
                BufferedReader br = new BufferedReader(new FileReader(infiles[i]));
                String line=br.readLine();
                int cpt=0;
                while(br.ready()){
                    line=br.readLine();
                    cpt++;
                }
                if (cpt>max) {
                    max=cpt;
                    file=i;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        maxlines=(max/2);
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(infiles[file]));
            int attributes=0;
            String line="";
            while(br.ready()){
                line=br.readLine();
                if (line.startsWith("@attribute")){
                    attributes++;
                }
            }
            br.close();
            features=attributes;            
        } catch (Exception e) {
        }
        System.out.println("Number of features = "+features);
    }
    
    private static void generateCharts(String[] infiles) {

        HashMap<String,Integer> totals = new HashMap<String,Integer>();
        Double tab[][]= new Double[maxlines][infiles.length];
        total=new int[infiles.length];
        for (int i = 0; i < features; i++) {
            String title="";
            for (int j = 0; j < infiles.length; j++) {
                try {
                    BufferedReader br = new BufferedReader(new FileReader(infiles[j]));
                    String line = br.readLine();
                    int cpt=0;
                    int linenumber=0;   //to get column title                    
                    while (br.ready()) {
                        line = br.readLine();                        
                        if (line!=null&&!line.startsWith("@")){    
                            String classification=line.substring(line.lastIndexOf(",")+1);                            
                            if (classification.equals("true")) {
                                try {
                                    tab[cpt][j] = Double.valueOf(line.split(",")[i]);
                                } catch (NumberFormatException n) {
                                    // case true/false
                                    if (line.split(",")[i].trim().equals("true")) {
                                        tab[cpt][j] = 1.0;
                                    } else if (line.split(",")[i].trim().equals("false")){
                                        tab[cpt][j] = 0.0;
                                    }  
                                    // case A U G C 
                                    else if (line.split(",")[i].trim().equals("A")) {
                                        tab[cpt][j] = 1.0;
                                    } else if (line.split(",")[i].trim().equals("U")) {
                                        tab[cpt][j] = 2.0;
                                    } else if (line.split(",")[i].trim().equals("G")) {
                                        tab[cpt][j] = 3.0;
                                    } else if (line.split(",")[i].trim().equals("C")) {
                                        tab[cpt][j] = 4.0;
                                    }
                                }
                                cpt++;
                            }
                        } else {
                            // get column title
                            if (linenumber==i){
                                title = line.split(" ")[1];
                                
                            }
                        }
                        linenumber++;
                    }
                    total[j]=cpt;
                    totals.put(infiles[j], cpt);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                
            }
            createGnuPlotFiles(infiles,title,tab,totals);
        }
       
    }


    /**
     * 
     * @param infiles
     * @param title
     * @param rawdatatab 
     */
    private static void createGnuPlotFiles(String[] infiles, String title, 
            Double[][] rawdatatab, HashMap<String,Integer> totals) {

        System.out.println("Processing "+title);
        try {
            // analyse tab            
            ArrayList<List<Entry<Double, Double>>> al = new ArrayList(infiles.length);                                    
            for (int j = 0; j < infiles.length; j++) {
                HashMap<Double,Double> hm = new HashMap<Double,Double>();
                for (int i = 0; i < rawdatatab.length; i++) {                    
                    if (rawdatatab[i][j]!=null) {
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

                List<Entry<Double, Double>> entries = new ArrayList<Entry<Double, Double>>(hm.entrySet());
                
                Collections.sort(entries, new Comparator<Entry<Double, Double>>() {                    
                    @Override
                    public int compare(final Entry<Double, Double> e1, final Entry<Double, Double> e2) {
                      return e1.getKey().compareTo(e2.getKey());
                    }
                    });                
                
                
                al.add(entries);
            }
            
            // get biggest hashmap
            int index=0;
            int numberOfKeys=0;
            for (int i = 0; i < al.size(); i++) {
                List<Entry<Double, Double>> entries = al.get(i);                
                if (entries.size()>numberOfKeys) {
                    numberOfKeys=entries.size();
                    index=i;
                }               
            }
            List<Entry<Double, Double>> biggestentries = al.get(index);
            // fill a new table
            Double tmptab[][]= new Double[1000][infiles.length+1];
                //fill 1st column
            int cpt=0;
            for (Entry entry : biggestentries) {
                tmptab[cpt][0]=(Double) entry.getKey();
                cpt++;
            }
                //fill others columns            
            for (int j = 1; j < infiles.length+1; j++) {
                List<Entry<Double, Double>> list =al.get(j-1);
                Map<Double, Double> map = new HashMap<Double, Double>();
                for (Entry<Double, Double> entry : list) {
                    map.put(entry.getKey(), entry.getValue());
                }
                for (int i = 0; i < numberOfKeys; i++) {                    
                    double key = tmptab[i][0];
                    if (map.get(key)==null) {
                        tmptab[i][j] = 0.0;
                    } else {
                        tmptab[i][j] = map.get(key);
                    }
                }                
            }
                //get totals   
            ArrayList<Double> altot= new ArrayList<Double>();
            for (int j = 1; j < infiles.length+1; j++) {
                double total=0.0;
                for (int i = 0; i < tmptab.length; i++) {
                    if (tmptab[i][j]!=null) {
                        total+=tmptab[i][j];
                    }                    
                }
                altot.add(total);
            }
            
            // transform table in percentages
            Double finaltab[][]= new Double[numberOfKeys][infiles.length+1];   
            double maxperc=0;
            for (int j = 0; j < infiles.length+1; j++) {
                for (int i = 0; i < finaltab.length; i++) {
                   if (j==0){
                       finaltab[i][j]=tmptab[i][j];
                   } else{
                       double value = tmptab[i][j];
                       double perc = value*100/altot.get(j-1);
                       finaltab[i][j]=perc;
                       if (perc>maxperc) {
                           maxperc=perc;
                       }
                   }                    
                }                
            }
            
            // Getting range
            String minRange=String.valueOf(finaltab[0][0]);
            String maxRange=String.valueOf(finaltab[finaltab.length-1][0]);
            
            if (title.toUpperCase().equals("LenghtOfBiggestBulge")){
                System.out.println("");
            }
            //improve table     
            boolean intervalled=false;
            boolean needRotation=false;
            String intervaltab[][]=null;
            if ((numberOfKeys>=16&&!title.equals("length")
                    &&!title.toUpperCase().equals("STARTOFPERFECT10MERBASEPAIR")
                    ||title.toUpperCase().equals("U.(.")                    
                    ||title.toUpperCase().equals("A.(.")
                    ||title.toUpperCase().equals("C.(."))
                    &&!title.equals("LenghtOfBiggestBulge")
                    ){
                intervalled=true;
                double minvalue=finaltab[1][0];
                double maxvalue=0;
                for (int i = 0; i < finaltab.length; i++) {
                    if (minvalue>finaltab[i][0]) minvalue=finaltab[i][0];
                    if (maxvalue<finaltab[i][0]) maxvalue=finaltab[i][0];                    
                }
                
                if (minvalue<-10) {
                    intervaltab= TransformIntervallesmfe(5,finaltab,minvalue,maxvalue,infiles.length+1);
                }
                //case where we have something like 1 to 15 with a lot of decimals (1.1, 2.3, ...)                                 
                else if ((maxvalue)<=10){
                    intervaltab= TransformIntervalles(1,finaltab,minvalue,maxvalue,infiles.length+1);
                } 
                else if ((maxvalue)<50){
                    intervaltab= TransformIntervalles(5,finaltab,minvalue,maxvalue,infiles.length+1);
                }
                //case 0 a plusieurs centaines
                else if (maxvalue>=50&&maxvalue<=100){
                    intervaltab= TransformIntervalles(10,finaltab,minvalue,maxvalue,infiles.length+1);
                }
                else {
                    intervaltab= TransformIntervalles(10,finaltab,minvalue,maxvalue,infiles.length+1);
                }
                
                //update maxperc
                maxperc=0;
                for (int j = 1; j < infiles.length+1; j++) {
                    for (int i = 0; i < intervaltab.length; i++) {
                        double d = 0;
                        try {
                            d = Double.parseDouble(intervaltab[i][j].replace(",", "."));
                        } catch (NumberFormatException numberFormatException) {
                            d=0;
                        }
                        if (d>maxperc) {
                            maxperc=d;
                        }                  
                    }                
                }
            }
            
            
            //print headers
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(folder+title+".dat", false)));
            pw.println("#"+title.toUpperCase()+"\n#");
            String headers=title;
            for (String infile : infiles) {
                int t = totals.get(infile);
                String h = infile.replace(".arff", "").replace("all", "miRbase");
                headers+="\t"+h+"("+t+")";
            }
            pw.println(headers+""+headers.substring(headers.indexOf("\t")));
            pw.flush();
            
            //add errorbars
            String finaltabWithErrors[][]=null;
            if (intervalled){
                finaltabWithErrors= new String[intervaltab.length][(infiles.length*2+1)];
                for (int i = 0; i < intervaltab.length; i++) {
                    for (int j = 0; j < infiles.length+1; j++) {                        
                        if (j>0) {
                            double v = 0;
                            try {
                                v = Double.parseDouble(intervaltab[i][j].replace(",", "."));
                            } catch (NumberFormatException numberFormatException) {
                                v=0;
                            }
                            double err = 1.96 * Math.sqrt((v / 100) * (1 - (v / 100)) / total[j-1]) * 100;
                            finaltabWithErrors[i][j + infiles.length] = dec.format(err);
                        } 
                        finaltabWithErrors[i][j]=intervaltab[i][j];                         
                    }
                }
            } else {
                finaltabWithErrors= new String[finaltab.length][(infiles.length*2+1)];
                for (int i = 0; i < finaltab.length; i++) {
                    for (int j = 0; j < infiles.length+1; j++) {                        
                        if (j>0) {
                            double v=finaltab[i][j];
                            double err = 1.96 * Math.sqrt((v / 100) * (1 - (v / 100)) / total[j-1]) * 100;
                            finaltabWithErrors[i][j + infiles.length] = dec.format(err);
                        }
                        finaltabWithErrors[i][j]=dec.format(finaltab[i][j]);                         
                    }
                }
            }
            
            //print table
            
            //case A U G C
            if (title.startsWith("nucleotide")){
                for (int i = 0; i < finaltabWithErrors.length; i++) {
                    if (i==0) pw.print("A\t");
                    if (i==1) pw.print("U\t");
                    if (i==2) pw.print("G\t");
                    if (i==3) pw.print("C\t");
                    for (int j = 1; j < infiles.length*2+1; j++) {
                        pw.print(finaltabWithErrors[i][j]+"\t");                    
                    }
                    pw.println();
                }
            } 
            //case true/false
            else if (finaltabWithErrors.length==2&&
                    !title.toUpperCase().startsWith("C")
                        &&!title.toUpperCase().startsWith("G")
                        &&!title.toUpperCase().startsWith("A")
                        &&!title.toUpperCase().startsWith("U")){                                
                for (int i = 0; i < finaltabWithErrors.length; i++) {
                    if (i == 0) {
                        pw.print("false\t");
                    }
                    if (i == 1) {
                        pw.print("true\t");
                    }
                    for (int j = 1; j < infiles.length * 2 + 1; j++) {
                        pw.print(finaltabWithErrors[i][j] + "\t");                            
                    }
                    pw.println();
                }                
            }
            // other cases
            else {
                for (int i = 0; i < finaltabWithErrors.length; i++) {
                    for (int j = 0; j < infiles.length*2+1; j++) {
                        pw.print(finaltabWithErrors[i][j].replace(",", ".") +"\t");                    
                    }
                    pw.println();
                }
            }
            pw.flush();pw.close();
            
            //print gnuplot script            
            PrintWriter pwplt = new PrintWriter(new BufferedWriter(new FileWriter(plt, true)));
            String rotation="";
            
            if(needRotation){
                rotation="rotate by -45";
            }
            String plot="plot '"+title+".dat' using 2:8:xtic(1)";
            for (int i = 1; i < infiles.length; i++) {
                plot+=", '' u "+(i+2)+":"+(i+8)+":xtic(1)";
                
            }
            
            //set xlabel
            String xlabel;
            if (title.equals("mfe")){
                xlabel="kcal/mol";
            } else if (title.equals("PercBasedPair")){
                xlabel="Percentage of miRNA bases pairs";
            } else if (title.equals("NumberOfBulges")){
                xlabel="number of bulges";
            } else if (title.toLowerCase().contains("perc")||title.length()<5){
                xlabel="Percentage of miRNA length";
            } else if (title.toLowerCase().contains("position")
                    ||title.toLowerCase().contains("included")
                    ||title.toLowerCase().contains("presence")){
                xlabel="false/true";
            } else xlabel="nucleotides";
            
            String range = "\\nRange: "+minRange+" to "+maxRange;
            if (xlabel.equals("false/true")){
                range="";
            }
            
            String newtitle=changeTitle(title);
            pwplt.println("#"+title.toUpperCase()+"\n"
                    + "set terminal pngcairo enhanced font \"arial,15\" fontscale 2.0 size 1600,900 \n"
                    + "set output '"+title+".png'\n"
                    + "set boxwidth 0.9 absolute\n"
                    + "set style fill solid 0.7 border lt -1\n"
                    + "set key inside right top vertical Right noreverse noenhanced autotitles columnhead nobox\n"
                    + "set grid ytics\n"
                    + "set nokey\n"
                    + "set style histogram errorbars linewidth 1 gap 3 title offset character 0, 0, 0\n"
                    + "set datafile missing '-'\n"
                    + "set style data histograms\n"
                    + "set xtics border in scale 1,0.5 nomirror "+rotation+" offset character 0, 0, 0 autojustify\n"
                    + "set xtics norangelimit font \",12\"\n"
                    + "set xtics ()\n"
                    + "set xlabel \""+xlabel+range+"\"\n"
                    + "set ylabel \"Percentage (%)\"\n"
                    + "set title \""+newtitle+"\" \n"
                    + "set bars 0.3 front\n"
                    + "set datafile separator \"\t\"\n"
                    + "set yrange [ 0 : "+(maxperc+5)+" ] noreverse nowriteback\n"
                    + plot+"\n");
            
           pwplt.flush();pwplt.close(); 
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static String[][] TransformIntervalles(int interval, Double[][] finaltab, double minvalue, double maxvalue, int columns) {
        int firstminimum = Math.round((float)minvalue);
        int lastmaximum = Math.round((float)Math.ceil(maxvalue));
        
        int numberOfLines=(Math.max(Math.abs(firstminimum), Math.abs(lastmaximum))/interval);        
        String intervaltab[][]= new String[numberOfLines][columns];

        //prepare 1st column
        ArrayList<String> al = new ArrayList<String>();
        for (int i = firstminimum; i < lastmaximum; i=i+interval) {
            al.add(i+";"+(i+interval));            
        }
        
        //fill 1st column
        for (int i = 0; i < intervaltab.length; i++) {
            intervaltab[i][0]="]"+al.get(i)+"]";            
        }
        
        //fill other columns
        for (int k = 0; k < intervaltab.length; k++) {
            String intervals= al.get(k);
            double intervala = Integer.valueOf(intervals.split(";")[0]);
            double intervalb = Integer.valueOf(intervals.split(";")[1]);
            Double line[] = new Double[columns];
            for (int j = 0; j < columns; j++) {
                line[j]=0.0;                
            }
            for (int i = 0; i < finaltab.length; i++) {
                if (finaltab[i][0]>intervala&&finaltab[i][0]<=intervalb){
                    for (int j = 1; j < columns; j++) {
                        line[j]+=finaltab[i][j];
                    }
                } 
            }
            for (int j = 1; j < columns; j++) {
                intervaltab[k][j]=String.valueOf(dec.format(line[j]));                
            }                        
        }
        
        return intervaltab;
    }
    
    
    private static String[][] TransformIntervallesmfe(int interval, Double[][] finaltab, double minvalue, double maxvalue, int columns) {
        int firstminimum = Math.round((float)minvalue);
        int lastmaximum = Math.round((float)Math.ceil(maxvalue));
        
        int numberOfLines=(Math.max(Math.abs(firstminimum), Math.abs(lastmaximum))/interval);        
        String intervaltab[][]= new String[numberOfLines][columns];

        //prepare 1st column

        ArrayList<String> al = new ArrayList<String>();
        for (int i = firstminimum; i < lastmaximum; i=i+interval) {
            al.add(i+";"+(i+interval));            
        }
        
        //fill 1st column
        for (int i = 0; i < intervaltab.length; i++) {
            intervaltab[i][0]="]"+al.get(i)+"]";           
        }
        
        //fill other columns
        for (int k = 0; k < intervaltab.length; k++) {
            String intervals= al.get(k);
            double intervala = Integer.valueOf(intervals.split(";")[0]);
            double intervalb = Integer.valueOf(intervals.split(";")[1]);
            Double line[] = new Double[columns];
            for (int j = 0; j < columns; j++) {
                line[j]=0.0;                
            }
            for (int i = 0; i < finaltab.length; i++) {
                if (finaltab[i][0]>intervala&&finaltab[i][0]<=intervalb){
                    for (int j = 1; j < columns; j++) {
                        line[j]+=finaltab[i][j];
                    }
                } 
            }
            for (int j = 1; j < columns; j++) {
                intervaltab[k][j]=String.valueOf(dec.format(line[j]));                
            }                        
        }        
        return intervaltab;
    }

    private static void executegnuplot() {
        String os=System.getProperty("os.name").toLowerCase();
        if (os.startsWith("win")){
            try {
                File f = new File(plt);
                System.out.println("executing command: gnuplot "+f.getCanonicalPath());
                Tools.executeWindowsCommand("gnuplot "+f.getCanonicalPath());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        } else{
            try {
                File f = new File(plt);
                System.out.println("executing command: gnuplot "+f.getCanonicalPath());
                Tools.executeLinuxCommand("gnuplot "+f.getCanonicalPath());
            } catch (Exception ex) {
                ex.printStackTrace();
            }
        }
    }


    public static String changeTitle(String oldtitle){
        String newtitle="";
        
        if (oldtitle.equals("length")){
            newtitle="Length of miRNAs";
        } else if (oldtitle.equals("mfe")){
            newtitle="Minimum free energy";
        } else if (oldtitle.equals("GCperc")){
            newtitle="GC content";
        } else if (oldtitle.equals("MaximumLengthWithoutBulgesPerc")){
            newtitle="Length of largest bulge-free region";
        } else if (oldtitle.equals("StartOfPerfect10MerBasePair")){
            newtitle="Position of first bulge-free region of more than 10bp";
        } else if (oldtitle.equals("DistanceFromTerminalLoop")){
            newtitle="Distance from terminal loop";
        } else if (oldtitle.equals("NumberOfBulges")){
            newtitle="Number of bulges";
        } else if (oldtitle.equals("LenghtOfBiggestBulge")){
            newtitle="Lenght of largest bulge";
        } else if (oldtitle.equals("nucleotideAt0")){
            newtitle="nucleotide at position 0(start of the miRNA)";
        }else if (oldtitle.equals("")){
            newtitle="";
        }else if (oldtitle.equals("")){
            newtitle="";
        }else if (oldtitle.equals("")){
            newtitle="";
        }else {
            newtitle=oldtitle;
        }
        
        return newtitle;
    }
            
    
}
