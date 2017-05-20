

/*
 * Kolmogorov smirnov TEST
 */
package paper;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;

/**
 *
 * @author fyrox
 */


public class KStest_wheat {
    static DecimalFormat dec = new DecimalFormat();
    
    public static void main(String[] args) {
        //String infile="wheat.mDup35Exp_90_1020_lgth_unique_perce_dist";
        //String infile="wheat.all_lib.smallRNAs_lgth_perce_unique_dist";
        //String infile="wheat.mDup35Exp_90_1020_deg_lgth_perce_dist";
        String infile="all.lib.reads.lengths.txt";
        
        String outfile=infile+"_kstest.txt";
        
        try {
        //check number of lines and columns
            PrintWriter pw = new PrintWriter(new FileWriter(outfile));
            int lines = 0;
            int columns = 0;
            try {
                
                BufferedReader br = new BufferedReader(new FileReader(infile));
                String line = br.readLine(); //jump header
                columns = line.split("\t").length - 1;
                while (br.ready()) {
                    line = br.readLine();
                    lines++;
                }                
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            //fill table
            Double table[][] = new Double[lines][columns];
            try {
                BufferedReader br = new BufferedReader(new FileReader(infile));
                String line = br.readLine(); //jump header
                int cpt = 0;
                while (br.ready()) {
                    line = br.readLine();
                    String tmp[] = line.split("\t");
                    for (int i = 1; i < tmp.length; i++) {
                        table[cpt][i - 1] = Double.valueOf(tmp[i]);                        
                    }
                    cpt++;
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

            //do KS
            //critical D
            double n = 0;
            for (int i = 0; i < table.length; i++) {
                n += table[i][1];                
            }
            double critD = 1.36 / Math.sqrt(n);

            //normalize
            Double normtab[][] = new Double[table.length][columns];
            for (int i = 0; i < table.length; i++) {
                for (int j = 0; j < columns; j++) {
                    double value = table[i][j];                    
                    normtab[i][j] = value / n;                    
                }
            }
            //cumulative
            Double cumultab[][] = new Double[table.length][columns];
            for (int j = 0; j < columns; j++) {
                for (int i = 0; i < normtab.length; i++) {                    
                    double normvalue = normtab[i][j];                    
                    try {
                        cumultab[i][j] = normvalue + cumultab[i - 1][j];
                    } catch (Exception e) {
                        cumultab[i][j] = normvalue + 0;
                    }
                    
                }
            }
            //compare columns between each other
            int totalComparisons = 0;
            for (int i = 1; i < columns; i++) {
                totalComparisons += i;
                
            }
            Double comparetab[][] = new Double[cumultab.length][totalComparisons];            
            int m = 0;
            for (int i = 0; i < columns - 1; i++) {                
                for (int j = i + 1; j < columns; j++) {
                    for (int k = 0; k < cumultab.length; k++) {                        
                        double comp = Math.abs(cumultab[k][i] - cumultab[k][j]);
                        comparetab[k][m] = comp;                        
                    }                    
                    m++;
                }                
            }
            dec.setMaximumFractionDigits(3);
            //print table comparetab
            //header
            for (int i = 0; i < columns - 1; i++) {
                for (int j = i + 1; j < columns; j++) {                    
                    System.out.print("lib" + (i + 1) + "-lib" + (j + 1) + "\t");    
                    pw.print("lib" + (i + 1) + "-lib" + (j + 1) + "\t"); 
                }                
            }
            System.out.println("");
            pw.println("");
            for (int j = 0; j < comparetab.length; j++) {
                for (int i = 0; i < totalComparisons; i++) {                    
                    double value = comparetab[j][i];
                    System.out.print(dec.format(value) + "\t");     
                    pw.print(dec.format(value) + "\t");
                }                
                System.out.println("");
                pw.println("");
            }
            System.out.println("\n\n");
            pw.println("\n\n");

            //check if a value is > to critD
            System.out.println("KS test. Comparative between cumulative values must be greater than critical D (1.36)"
                    + "\nTest is positive in : ");
            pw.println("KS test. Comparative between cumulative values must be greater than critical D (1.36)"
                    + "\nTest is positive in : ");
            boolean ks = false;
            for (int i = 0; i < columns; i++) {
                for (int j = 0; j < comparetab.length; j++) {
                    double value = comparetab[j][i];
                    if (value > critD) {
                        System.out.println("lib" + (i+1) + " vs lib" + (j+1) +  " " + dec.format(value) + ">" + dec.format(critD));
                        pw.println("lib" + (i+1) + " vs lib" + (j+1)  +  " " + dec.format(value) + ">" + dec.format(critD));
                        ks = true;
                    }                    
                }                
            }
            pw.flush();pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        
    }
    
    
    
}
