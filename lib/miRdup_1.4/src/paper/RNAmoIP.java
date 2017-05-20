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
 * RNAmoIP
 */
package paper;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import miRdup.Tools;

/**
 *
 * @author mycky
 */
public class RNAmoIP {
    
    public static void main(String[] args) {
        try {
            BufferedReader br = new BufferedReader(new FileReader("all.txt"));
            PrintWriter pw = new PrintWriter(new FileWriter("all.rnamoip.txt"));
            String line;
            while (br.ready()){
                line = br.readLine();
                String name = line.split("\t")[0].trim();
                System.out.println(name);
                String seq = line.split("\t")[3].trim();
                String struc = line.split("\t")[4].trim();
                String strucMotifs = executeRNAmoIP(seq, struc);
                pw.println(line+"\t"+strucMotifs);
                pw.flush();
            }
            pw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        
    }

    
    public static String executeRNAmoIP(String sequence, String structure){
        System.out.println(structure);
        String optSol="";
        StringBuilder struct=new StringBuilder(structure);
        try {
            String output=Tools.executeLinuxCommand("export GUROBI_HOME=\"/home/mycky/tools/gurobi500/linux32\";"
                    + "export PATH=\"${PATH}:${GUROBI_HOME}/bin\";"
                    + "export LD_LIBRARY_PATH=\"${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib\";"
                    + "export GRB_LICENSE_FILE=\"/home/mycky/tools/gurobi500/gurobi.lic\";"
                    + "python RNAMoIP.py "+sequence+" \""+structure+"\" tools/NO_RED_DESC/ 0.3 4"); 
            String out[]=output.split("\n");
            for (int i = 0; i < out.length; i++) {
                if (out[i].trim().startsWith("Optimal solution nb:")){
                    for (int j = i; j < out.length; j++) {
                        if (out[j].trim().startsWith("D")){
                            int a = Integer.valueOf(out[j].split("-")[1]);
                            int b = Integer.valueOf(out[j].split("-")[2]);
                            struct.setCharAt(a, '.');
                            struct.setCharAt(b, '.');
                        }
                        
                    }
                }                
            }
        } catch (Exception e) {
            e.printStackTrace();
        }        
        
        structure = struct.toString();
        System.out.println(structure+"\n");
        
        return structure;
    }        
}
