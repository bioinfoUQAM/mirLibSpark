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
 * To parralise the calculs
 */
package paper;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

/**
 *
 * @author Mickael
 */
public class splitFiles {
    
    public static void main(String[] args) {
        String infile="all.txt";
        String infileName=infile.substring(0, infile.lastIndexOf("."));
        int linesPerFile=100;
        
        try {
            BufferedReader br = new BufferedReader(new FileReader(infile));
            PrintWriter pw = null ;
            int cpt=0;
            while (br.ready()){                
                if (cpt%linesPerFile==0){
                    if (cpt>0){
                        pw.flush();pw.close();
                    }
                    pw = new PrintWriter(new FileWriter("/ibrixfs1/Data/mik/mirdup/miRdup/meanMethod.noExp/"+infileName+"_"+cpt+".splitted.txt"));
                    
                }
                cpt++;
                pw.println(br.readLine());
            }
            pw.flush();pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
