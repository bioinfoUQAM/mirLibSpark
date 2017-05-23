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
 * RNAfold executable
 */
package miRdup;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.Calendar;
import java.util.Date;

/**
 *
 * @author Mickael Leclercq
 */
public class Vienna {

    public static String executeRNAmoIP(String sequence, String structure) {
        System.out.println(structure);
        StringBuilder struct = new StringBuilder(structure);
        try {
            String output = Tools.executeLinuxCommand("export GUROBI_HOME=\"/home/mycky/tools/gurobi500/linux32\";"
                    + "export PATH=\"${PATH}:${GUROBI_HOME}/bin\";"
                    + "export LD_LIBRARY_PATH=\"${LD_LIBRARY_PATH}:${GUROBI_HOME}/lib\";"
                    + "export GRB_LICENSE_FILE=\"/home/mycky/tools/gurobi500/gurobi.lic\";"
                    + "python tools/RNAMoIP.py " + sequence + " \"" + structure + "\" tools/NO_RED_DESC/ 0.3 4");
            String out[] = output.split("\n");
            for (int i = 0; i < out.length; i++) {
                if (out[i].trim().startsWith("Optimal solution nb:")) {
                    for (int j = i; j < out.length; j++) {
                        if (out[j].trim().startsWith("D")) {
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
        System.out.println(structure);

        return structure;
    }

    /**
     *
     * @param sequence
     * @return structure
     */
    public static String GetSecondaryStructureTmp(String sequence) {
        if (Main.debug) {
            System.out.println("Getting secondary structure of sequence " + sequence);
        }
        String struct = "";
        String output = "";
        String mfe = "";
        String os = System.getProperty("os.name").toLowerCase();
        File f = null;
        if (os.startsWith("win")) {
            f = new File(/*c:/*/"tmpfold" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
        } else {
            f = new File("tmpfold" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
        }
        try {
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(f, false)));
            pw.println(sequence);
            pw.flush();
            pw.close();
            output = "";

            if (os.startsWith("win")) {
                output = Tools.executeWindowsCommand("RNAfold.exe --noPS < " + f.getAbsolutePath());
            } else {
                output = Tools.executeLinuxCommand(Main.rnafoldlinux + "RNAfold --noPS < " + f.getAbsolutePath());
            }
            //struct=output.replace("A", "").replace("U", "").replace("G", "").replace("C", "").trim();
            struct = output.split("\n")[1];
            mfe = struct.substring(struct.indexOf(" (") + 2, struct.length() - 1);
            struct = struct.replace("(" + mfe + ")", "");
        } catch (Exception e) {
            e.printStackTrace();
        }
        f.delete();
        return struct.trim();
    }

    /**
     *
     * @param sequence
     * @return structure
     */
    public static String GetSecondaryStructure(String sequence) {
        if (Main.debug) {
            System.out.println("Getting secondary structure of sequence " + sequence);
        }
        String struct = "";
        String output = "";
        String mfe = "";
        String os = System.getProperty("os.name").toLowerCase();

        try {
            String cmd = "";
            output = "";
            if (os.startsWith("win")) {
                //cmd="cmd /C RNAfold.exe --noPS";
                cmd = Main.rnafoldlinux + "\\RNAfold.exe --noPS";
            } else {
                cmd = Main.rnafoldlinux + "RNAfold --noPS";
            }
            Process p = Runtime.getRuntime().exec(cmd);
            new PrintWriter(p.getOutputStream(), true).println(sequence);
            BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
            new PrintWriter(p.getOutputStream(), true).println("@");
            p.waitFor();
            String line = "";
            while (bri.ready()) {
                line = bri.readLine();
                output += line + "\n";
            }
            struct = output.split("\n")[1];
            mfe = struct.substring(struct.indexOf(" (") + 2, struct.length() - 1);
            struct = struct.replace("(" + mfe + ")", "");
            bri.close();
            p.destroy();
        } catch (Exception e) {
            e.printStackTrace();
            return "0";
        }

        return struct.trim();
    }

    /**
     *
     * @param sequence
     * @return structure
     */
    public static File GetSecondaryStructuresFromFile(File infile) {
        if (Main.debug) {
            System.out.println("Getting secondary structure of file " + infile);
        }
        File outfile = new File(infile + ".folded");
        try {
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outfile, false)));
            pw.flush();
            pw.close();
            String os = System.getProperty("os.name").toLowerCase();
            if (os.startsWith("win")) {
                Tools.executeWindowsCommandToFile("RNAfold.exe --noPS < " + infile.getAbsolutePath() + ">" + outfile.getAbsolutePath());
            } else {
                //
                Tools.executeLinuxCommandToFile(Main.rnafoldlinux + "RNAfold --noPS < " + infile.getAbsolutePath() + ">" + outfile.getAbsolutePath());
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return outfile;
    }

    /**
     *
     * @param sequence
     * @return structure
     */
    public static File GetMfeDuplexFromFile(File infile) {
        if (Main.debug) {
            System.out.println("Getting secondary structure of file " + infile);
        }
        File outfile = new File(infile + ".folded");
        try {
            PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(outfile, false)));
            pw.flush();
            pw.close();
            String os = System.getProperty("os.name").toLowerCase();
            if (os.startsWith("win")) {
                Tools.executeWindowsCommandToFile("RNAduplex.exe < " + infile.getAbsolutePath() + ">" + outfile.getAbsolutePath());
            } else {
                Tools.executeLinuxCommandToFile(Main.rnafoldlinux + "RNAduplex < " + infile.getAbsolutePath() + ">" + outfile.getAbsolutePath());
            }

        } catch (Exception e) {
            e.printStackTrace();
        }
        return outfile;
    }

    public static String GetMfeDuplexRNAduplex(String mirna, String mirnastar) {
        String mfe = "";
        String output = "";
        try {
            if (mirnastar.equals("loop") || mirnastar.equals("error")) {
                return "0";
            }
            String os = System.getProperty("os.name").toLowerCase();
            String cmd = "";
            if (os.startsWith("win")) {
                cmd = "RNAduplex.exe";
            } else {
                cmd = Main.rnafoldlinux + "RNAduplex";
            }

            Process p = Runtime.getRuntime().exec(cmd);
            new PrintWriter(p.getOutputStream(), true).println(mirna);
            new PrintWriter(p.getOutputStream(), true).println(mirnastar);
            BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
            new PrintWriter(p.getOutputStream(), true).println("@");
            p.waitFor();
            String line = "";
            while (bri.ready()) {
                line = bri.readLine();
                output += line + "\n";
            }
            output = output.trim();
            mfe = output.substring(output.indexOf(" (") + 2, output.length() - 1);
            if (Double.valueOf(mfe) > 0) {
                mfe = "0";
            }
            bri.close();
            p.destroy();
            return mfe;
        } catch (Exception e) {
            return "0";
        }
    }

    public static ViennaObject GetInfosDuplexRNAcofold(String mirna, String mirnastar) {
        ViennaObject vo = new ViennaObject();
        String output = "";
        try {
            if (mirnastar.equals("loop") || mirnastar.equals("error")) {
                vo.error = true;
                return vo;
            }
            String os = System.getProperty("os.name").toLowerCase();
            String cmd = "";
            if (os.startsWith("win")) {
                cmd = "RNAcofold.exe --noPS -a";
            } else {
                cmd = Main.rnafoldlinux + "RNAcofold --noPS -a";
            }

            Process p = Runtime.getRuntime().exec(cmd);
            new PrintWriter(p.getOutputStream(), true).println(mirna + "&" + mirnastar);
            BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
            new PrintWriter(p.getOutputStream(), true).println("@");
            p.waitFor();
            String line = "";
            while (bri.ready()) {
                line = bri.readLine();
                output += line + "\n";
            }
            bri.close();
            p.destroy();
            vo.parseRNAcofold(output);
            return vo;
        } catch (Exception e) {
            vo.error = true;
            return vo;
        }
    }

    public static ViennaObject GetInfosDuplexRNAcofoldConstraint(String mirna, String mirnastar, String constraint) {
        ViennaObject vo = new ViennaObject();
        vo.calculateConstraint = false;
        String output = "";
        try {
            if (mirnastar.equals("loop") || mirnastar.equals("error")) {
                vo.error = true;
                return vo;
            }
            String os = System.getProperty("os.name").toLowerCase();
            String cmd = "";
            if (os.startsWith("win")) {
                cmd = "RNAcofold.exe --noPS -a";
            } else {
                cmd = Main.rnafoldlinux + "RNAcofold --noPS -C -a";
            }

            Process p = Runtime.getRuntime().exec(cmd);
            new PrintWriter(p.getOutputStream(), true).println(mirna + "&" + mirnastar);
            new PrintWriter(p.getOutputStream(), true).println(constraint);
            BufferedReader bri = new BufferedReader(new InputStreamReader(p.getInputStream()));
            new PrintWriter(p.getOutputStream(), true).println("@");
            p.waitFor();
            String line = "";
            while (bri.ready()) {
                line = bri.readLine();
                output += line + "\n";
            }
            bri.close();
            p.destroy();
            vo.parseRNAcofold(output);
            return vo;
        } catch (Exception e) {
            vo.error = true;
            return vo;
        }
    }

    public static String GetMfeDuplexTmp(String mirna, String mirnastar) {
        String mfe = "";
        File f = null;
        String out = "";
        try {

            if (mirnastar.equals("loop") || mirnastar.equals("error")) {
                return "0";
            }
            String os = System.getProperty("os.name").toLowerCase();
            if (os.startsWith("win")) {
                f = new File(/*c:/*/"tmpfolddup" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
                PrintWriter pw = new PrintWriter(new FileWriter(f));
                pw.write(mirna + "\n" + mirnastar);
                pw.close();
                //RNAfold.GetMfeDuplex(f);
                out = Tools.executeWindowsCommand("RNAduplex.exe < " + f.getAbsolutePath()).trim();
                mfe = out.substring(out.indexOf(" (") + 2, out.length() - 1);
            } else {
                f = new File("tmpfolddup" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
                PrintWriter pw = new PrintWriter(new FileWriter(f));
                pw.write(mirna + "\n" + mirnastar);
                pw.close();
                ///Users/mickael/Downloads/ViennaRNA-1.8.5/Progs/
                out = Tools.executeLinuxCommandGetLine(Main.rnafoldlinux + "RNAduplex < " + f.getAbsolutePath());
                mfe = out.substring(out.indexOf(" (") + 2, out.length() - 1);
            }
            if (Double.valueOf(mfe) > 0) {
                mfe = "0";
            }
            f.delete();
            return mfe;
        } catch (Exception e) {
            f.delete();
            //System.err.println("Error with RNAduplexe: "+out);
            //System.exit(0);
            return "0";
        }
    }

    public static boolean checkViennaTools() {
        //RNAfold

        if (GetSecondaryStructure("AAAAAAAAAAAUUUUUUUUUUUUUUUUAAAAAAAAAAAAAAUUUUUUUUUUUUUUAAAAAAAA").equals("0")) {
            System.err.println("checking RNAfold... FAIL");
            return false;
        }
        System.out.println("checking RNAfold... OK");
        //RNAcofold
//        if (GetInfosDuplexRNAcofold("AAAAAAAAAAAAAAAA","UUUUUUUUUUUUUUUU").hasError()){
//            System.err.println("checking RNAcofold... FAIL");
//            return false;
//        }
        //RNAduplex
        if (GetMfeDuplexRNAduplex("AAAAAAAAAAAAAAAA", "UUUUUUUUUUUUUUUU").equals("0")) {
            System.err.println("checking RNAduplex... FAIL");
            return false;
        }
        System.out.println("checking RNAduplex... OK");
        return true;
    }
}
