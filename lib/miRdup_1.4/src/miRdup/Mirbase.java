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
 * This class connect to miRbase to download mature miRNAs and their precursors.
 * Only choosen species based on a keyword are kept to train the model.
 * Two fasta files are created, one for the matures sequences, another one for precursors.
 * Data is also stored in an ArrayList of miRNAs objects
 */
package miRdup;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;

/**
 *
 * @author Mickael Leclercq
 */
public class Mirbase {

    private boolean debug = true;
    private static String outfileMirbaseMatures = "mirbase.matures.fasta";
    private static String outfileMirbasePrecursors = "mirbase.hairpins.fasta";
    private static String outfileMirbaseEMBL = "miRNA.dat";

    public void mirbase() {

    }

    public ArrayList getSequencesFromMirbaseEMBL(String keyword) {
        System.out.println("Downloading mirbase...");
        getEMBL();

        // get species list and add them to a hashmap
        if (Main.debug) {
            System.out.println("adding species in memory...");
        }
        if (keyword == null || keyword.isEmpty()) {
            keyword = "all";
        }
        HashMap<String, String> hmsp = getSpecies(keyword);

        // Add precursors to Hashmap
        if (Main.debug) {
            System.out.println("adding precursors in memory...");
        }
        HashMap<String, String> hmPrec = new HashMap<String, String>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(outfileMirbaseEMBL));
            String line = "";
            line = br.readLine();
            while (br.ready()) {
                if (line.startsWith("ID")) {
                    String name = line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                    String seq = "";
                    while (!line.startsWith("SQ")) {
                        line = br.readLine();
                    }
                    line = br.readLine();
                    while (!line.startsWith("//")) {
                        seq += line.substring(0, 75).replace(" ", "").toUpperCase();
                        line = br.readLine();
                    }
                    hmPrec.put(name.toLowerCase(), seq);
                } else {
                    line = br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // add mirnas in objects Arraylist
        ArrayList<MirnaObject> alobj = new ArrayList<MirnaObject>();
        int cpt = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(outfileMirbaseEMBL));
            String line = "";

            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith("ID")) {
                    String species = line.substring(5, 8);
                    if (hmsp.containsKey(species)) {
                        String shortname = line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                        while (!line.startsWith("FT")) {
                            line = br.readLine();
                        }
                        while (!line.startsWith("XX")) {
                            if (line.contains("FT   miRNA")) {
                                MirnaObject m = new MirnaObject();
                                m.setId(cpt++);
                                m.setShortName(shortname);
                                //mirna
                                String positions = line.substring(21);
                                int start = Integer.valueOf(positions.split("\\.")[0]) - 1;
                                int end = Integer.valueOf(positions.split("\\.")[2]);
                                String prec = hmPrec.get(m.getShortName());

                                m.setPrecursorSequence(prec);
                                String mirna = prec.substring(start, end);
                                m.setMatureSequence(mirna);

                                br.readLine(); // accession
                                line = br.readLine();// product
                                String product = line.substring(line.indexOf("=") + 1).replace("\"", "");
                                m.setFullName(product); // set star in the same time
                                line = br.readLine();// evidence
                                boolean experimental = true;
                                if (line.contains("not")) {
                                    experimental = false;
                                }
                                m.setExperimental(experimental);

                                alobj.add(m);
                            }
                            line = br.readLine();
                        }
                    }
                }
            }
            System.out.println("Total iterations (miRNAs with their corresponding precursor): " + alobj.size());
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        foldSequences(alobj);

        // print sequences that will be trained
        try {
            System.out.println("Store requested sequences in " + keyword + ".txt");
            PrintWriter pw = new PrintWriter(new FileWriter(keyword + ".txt"));
            PrintWriter pw2 = new PrintWriter(new FileWriter(keyword + ".notExperimental.txt"));
            int exp = 0;
            int notexp = 0;
            for (MirnaObject o : alobj) {
                if (o.isExperimental()/*&&!o.isStar()*/) {
                    pw.println(o.toStringTXT());
                    exp++;
                } else {
                    pw2.println(o.toStringTXT());
                    notexp++;
                }
            }
            pw.flush();
            pw.close();
            pw2.flush();
            pw2.close();
            System.out.println("Total experimental sequences: " + exp + "\n"
                    + "Total non experimental sequences "
                    + "(discarded in " + keyword + ".notExperimental.txt): " + notexp);
        } catch (Exception e) {
            e.printStackTrace();
        }

        //removing non experimental
        ArrayList<MirnaObject> alexpobj = new ArrayList<MirnaObject>();
        for (MirnaObject o : alobj) {
            if (o.isExperimental()/*&&!o.isStar()*/) {
                alexpobj.add(o);
            }
        }

        return alexpobj;
    }

    public ArrayList getSequencesFromMirbaseEMBLFile(String keyword, String emblFile, String structures, String organisms) {
        if (organisms.trim().isEmpty()) {
            getOrganisms();
            organisms = "organisms.txt";
        }

        System.out.println("Get informations from submitted files...:\n"
                + "\tEMBL File:" + emblFile + "\n"
                + "\tOrganisms:" + organisms + "\n"
                + "\tStructures of hairpins:" + structures + "\n");

        // get species list and add them to a hashmap
        if (Main.debug) {
            System.out.println("adding organisms in memory...");
        }
        HashMap<String, String> hmsp = new HashMap<String, String>();
        if (keyword == null || keyword.isEmpty()) {
            keyword = "all";
            hmsp = getSpecies(keyword);
        } else {
            hmsp = getSpecies(keyword, organisms);
        }

        // Add precursors to Hashmap
        if (Main.debug) {
            System.out.println("adding precursors in memory...");
        }
        HashMap<String, String> hmPrec = new HashMap<String, String>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(emblFile));
            String line = "";
            line = br.readLine();
            while (br.ready()) {
                if (line.startsWith("ID")) {
                    String name = line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                    String seq = "";
                    while (!line.startsWith("SQ")) {
                        line = br.readLine();
                    }
                    line = br.readLine();
                    while (!line.startsWith("//")) {
                        seq += line.substring(0, 75).replace(" ", "").toUpperCase();
                        line = br.readLine();
                    }
                    hmPrec.put(name.toLowerCase(), seq);
                } else {
                    line = br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // add mirnas in objects Arraylist
        ArrayList<MirnaObject> alobj = new ArrayList<MirnaObject>();
        int cpt = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(emblFile));
            String line = "";

            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith("ID")) {
                    String species = line.substring(5, 8);
                    if (hmsp.containsKey(species)) {
                        String shortname = line.substring(4, line.indexOf("stand")).trim().toLowerCase();
                        while (!line.startsWith("FT")) {
                            line = br.readLine();
                        }

                        while (!line.startsWith("XX")) {
                            if (line.substring(2).trim().startsWith("miRNA")) {
                                MirnaObject m = new MirnaObject();
                                m.setId(cpt++);
                                m.setShortName(shortname);
                                //mirna
                                String positions = line.substring(21);
                                int start = Integer.valueOf(positions.split("\\.")[0]) - 1;
                                int end = Integer.valueOf(positions.split("\\.")[2]);
                                String prec = hmPrec.get(m.getShortName());

                                m.setPrecursorSequence(prec);
                                String mirna = prec.substring(start, end);
                                m.setMatureSequence(mirna);

                                br.readLine(); // accession
                                line = br.readLine();// product
                                String product = line.substring(line.indexOf("=") + 1).replace("\"", "");
                                m.setFullName(product);
                                line = br.readLine();// evidence
                                boolean experimental = true;
                                if (line.contains("not")) {
                                    experimental = false;
                                }
                                m.setExperimental(experimental);

                                alobj.add(m);
                            }
                            line = br.readLine();
                        }
                    }
                }
            }
            System.out.println("Total iterations: " + alobj.size());
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.err.println(cpt);
        }

        // Get structures
        // must be an output of RNAfold
        if (structures == null || structures.isEmpty()) {
            foldSequences(alobj);
        } else {
            foldSequences(alobj, structures);
        }

        // print sequences that will be trained
        try {
            System.out.println("Print requested sequences in " + keyword + ".txt");
            PrintWriter pw = new PrintWriter(new FileWriter(keyword + ".txt"));
            PrintWriter pw2 = new PrintWriter(new FileWriter(keyword + ".notExperimental.txt"));
            int exp = 0;
            int notexp = 0;
            for (MirnaObject o : alobj) {
                if (o.isExperimental()/*&&!o.isStar()*/) {
                    pw.println(o.toStringTXT());
                    exp++;
                } else {
                    pw2.println(o.toStringTXT());
                    notexp++;
                }
            }
            pw.flush();
            pw.close();
            pw2.flush();
            pw2.close();
            System.out.println("Total experimental sequences: " + exp + ""
                    + "\nTotal non experimental sequences "
                    + "(discarded in " + keyword + ".notExperimental.txt): " + notexp);

        } catch (Exception e) {
            e.printStackTrace();
        }

        //removing non experimental and stars
        ArrayList<MirnaObject> alexpobj = new ArrayList<MirnaObject>();
        for (MirnaObject o : alobj) {
            if (o.isExperimental()/*&&!o.isStar()*/) {
                alexpobj.add(o);
            }
        }

        return alexpobj;
    }

    /**
     * Get sequences from miRbase based on a keyword
     *
     * @param keyword
     */
    public ArrayList getSequencesFromMirbase(String keyword) {
        System.out.println("Downloading mirbase...");
        getMatures();
        getPrecursors();

        // get species list and add them to a hashmap
        if (Main.debug) {
            System.out.println("adding species in memory...");
        }
        if (keyword == null || keyword.isEmpty()) {
            keyword = "all";
        }
        HashMap<String, String> hmsp = getSpecies(keyword);

        // Add precursors to Hashmap
        if (Main.debug) {
            System.out.println("adding precursors in memory...");
        }
        HashMap<String, String> hmPrec = new HashMap<String, String>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(outfileMirbasePrecursors));
            String line = "";
            line = br.readLine();
            while (br.ready()) {
                if (line.startsWith(">")) {
                    String name = line.substring(1, line.indexOf(" ")).trim();
                    String seq = "";
                    line = br.readLine();
                    while (line != null && !line.startsWith(">")) {
                        seq += line;
                        line = br.readLine();
                    }
                    hmPrec.put(name.toLowerCase(), seq);
                } else {
                    line = br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // add mirnas in objects Arraylist
        ArrayList<MirnaObject> alobj = new ArrayList<MirnaObject>();
        int cpt = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(outfileMirbaseMatures));
            String line = "";

            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith(">")) {
                    if (hmsp.containsKey(line.substring(1, line.indexOf("-")))) {
                        MirnaObject m = new MirnaObject();
                        m.setId(cpt++);
                        m.setFullName(line.substring(1, line.indexOf(" ")).toLowerCase());

                        String shortname = line.split(" ")[0].replaceAll("\\.[1-9]", "")
                                .replaceAll("\\.[1-9][1-9]", "")
                                .replaceAll("-[1-9] ", "")
                                .replaceAll("-[1-9][1-9] ", "")
                                .replace("-3p", "")
                                .replace("-5p", "")
                                .replace("*", "")
                                .replace(">", "").toLowerCase();
                        m.setShortName(shortname);
                        line = br.readLine();
                        m.setMatureSequence(line);

                        String precSeqs = hmPrec.get(shortname);
                        if (precSeqs == null) {
                            for (String s : hmPrec.keySet()) {
                                if (s.contains(shortname) && hmPrec.get(s).contains(m.getMatureSequence())) {
                                    alobj.add(new MirnaObject(m.getFullName(),
                                            m.getShortName(), m.getMatureSequence(), hmPrec.get(s)));
                                }
                            }
                        } else {
                            m.setPrecursorSequence(hmPrec.get(shortname));
                            if (m.getPrecursorSequence().contains(m.getMatureSequence())) {
                                alobj.add(m);
                            }
                        }
                    }
                }
            }
            System.out.println("Total iterations: " + alobj.size());
        } catch (Exception e) {
            e.printStackTrace();
        }

        foldSequences(alobj);

        try {
            System.out.println("Print requested sequences in " + keyword + ".txt");
            PrintWriter pw = new PrintWriter(new FileWriter(keyword + ".txt"));
            for (MirnaObject o : alobj) {
                pw.println(o.toStringTXT());
            }
            pw.flush();
            pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return alobj;
    }

    /**
     * Get sequences from miRbase based on a keyword
     *
     * @param keyword
     */
    public ArrayList getSequencesFromFiles(String matures, String hairpins, String organisms, String structures, String keyword, String filename) {
        if (organisms.trim().isEmpty() && !keyword.trim().isEmpty()) {
            getOrganisms();
            organisms = "organisms.txt";
        }

        if (keyword.trim().isEmpty()) {
            System.out.println("Get informations from submitted files...:\n"
                    + "\tMatures:" + matures + "\n"
                    + "\tHairpins:" + hairpins + "\n"
                    + "\tOrganisms:No keyword given, we assume sequences come from another source "
                    + "than miRbase\n"
                    + "\tStructures of hairpins:" + structures + "\n");
        } else {
            System.out.println("Get informations from submitted files...:\n"
                    + "\tMatures:" + matures + "\n"
                    + "\tHairpins:" + hairpins + "\n"
                    + "\tOrganisms:" + organisms + "\n"
                    + "\tStructures of hairpins:" + structures + "\n");
        }

        // get species list and add them to a hashmap
        HashMap<String, String> hmsp = null;
        if (!keyword.trim().isEmpty()) {
            if (Main.debug) {
                System.out.println("adding organisms in memory...");
            }
            hmsp = new HashMap<String, String>();
            if (keyword.isEmpty()) {
                keyword = "all";
                hmsp = getSpecies(keyword);
            } else {
                hmsp = getSpecies(keyword, organisms);
            }
        }

        // Add precursors to Hashmap
        if (Main.debug) {
            System.out.println("adding hairpins in memory...");
        }
        HashMap<String, String> hmPrec = new HashMap<String, String>();

        try {
            BufferedReader br = new BufferedReader(new FileReader(hairpins));
            String line = "";
            line = br.readLine();
            while (br.ready()) {
                if (line != null && line.startsWith(">")) {
                    String name = null;
                    try {
                        name = line.substring(1, line.indexOf(" ")).trim();
                    } catch (Exception e) {
                        name = line.substring(1).trim();
                    }
                    String seq = "";
                    line = br.readLine();
                    while (line != null && !line.startsWith(">")) {
                        seq += line;
                        line = br.readLine();
                    }
                    hmPrec.put(name.toLowerCase(), seq);
                } else {
                    line = br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

        // add mirnas in objects Arraylist
        ArrayList<MirnaObject> alobj = new ArrayList<MirnaObject>();
        int cpt = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(matures));
            String line = "";
            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith(">")) {
                    if (keyword.trim().isEmpty() || (hmsp != null && hmsp.containsKey(line.substring(1, line.indexOf("-"))))) {
                        MirnaObject m = new MirnaObject();
                        m.setId(cpt++);
                        try {
                            m.setFullName(line.substring(1, line.indexOf(" ")).toLowerCase());
                        } catch (Exception e) {
                            m.setFullName(line.substring(1));
                        }

                        String shortname = line.split(" ")[0].replaceAll("\\.[1-9]", "")
                                .replaceAll("\\.[1-9][1-9]", "")
                                .replaceAll("-[1-9] ", "")
                                .replaceAll("-[1-9][1-9] ", "")
                                .replace("-3p", "")
                                .replace("-5p", "")
                                .replace("*", "")
                                .replace(">", "").toLowerCase();
                        m.setShortName(shortname);
                        line = br.readLine();
                        m.setMatureSequence(line);

                        String precSeqs = hmPrec.get(shortname);
                        if (precSeqs == null) {
                            for (String s : hmPrec.keySet()) {
                                if (s.contains(shortname) && hmPrec.get(s).contains(m.getMatureSequence())) {
                                    alobj.add(new MirnaObject(m.getFullName(),
                                            m.getShortName(), m.getMatureSequence(), hmPrec.get(s)));
                                }
                            }
                        } else {
                            m.setPrecursorSequence(hmPrec.get(shortname));
                            if (m.getPrecursorSequence().contains(m.getMatureSequence())) {
                                alobj.add(m);
                            }
                        }
                    }
                }
            }
            System.out.println("Total iterations: " + alobj.size());
        } catch (Exception e) {
            e.printStackTrace();
        }

        // Get structures
        // must be an output of RNAfold
        if (structures == null || structures.isEmpty()) {
            foldSequences(alobj);
        } else {
            foldSequences(alobj, structures);
        }

        try {
            System.out.println("Print requested sequences in " + filename + ".txt");
            PrintWriter pw = new PrintWriter(new FileWriter(filename + ".txt"));
            for (MirnaObject o : alobj) {
                pw.println(o.toStringTXT());
            }
            pw.flush();
            pw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return alobj;
    }

    /**
     * Fold in batch with rnafold
     *
     * @param alobj
     * @return
     */
    public static ArrayList foldSequences(ArrayList<MirnaObject> alobj) {
        System.out.println("Folding sequences...");
        String os = System.getProperty("os.name").toLowerCase();
        File tmpFile = null;
        if (os.startsWith("win")) {
            tmpFile = new File(/*c:/*/"tmpfold" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
        } else {
            tmpFile = new File("tmpfold" + (Math.abs((int) (Calendar.getInstance().hashCode() / (new Date()).getTime() + Math.random() * 10000000) * 10000)));
        }

        try {
            //save sequences in a file
            PrintWriter pw = new PrintWriter(new FileWriter(tmpFile));
            for (MirnaObject o : alobj) {
                pw.println(">" + o.getPrecursorSequence() + "\n" + o.getPrecursorSequence());
            }
            pw.flush();
            pw.close();

            //read folded sequences
            File folded = Vienna.GetSecondaryStructuresFromFile(tmpFile);
            BufferedReader br = new BufferedReader(new FileReader(folded));
            String line = "";
            HashMap<String, String> hm = new HashMap<String, String>();
            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith(">")) {
                    String seq = br.readLine();
                    String struct = br.readLine().trim();
                    String mfe = struct.substring(struct.indexOf(" (") + 2, struct.length() - 1);
                    struct = struct.replace("(" + mfe + ")", "");
                    hm.put(seq, struct);
                }
            }

            // replace sequences in arraylist
            for (MirnaObject o : alobj) {
                o.setStructure(hm.get(o.getPrecursorSequence()));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return alobj;
    }

    /**
     * Get structures from a file
     *
     * @param alobj
     * @param structures
     * @return
     */
    public static ArrayList foldSequences(ArrayList<MirnaObject> alobj, String structures) {
        System.out.println("Folding sequences...");
        try {
            BufferedReader br = new BufferedReader(new FileReader(structures));
            String line = "";
            HashMap<String, String> hm = new HashMap<String, String>();
            while (br.ready()) {
                line = br.readLine();
                if (line.startsWith(">")) {
                    String seq = br.readLine();
                    String struct = br.readLine().trim();
                    if (!struct.contains("(")) {
                        Tools.error("The file " + structures + " does not contains a good structure for sequence " + seq);
                    }
                    String mfe = struct.substring(struct.indexOf(" (") + 2, struct.length() - 1);
                    struct = struct.replace("(" + mfe + ")", "");
                    hm.put(seq, struct);
                }
            }

            // replace sequences in arraylist
            for (MirnaObject o : alobj) {
                o.setStructure(hm.get(o.getPrecursorSequence()));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return alobj;
    }

    /**
     * Get species informations depending a keyword If keyword=all, all mirbase
     * is used
     *
     * @param keyword
     */
    public static HashMap getSpecies(String keyword) {
        if (Main.debug) {
            System.out.println("Downloading species...");
        }
        System.out.println("Requested keyword: " + keyword);
        HashMap<String, String> hmsp = new HashMap<String, String>();
        keyword = keyword.toLowerCase();
        getOrganisms();
        ArrayList<String> lines = new ArrayList<String>();
        try {
            BufferedReader br = new BufferedReader(new FileReader("organisms.txt"));
            while (br.ready()) {
                lines.add(br.readLine());
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        System.out.println("Included species: ");
        if (keyword.equals("all")) {
            for (String line : lines) {
                System.out.println(line);
                hmsp.put(line.split("\t")[0], "");
            }
        } else {
            for (String line : lines) {
                if (line.toLowerCase().contains(keyword)) {
                    System.out.println(line);
                    hmsp.put(line.split("\t")[0], "");
                }
            }
        }
        return hmsp;
    }

    /**
     * Get species informations depending a keyword If keyword=all, all mirbase
     * is used
     *
     * @param keyword
     */
    public static HashMap getSpecies(String keyword, String organimsFile) {
        if (Main.debug) {
            System.out.println("Get species...");
        }
        System.out.println("Requested keyword: " + keyword);
        HashMap<String, String> hmsp = new HashMap<String, String>();
        keyword = keyword.toLowerCase();

        System.out.println("Included species: ");
        try {
            BufferedReader br = new BufferedReader(new FileReader(organimsFile));
            String line = "";
            while (br.ready()) {
                line = br.readLine();
                if (keyword.equals("all")) {
                    System.out.println(line);
                    hmsp.put(line.split("\t")[0], "");
                } else if (line.toLowerCase().contains(keyword)) {
                    System.out.println(line);
                    hmsp.put(line.split("\t")[0], "");
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        return hmsp;
    }

    /**
     * get mature miRNAs from miRbas
     */
    public static void getMatures() {
        if (Main.debug) {
            System.out.println("Downloading mature miRNAs...");
        }
        Tools.downloadFile("ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.zip");
        Tools.unzipFile("mature.fa.zip");
        Tools.deleteFile("mature.fa.zip");
        Tools.renameFile("mature.fa", outfileMirbaseMatures);

    }

    /**
     * Get precursors file from miRbase
     */
    public static void getPrecursors() {
        if (Main.debug) {
            System.out.println("Downloading miRNAs precursors...");
        }
        Tools.downloadFile("ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.zip");
        Tools.unzipFile("hairpin.fa.zip");
        Tools.deleteFile("hairpin.fa.zip");
        Tools.renameFile("hairpin.fa", outfileMirbasePrecursors);
    }

    /**
     * get mature miRNAs from miRbas
     */
    public static void getEMBL() {
        if (Main.debug) {
            System.out.println("Downloading mature miRNAs...");
        }
        Tools.downloadFile("ftp://mirbase.org/pub/mirbase/CURRENT/miRNA.dat.zip");
        Tools.unzipFile("miRNA.dat.zip");
        Tools.deleteFile("miRNA.dat.zip");
        Tools.renameFile("miRNA.dat", outfileMirbaseEMBL);

    }

    /**
     * get organism list from miRbas
     */
    public static void getOrganisms() {
        if (Main.debug) {
            System.out.println("Downloading mature miRNAs...");
        }
        Tools.downloadFile("ftp://mirbase.org/pub/mirbase/CURRENT/organisms.txt.zip");
        Tools.unzipFile("organisms.txt.zip");
        Tools.deleteFile("organisms.txt.zip");
    }

}
