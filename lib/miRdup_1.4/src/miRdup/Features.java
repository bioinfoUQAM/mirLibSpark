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
 * Get features from mirnas, precursors and structures
 */
package miRdup;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;

/**
 *
 * @author Mickael Leclercq
 */
public final class Features {

    String mirna;
    String prec;
    String precStruc;
    String mirnaStruc;
    private String complementaritySequence;
    String mirnastarStruc;
    Boolean positive;
    public static DecimalFormat df = new DecimalFormat();

    public Features(String miRNA, String precursor, String secondaryStructure, Boolean positiveDataset) {
        mirna = miRNA;
        prec = precursor;
        precStruc = secondaryStructure;
        positive = positiveDataset;
        mirnaStruc = precStruc.substring(prec.indexOf(mirna), prec.indexOf(mirna) + mirna.length());
        computeComplementaritySequence();
        mirnastarStruc = precStruc.substring(prec.indexOf(complementaritySequence), prec.indexOf(complementaritySequence) + complementaritySequence.length());
        df.setMaximumFractionDigits(2);
    }

    Features() {

    }

    //useless since RNAcofold calculate it
    public String getMFE() {
        computeComplementaritySequence();
        return Vienna.GetMfeDuplexRNAduplex(mirna, getComplementaritySequence());
    }

    //partition function and base pairing probabilities
    public String getCoFoldInfos() {
        computeComplementaritySequence();
        ViennaObject vo = Vienna.GetInfosDuplexRNAcofold(mirna, getComplementaritySequence());
        if (vo.hasError()) {
            return vo.toStringError();
        } else {
            return vo.toStringRNAcofold();
        }
    }

    public int getLength() {
        return mirna.length();
    }

    //Maximum length on the miRNA without bulges
    public double getGCperc() {
        int l = mirna.length();
        int gc = 0;
        for (int i = 0; i < l; i++) {
            if (mirna.charAt(i) == 'C' || mirna.charAt(i) == 'G') {
                gc++;
            }
        }
        double perc = (gc * 100) / l;
        return perc;
    }

    //Maximum length on the miRNA without bulges
    public double getGCpercNormalized() {
        int l = mirna.length();
        int gcmirna = 0;
        for (int i = 0; i < l; i++) {
            if (mirna.charAt(i) == 'C' || mirna.charAt(i) == 'G') {
                gcmirna++;
            }
        }
        double percMirna = (gcmirna * 100) / l;

        l = prec.length();
        int gcprec = 0;
        for (int i = 0; i < l; i++) {
            if (prec.charAt(i) == 'C' || prec.charAt(i) == 'G') {
                gcprec++;
            }
        }
        double percPrec = (gcprec * 100) / l;

        double norm = percMirna / percPrec;

        return norm;
    }

    //Maximum length on the miRNA without bulges
    public int getMaximumLengthWithoutBulges() {
        int lenght = 0;
        int max = 0;
        for (int i = 0; i < mirnaStruc.length(); i++) {
            if (mirnaStruc.charAt(i) == '(' || mirnaStruc.charAt(i) == ')') {
                lenght++;
            } else {
                lenght = 0;
            }
            if (max < lenght) {
                max = lenght;
            }
        }
        return max;
    }

    //Maximum length on the miRNA without bulges in percentage
    public double getMaximumLengthWithoutBulgesPerc() {
        int b = getMaximumLengthWithoutBulges();
        double perc = (b * 100) / mirna.length();
        return perc;
    }

    //Start Length without bulges
    public int getStartLengthWithoutBulges() {
        int cpt = 0;
        try {
            if (mirnaStruc.replace("(", "").trim().length() == 0
                    || mirnaStruc.replace(")", "").trim().length() == 0) {
                return mirnaStruc.length();
            } else {

                char c = mirnaStruc.charAt(0);
                int i = 0;
                while (c == '(' || c == ')') {
                    c = mirnaStruc.charAt(i);
                    cpt++;
                    i++;
                }
                return cpt - 1;
            }
        } catch (Exception e) {
            return cpt - 1;
        }
    }

    //Base pairs in duplexe miRNA-miRNA*
    public int getBasePairsInDuplex() {
        int p = 0;
        for (int i = 0; i < mirna.length(); i++) {
            if (mirnaStruc.charAt(i) == '(' || mirnaStruc.charAt(i) == ')') {
                p++;
            }
        }
        return p;
    }

    //20mer base paired
    public boolean getPresenceOfPerfect20MerBasePair() {
        if (mirnaStruc.contains("((((((((((((((((((((") || mirnaStruc.contains("))))))))))))))))))))")) {
            return true;
        } else {
            return false;
        }
    }

    // start of 20mer base paired
    public int getStartOfPerfect20MerBasePair() {
        if (getPresenceOfPerfect20MerBasePair()) {
            int s = mirnaStruc.indexOf("((((((((((((((((((((");
            if (s == -1) {
                s = mirnaStruc.indexOf("))))))))))))))))))))");
            }
            return s + 1;
        } else {
            return -1;
        }
    }

    //10mer base paired
    public boolean getPresenceOfPerfect10MerBasePair() {
        if (mirnaStruc.contains("((((((((((") || mirnaStruc.contains("))))))))))")) {
            return true;
        } else {
            return false;
        }
    }

    // start of 10mer base paired
    public int getStartOfPerfect10MerBasePair() {
        if (getPresenceOfPerfect10MerBasePair()) {
            int s = mirnaStruc.indexOf("((((((((((");
            if (s == -1) {
                s = mirnaStruc.indexOf(")))))))))))");
            }
            return s + 1;
        } else {
            return -1;
        }
    }

    //5mer base paired
    public boolean getPresenceOfPerfect5MerBasePair() {
        if (mirnaStruc.contains("(((((") || mirnaStruc.contains(")))))")) {
            return true;
        } else {
            return false;
        }
    }

    // start of 5mer base paired
    public int getStartOfPerfect5MerBasePair() {
        if (getPresenceOfPerfect5MerBasePair()) {
            int s = mirnaStruc.indexOf("(((((");
            if (s == -1) {
                s = mirnaStruc.indexOf(")))))");
            }
            return s + 1;
        } else {
            return -1;
        }
    }

    //Start of perfect 8mer base paired followed by a 95% 12mer base paired
    public int getStartOfPerfect8And12MerBasePair() {

        return 0;
    }

    //Presence of A
    public boolean getPresenceOfA() {
        if (mirna.contains("A")) {
            return true;
        } else {
            return false;
        }
    }

    //Percentage of A
    public double getPercOfA() {
        int nt = 0;
        for (int i = 0; i < mirna.length(); i++) {
            if (mirna.charAt(i) == 'A') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    //Presence of U
    public boolean getPresenceOfU() {
        if (mirna.contains("U")) {
            return true;
        } else {
            return false;
        }
    }

    //Percentage of U
    public double getPercOfU() {
        int nt = 0;
        for (int i = 0; i < mirna.length(); i++) {
            if (mirna.charAt(i) == 'U') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    //Presence of G
    public boolean getPresenceOfG() {
        if (mirna.contains("G")) {
            return true;
        } else {
            return false;
        }
    }

    //Percentage of G
    public double getPercOfG() {
        int nt = 0;
        for (int i = 0; i < mirna.length(); i++) {
            if (mirna.charAt(i) == 'G') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    //Presence of C
    public boolean getPresenceOfC() {
        if (mirna.contains("C")) {
            return true;
        } else {
            return false;
        }
    }

    //Percentage of C
    public double getPercOfC() {
        int nt = 0;
        for (int i = 0; i < mirna.length(); i++) {
            if (mirna.charAt(i) == 'C') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    //Distance from terminal loop
    public int getDistanceFromTerminalLoop() {
        String arm = getArm();
        try {
            if (arm.equals("5'")) {
                int mirnaloc = prec.indexOf(mirna);
                int closestLoopEndPar = precStruc.substring(mirnaloc).indexOf(")") + mirnaloc;
                int loopStart = precStruc.substring(mirnaloc, closestLoopEndPar).lastIndexOf("(") + mirnaloc;
                return loopStart - (prec.indexOf(mirna) + mirna.length() - 1);
            } else if (arm.equals("3'")) {
                int mirnaloc = prec.indexOf(mirna);
                int closestLoopOpenPar = precStruc.substring(0, mirnaloc).lastIndexOf("(");
                int loopEnd = precStruc.substring(closestLoopOpenPar, mirnaloc).indexOf(")") + closestLoopOpenPar;
                return mirnaloc - loopEnd;
            } else {//case where arm=loop or error. Return a negative value
                int inloop = 0;
                try {
                    inloop = -mirnaStruc.substring(mirnaStruc.lastIndexOf("("), mirnaStruc.indexOf(")")).length() + 1;
                } catch (Exception e) {
                    inloop = -1;
                }
                return inloop;
            }
        } catch (Exception e) {
            return -1;
        }

    }

    //Distance from terminal loop. Old
    public int getDistanceFromTerminalLoopOld() {
        int loopStart = precStruc.lastIndexOf("(");
        int loopEnd = precStruc.indexOf(")");
        String arm = getArm();
        if (arm.equals("5'")) {
            return loopStart - (prec.indexOf(mirna) + mirna.length());
        } else if (arm.equals("3'")) {
            return (prec.indexOf(mirna) - loopEnd);
        } else { //case where arm=loop or error. Return a negative value

            if (prec.indexOf(mirna) + mirna.length() > loopStart) {
                return loopStart - (prec.indexOf(mirna) + mirna.length());
            } else if (prec.indexOf(mirna) < loopEnd) {
                return (prec.indexOf(mirna) - loopEnd);
            }
            return 0;
        }

    }

    //Distance from start of the hairpin
    public int getDistanceFromHairpinStart() {
        String arm = getArm();
        if (arm.equals("5'")) {
            return prec.indexOf(mirna);
        } else if (arm.equals("3'")) {
            return prec.length() - (prec.indexOf(mirna) + mirna.length());
        } else {
            // Case where arm=loop or error
            if (!mirnaStruc.contains(")")) {
                return prec.indexOf(mirna);
            } else if (!mirnaStruc.contains("(")) {
                return prec.length() - (prec.indexOf(mirna) + mirna.length());
            } else {
                return prec.indexOf(mirna);
            }
        }
    }

    //Length overlap in loop region
    public int getLengthOfOverlapInLoop() {
        int tl = getDistanceFromTerminalLoop();
        if (tl < 0) {
            return Math.abs(tl);
        } else {
            return 0;
        }
    }

    //Average number of paired bases in window of 7
    public double getAverageNumberOfPairedBasesInWindow7() {
        ArrayList<Integer> al = new ArrayList<Integer>();
        for (int i = 0; i <= mirna.length() - 7; i++) {
            String s = mirnaStruc.substring(i, i + 7);
            int pair = 0;
            for (int j = 0; j < s.length(); j++) {
                if (s.charAt(j) == '(' || s.charAt(j) == ')') {
                    pair++;
                }
            }
            al.add(pair);
        }
        double tot = 0;
        for (Integer i : al) {
            tot = tot + i;
        }
        return (tot / al.size());
    }

    //Average number of paired bases in window of 5
    public double getAverageNumberOfPairedBasesInWindow5() {
        ArrayList<Integer> al = new ArrayList<Integer>();
        for (int i = 0; i <= mirna.length() - 5; i++) {
            String s = mirnaStruc.substring(i, i + 5);
            int pair = 0;
            for (int j = 0; j < s.length(); j++) {
                if (s.charAt(j) == '(' || s.charAt(j) == ')') {
                    pair++;
                }
            }
            al.add(pair);
        }
        double tot = 0;
        for (Integer i : al) {
            tot = tot + i;
        }
        return (tot / al.size());
    }

    //Average number of paired bases in window of 3
    public double getAverageNumberOfPairedBasesInWindow3() {
        ArrayList<Integer> al = new ArrayList<Integer>();
        for (int i = 0; i <= mirna.length() - 3; i++) {
            String s = mirnaStruc.substring(i, i + 3);
            int pair = 0;
            for (int j = 0; j < s.length(); j++) {
                if (s.charAt(j) == '(' || s.charAt(j) == ')') {
                    pair++;
                }
            }
            al.add(pair);
        }
        double tot = 0;
        for (Integer i : al) {
            tot = tot + i;
        }
        return (tot / al.size());
    }

    //Bulge at position 0
    public boolean getBulgeAtPosition0() {
        if (mirnaStruc.charAt(0) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at position -1
    public boolean getBulgeAtPositionMinus1() {
        try {
            int mirnaPos = prec.indexOf(mirna);
            if (precStruc.charAt(mirnaPos - 1) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) {
            return false;
        }
    }

    //Bulge at position 1
    public boolean getBulgeAtPosition1() {
        if (mirnaStruc.charAt(1) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at position -2
    public boolean getBulgeAtPositionMinus2() {
        try {
            int mirnaPos = prec.indexOf(mirna);
            if (precStruc.charAt(mirnaPos - 2) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) { //happening when mirna start at the begining
            return false;
        }
    }

    //Bulge at position 2
    public boolean getBulgeAtPosition2() {
        if (mirnaStruc.charAt(2) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at position -3
    public boolean getBulgeAtPositionMinus3() {
        try {
            int mirnaPos = prec.indexOf(mirna);
            if (precStruc.charAt(mirnaPos - 3) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) { //happening when mirna start at the begining
            return false;
        }
    }

    //Bulge at position 3
    public boolean getBulgeAtPosition3() {
        if (mirnaStruc.charAt(3) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at position -4
    public boolean getBulgeAtPositionMinus4() {
        try {
            int mirnaPos = prec.indexOf(mirna);
            if (precStruc.charAt(mirnaPos - 4) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) { //happening when mirna start at the begining
            return false;
        }
    }

    //Bulge at end position 0
    public boolean getBulgeAtEndPosition0() {
        if (mirnaStruc.charAt(mirnaStruc.length() - 1) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at end position -1
    public boolean getBulgeAtEndPositionMinus1() {
        if (mirnaStruc.charAt(mirnaStruc.length() - 2) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at end position -2
    public boolean getBulgeAtEndPositionMinus2() {
        if (mirnaStruc.charAt(mirnaStruc.length() - 3) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at end position -3
    public boolean getBulgeAtEndPositionMinus3() {
        if (mirnaStruc.charAt(mirnaStruc.length() - 4) == '.') {
            return true;
        } else {
            return false;
        }
    }

    //Bulge at end position +1
    public boolean getBulgeAtEndPositionPlus1() {
        try {
            int mirnaPos = prec.indexOf(mirna) + mirna.length();
            if (precStruc.charAt(mirnaPos) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) {
            return false;
        }
    }

    //Bulge at end position +2
    public boolean getBulgeAtEndPositionPlus2() {
        try {
            int mirnaPos = prec.indexOf(mirna) + mirna.length();
            if (precStruc.charAt(mirnaPos + 1) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) {
            return false;
        }
    }

    //Bulge at end position +1
    public boolean getBulgeAtEndPositionPlus3() {
        try {
            int mirnaPos = prec.indexOf(mirna) + mirna.length();
            if (precStruc.charAt(mirnaPos + 2) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) {
            return false;
        }
    }

    //Bulge at end position +1
    public boolean getBulgeAtEndPositionPlus4() {
        try {
            int mirnaPos = prec.indexOf(mirna) + mirna.length();
            if (precStruc.charAt(mirnaPos + 3) == '.') {
                return true;
            } else {
                return false;
            }
        } catch (Exception e) {
            return false;
        }
    }

    //Number of bulges
    public int getNumberOfBulges() {
        int bulges = 0;
        for (int i = 0; i <= mirnaStruc.length() - 2; i++) {
            String s = mirnaStruc.substring(i, i + 2);
            if (s.equals(".)") || s.equals(".(")) {
                bulges++;
            }
        }
        return bulges;
    }

    //length in the biggest bulge
    public int getLenghtOfBiggestBulge() {
        int l = 0;
        int max = 0;
        for (int i = 0; i < mirnaStruc.length(); i++) {
            if (mirnaStruc.charAt(i) == '.') {
                l++;
            } else {
                l = 0;
            }
            if (max < l) {
                max = l;
            }
        }
        return max;
    }

    //length on the miRNA biggest bulge in percentage
    public double getLengthBiggestBulgesPerc() {
        int b = getLenghtOfBiggestBulge();
        double perc = (b * 100) / mirna.length();
        return perc;
    }

    // Pairs AU, GC, GU
    public String getPercOfbasepairs() {

        try {
            if (getMirnaIncludedInLoop()) {
                return "0.0,0.0,0.0";
            }

            String sequence;
            String struc;
            int mirnaStart;
            int mirnaEnd;
            int starStart;
            int starEnd;
            if (getArm().equals("5'")) {
                sequence = mirna + getComplementaritySequence();
                struc = mirnaStruc + mirnastarStruc;
                mirnaStart = 0;
                mirnaEnd = mirna.length() - 1;
                starStart = mirna.length();
                starEnd = sequence.length();
            } else {
                sequence = getComplementaritySequence() + mirna;
                struc = mirnastarStruc + mirnaStruc;
                mirnaStart = 0;
                mirnaEnd = getComplementaritySequence().length() - 1;
                starStart = getComplementaritySequence().length();
                starEnd = sequence.length();
            }

            HashMap<String, Integer> hm = new HashMap<String, Integer>();
            hm.put("AU", 0);
            hm.put("GC", 0);
            hm.put("GU", 0);
            String pairs[] = {"AU", "GC", "GU"};

            int posOnStar = starEnd;
            for (int i = mirnaStart; i <= mirnaEnd; i++) {
                posOnStar--;
                char a = struc.charAt(i);
                char b = struc.charAt(posOnStar);
//            System.out.println(i+"\t"+a+" "+sequence.charAt(i));
//            System.out.println(posOnStar+"\t"+b+" "+sequence.charAt(posOnStar));
                if (isParenthese(a) && isParenthese(b)) {
                    char c = sequence.charAt(i);
                    char d = sequence.charAt(posOnStar);
                    String pair = String.valueOf(c + "" + d);
                    if (pair.equals("AU") || pair.equals("UA")) {
                        int tmp = hm.get("AU");
                        tmp++;
                        hm.put("AU", tmp);
                    }
                    if (pair.equals("GC") || pair.equals("CG")) {
                        int tmp = hm.get("GC");
                        tmp++;
                        hm.put("GC", tmp);
                    }
                    if (pair.equals("GU") || pair.equals("UG")) {
                        int tmp = hm.get("GU");
                        tmp++;
                        hm.put("GU", tmp);
                    }
                } else if (a == '.' && isParenthese(b)) {
                    posOnStar++;
                } else if (b == '.' && isParenthese(a)) {
                    i--;
                }
            }

            int parenthesis = 0;
            for (int i = 0; i < mirnaStruc.length(); i++) {
                if (mirnaStruc.charAt(i) == '(' || mirnaStruc.charAt(i) == ')') {
                    parenthesis++;
                }
            }

            String pairsvalues = "";
            for (String p : pairs) {
                if (hm.get(p) == null) {
                    pairsvalues += ",0.0";
                } else {
                    double a = hm.get(p);
                    double l = parenthesis;
                    if (l == 0.0) {
                        return "0.0,0.0,0.0";
                    }
                    double d = (a / l) * 100;
                    pairsvalues += "," + df.format(d).replace(",", ".");
                }
            }
            return pairsvalues.substring(1);
        } catch (Exception e) {
//            return "0.0,0.0,0.0";
            return "?,?,?";
        }
    }

    // Pairs AU, GC, GU norm
    public String getPercOfbasepairsNormalized() {
        String basePairs = getPercOfbasepairs();
        if (basePairs.equals("0.0,0.0,0.0")) {
            return "0.0,0.0,0.0";
        }
        try {
            double AUnom = Double.valueOf(basePairs.split(",")[0]);
            double GCnom = Double.valueOf(basePairs.split(",")[0]);
            double GUnom = Double.valueOf(basePairs.split(",")[0]);

            double AUden = ((getNucleotideNumber(prec, 'A') * getNucleotideNumber(prec, 'U')) / prec.length()) * mirna.length();
            double GCden = ((getNucleotideNumber(prec, 'G') * getNucleotideNumber(prec, 'C')) / prec.length()) * mirna.length();
            double GUden = ((getNucleotideNumber(prec, 'G') * getNucleotideNumber(prec, 'U')) / prec.length()) * mirna.length();

            double AUnorm = AUnom / AUden;
            double GCnorm = GCnom / GCden;
            double GUnorm = GUnom / GUden;

            return df.format(AUnorm).replace(",", ".") + ","
                    + df.format(GCnorm).replace(",", ".") + ","
                    + df.format(GUnorm).replace(",", ".");
        } catch (Exception e) {
            return "?,?,?";
        }

    }

    public double getNucleotideNumber(String sequence, char nt) {
        double count = 0;
        for (int i = 0; i < sequence.length(); i++) {
            if (sequence.charAt(i) == nt) {
                count++;
            }
        }
        return count;
    }

    // Pairs AU, GC
    public String getPercOfbasedpairsBestFeatures() {

        try {
            if (getMirnaIncludedInLoop()) {
                return "0.0,0.0";
            }

            String sequence;
            String struc;
            int mirnaStart;
            int mirnaEnd;
            int starStart;
            int starEnd;
            if (getArm().equals("5'")) {
                sequence = mirna + getComplementaritySequence();
                struc = mirnaStruc + mirnastarStruc;
                mirnaStart = 0;
                mirnaEnd = mirna.length() - 1;
                starStart = mirna.length();
                starEnd = sequence.length();
            } else {
                sequence = getComplementaritySequence() + mirna;
                struc = mirnastarStruc + mirnaStruc;
                mirnaStart = 0;
                mirnaEnd = getComplementaritySequence().length() - 1;
                starStart = getComplementaritySequence().length();
                starEnd = sequence.length();
            }

            HashMap<String, Integer> hm = new HashMap<String, Integer>();
            hm.put("AU", 0);
            hm.put("GC", 0);
            String pairs[] = {"AU", "GC"};

            int posOnStar = starEnd;
            for (int i = mirnaStart; i <= mirnaEnd; i++) {
                posOnStar--;
                char a = struc.charAt(i);
                char b = struc.charAt(posOnStar);
//            System.out.println(i+"\t"+a+" "+sequence.charAt(i));
//            System.out.println(posOnStar+"\t"+b+" "+sequence.charAt(posOnStar));
                if (isParenthese(a) && isParenthese(b)) {
                    char c = sequence.charAt(i);
                    char d = sequence.charAt(posOnStar);
                    String pair = String.valueOf(c + "" + d);
                    if (pair.equals("AU") || pair.equals("UA")) {
                        int tmp = hm.get("AU");
                        tmp++;
                        hm.put("AU", tmp);
                    }
                    if (pair.equals("GC") || pair.equals("CG")) {
                        int tmp = hm.get("GC");
                        tmp++;
                        hm.put("GC", tmp);
                    }
                } else if (a == '.' && isParenthese(b)) {
                    posOnStar++;
                } else if (b == '.' && isParenthese(a)) {
                    i--;
                }
            }

            int parenthesis = 0;
            for (int i = 0; i < mirnaStruc.length(); i++) {
                if (mirnaStruc.charAt(i) == '(' || mirnaStruc.charAt(i) == ')') {
                    parenthesis++;
                }
            }

            String pairsvalues = "";
            for (String p : pairs) {
                if (hm.get(p) == null) {
                    pairsvalues += ",0.0";
                } else {
                    double a = hm.get(p);
                    double l = parenthesis;
                    if (l == 0.0) {
                        return "0.0,0.0";
                    }
                    double d = (a / l) * 100;
                    pairsvalues += "," + df.format(d).replace(",", ".");
                }
            }
            return pairsvalues.substring(1);
        } catch (Exception e) {
//            return "0.0,0.0,0.0";
            return "?,?";
        }

    }

    // dinucleotides
    // Pairs A[ACGU]
    public double getPercOfAA() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'A' && mirna.charAt(i + 1) == 'A') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfAC() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'A' && mirna.charAt(i + 1) == 'C') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfAG() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'A' && mirna.charAt(i + 1) == 'G') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfAU() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'A' && mirna.charAt(i + 1) == 'U') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    // Pairs C[ACGU]
    public double getPercOfCA() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'C' && mirna.charAt(i + 1) == 'A') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfCC() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'C' && mirna.charAt(i + 1) == 'C') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfCG() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'C' && mirna.charAt(i + 1) == 'G') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfCU() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'C' && mirna.charAt(i + 1) == 'U') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    // Pairs G[ACGU]
    public double getPercOfGA() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'G' && mirna.charAt(i + 1) == 'A') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfGC() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'G' && mirna.charAt(i + 1) == 'C') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfGG() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'G' && mirna.charAt(i + 1) == 'G') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfGU() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'G' && mirna.charAt(i + 1) == 'U') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    // Pairs U[ACGU]
    public double getPercOfUA() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'U' && mirna.charAt(i + 1) == 'A') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfUC() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'U' && mirna.charAt(i + 1) == 'C') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfUG() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'U' && mirna.charAt(i + 1) == 'G') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    public double getPercOfUU() {
        int nt = 0;
        for (int i = 0; i < mirna.length() - 1; i++) {
            if (mirna.charAt(i) == 'U' && mirna.charAt(i + 1) == 'U') {
                nt++;
            }
        }
        double perc = (nt * 100) / mirna.length();
        return perc;
    }

    // get complementarity region //not working for every cases....
    public void getComplementaritySequenceOld() {
        if (getMirnaIncludedInLoop()) {
            setComplementaritySequence("loop");
        }

        String star = "";
        int taille = mirna.length();
        int mirnaStart = prec.indexOf(mirna);
        int mirnaEnd = mirnaStart + taille;
        int loopStart = precStruc.lastIndexOf("(");
        int loopEnd = precStruc.indexOf(")");

        String arm = null;
        try {
            if (precStruc.substring(mirnaStart, mirnaEnd).contains("(")) {
                arm = "5'";
            } else {
                arm = "3'";
            }
        } catch (Exception e) {
            //System.err.println("error at "+mirna);
            setComplementaritySequence("error");
        }

        /////////5'
        if (arm.equals("5'")) {

            // Distance du mirna a la loop effecuté en comptant le nombre de parentheses
            int nbrParentheses = 0;
            String intervale = null;
            try {
                intervale = precStruc.substring(mirnaEnd, loopStart);
            } catch (Exception e) {
                intervale = "";
            }
            for (char c : intervale.toCharArray()) {
                if (c == '(') {
                    nbrParentheses++;
                }
            }
            // Determination du départ du mirna star en fonction du nombre de parentheses
            // après la fin de la loop
            int cpt = 0;
            int starStart = loopEnd;
            while (cpt != nbrParentheses) {
                if (precStruc.charAt(starStart) == ')') {
                    cpt++;
                    starStart++;
                } else {
                    starStart++;
                }
            }

            //recupération du mirna star
            int posOnStar = starStart;
            boolean dec = false;
            for (int i = mirnaEnd; i > mirnaStart; i--) {
                char a = precStruc.charAt(i);
                char b = precStruc.charAt(posOnStar++);
                if (isParenthese(a) && isParenthese(b) || a == b) {
                    star += prec.charAt(posOnStar);
                } else if (a == '.' && isParenthese(b)) {
                    i--;
                    posOnStar--;
                    dec = true;
                } else if (b == '.' && isParenthese(a)) {
                    star += prec.charAt(posOnStar);
                }
            }
            if (dec) {
                star += prec.charAt(posOnStar - 1);
            }
        } ///////////////3'
        else {

            // Distance du mirna a la loop effecuté en comptant le nombre de parentheses
            int nbrParentheses = 0;
            String intervale = null;
            try {
                intervale = precStruc.substring(loopEnd, mirnaStart);
            } catch (Exception e) {
                intervale = "";
            }
            for (char c : intervale.toCharArray()) {
                if (c == ')') {
                    nbrParentheses++;
                }
            }

            // Determination du départ du mirna star en fonction du nombre de parentheses
            // avant le début de la loop
            int cpt = 0;
            int starEnd = loopStart;
            while (cpt != nbrParentheses) {
                if (precStruc.charAt(starEnd) == '(') {
                    cpt++;
                    starEnd--;
                } else {
                    starEnd--;
                }
            }

            //recupération du mirna star
            int posOnStar = starEnd;
            boolean dec = false;
            for (int i = mirnaStart; i < mirnaEnd; i++) {
                if (posOnStar >= 0) {
                    char a = precStruc.charAt(i);
                    char b = precStruc.charAt(posOnStar);
                    if (isParenthese(a) && isParenthese(b) || a == b) {
                        star = prec.charAt(posOnStar) + star;
                        posOnStar--;
                    } else if (a == '.' && isParenthese(b)) {
                        i++;
                        dec = true;
                    } else if (b == '.' && isParenthese(a)) {
                        star = prec.charAt(posOnStar) + star;
                        i--;
                        posOnStar--;
                    }
                }
            }
            if (dec) {
                star = prec.charAt(posOnStar) + star;
            }
        }
        //System.out.println(star);
        setComplementaritySequence(star);
    }

    // get complementarity region
    public void computeComplementaritySequence() {
        if (getMirnaIncludedInLoop()) {
            setComplementaritySequence("loop");
        }

        int size = mirna.length();
        int mirnaStart = prec.indexOf(mirna);
        int mirnaEnd = prec.indexOf(mirna) + size;
        String arm = getArm();

        /////////5'
        if (arm.equals("5'")) {
            //analysing before mirna if presence of loop
            String startToMirna = precStruc.substring(0, mirnaStart);
            // get open parenthesis minus close parenthesis between precursor start and mirna
            int openPar = 0;
            for (int i = 0; i < startToMirna.length(); i++) {
                if (startToMirna.charAt(i) == '(') {
                    openPar++;
                }
                if (startToMirna.charAt(i) == ')') {
                    openPar--;
                }
            }
            //analysing from the end of precursor if presence of loop

            int openParEnd = 0;
            int i = precStruc.length() - 1;
            while (openPar != openParEnd) {
                if (precStruc.charAt(i) == ')') {
                    openParEnd++;
                }
                if (precStruc.charAt(i) == '(') {
                    openParEnd--;
                }
                if (openPar != openParEnd) {
                    i--;
                }

            }

            // get complement
            String star = "";
            int posOnStar = i;
            for (i = mirnaStart; i < mirnaEnd; i++) {
                posOnStar--;
                char a = precStruc.charAt(i);
                char b = precStruc.charAt(posOnStar);
                if (isParenthese(a) && isParenthese(b) || a == b) {
                    star += prec.charAt(posOnStar);
                } else if (a == '.' && isParenthese(b)) {
                    posOnStar++;
                } else if (b == '.' && isParenthese(a)) {
                    star += prec.charAt(posOnStar);
                    i--;
                }
            }
            setComplementaritySequence(Tools.Reverse(star));

        } ///////////////3'
        else {
            String EndToMirna = precStruc.substring(mirnaStart + size, precStruc.length());
            // get open parenthesis minus close parenthesis between precursor end and mirna end
            int closePar = 0;
            for (int i = 0; i < EndToMirna.length(); i++) {
                if (EndToMirna.charAt(i) == ')') {
                    closePar++;
                }
                if (EndToMirna.charAt(i) == '(') {
                    closePar--;
                }
            }

            //analysing from the end of precursor if presence of loop
            int openPar = 0;
            int i = 0;
            while (openPar != closePar) {
                if (precStruc.charAt(i) == '(') {
                    openPar++;
                }
                if (precStruc.charAt(i) == ')') {
                    openPar--;
                }
                if (openPar != closePar) {
                    i++;
                }
            }

            // get complement
            String star = "";
            int posOnStar = i - 1;
            for (i = mirnaEnd - 1; i >= mirnaStart; i--) {
                posOnStar++;
                char a = precStruc.charAt(i);
                prec.charAt(posOnStar);
                char b = precStruc.charAt(posOnStar);
                if (isParenthese(a) && isParenthese(b) || a == b) {
                    star += prec.charAt(posOnStar);
                } else if (a == '.' && isParenthese(b)) {
                    posOnStar--;
                } else if (b == '.' && isParenthese(a)) {
                    star += prec.charAt(posOnStar);
                    i++;
                }

            }
            setComplementaritySequence(star);
        }
    }

    public String getMirnaStar() {
        int csStart = prec.indexOf(getComplementaritySequence());
        int csEnd = prec.indexOf(getComplementaritySequence()) + getComplementaritySequence().length();

        try {
            csStart = csStart + 2;
            csEnd = csEnd + 2;
            String star = prec.substring(csStart, csEnd);
            return star;
        } catch (Exception e) {
            String star = getComplementaritySequence();
            return star;
        }
    }

    public String getArm() {
        int taille = mirna.length();
        int mirnaStart = prec.indexOf(mirna);
        int mirnaEnd = mirnaStart + taille;
        if (getMirnaIncludedInLoop()) {
            return "loop";
        }
        String arm = null;
        try {
            if (precStruc.substring(mirnaStart, mirnaEnd).contains("(")) {
                arm = "5'";
            } else {
                arm = "3'";
            }
        } catch (Exception e) {
            return "error";
        }
        return arm;
    }

    public boolean isParenthese(char a) {
        if (a == '(' || a == ')') {
            return true;
        } else {
            return false;
        }
    }

    /**
     * Vérifie que le miRNA n'est pas complètement inclu dans la loop
     *
     * @param prec
     * @param struct
     * @param mirna
     * @return true si le miRNA est dans la loop
     */
    public boolean getMirnaIncludedInLoop() {
        try {
            int start = prec.indexOf(mirna);
            int end = start + mirna.length();
            if (end > prec.length() || start == -1) {
                return false;
            }
            if (precStruc.substring(start, end).contains("(") && precStruc.substring(start, end).contains(")")) {
                return true;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return false;
    }

    public boolean isMultipleLoop() {
        int firstclosePar = precStruc.indexOf(")");
        if (precStruc.substring(firstclosePar).contains("(")) {
            return true;
        } else {
            return false;
        }
    }

    //mipred
    public String getAllTriplets() {
        String triplets[] = {
            "A...",
            "C...",
            "G...",
            "U...",
            "A(..",
            "C(..",
            "G(..",
            "U(..",
            "A((.",
            "C((.",
            "G((.",
            "U((.",
            "A.((",
            "C.((",
            "G.((",
            "U.((",
            "A(((",
            "C(((",
            "G(((",
            "U(((",
            "A.(.",
            "C.(.",
            "G.(.",
            "U.(.",
            "A..(",
            "C..(",
            "G..(",
            "U..(",
            "A(.(",
            "C(.(",
            "G(.(",
            "U(.(",};

        HashMap<String, Integer> hm = new HashMap<String, Integer>();
        for (int i = 1; i < mirna.length() - 1; i++) {
            String triplet = String.valueOf(mirna.charAt(i));
            triplet += mirnaStruc.substring(i - 1, i + 2);
            if (hm.get(triplet) == null) {
                hm.put(triplet, 1);
            } else {
                int tmp = hm.get(triplet);
                tmp++;
                hm.put(triplet, tmp);
            }
        }
        String tripletsvalues = "";
        for (String t : triplets) {
            if (hm.get(t) == null) {
                tripletsvalues += ",0.0";
            } else {
                double a = hm.get(t);
                double l = mirna.length();
                double d = (a / l) * 100;
                tripletsvalues += "," + df.format(d).replace(",", ".");
            }
        }

        return tripletsvalues.substring(1);

    }

    public String getBestTriplets() {
        String triplets[] = {
            "A...",
            "C...",
            "G...",
            "U...",
            "G(((",};

        HashMap<String, Integer> hm = new HashMap<String, Integer>();
        for (int i = 1; i < mirna.length() - 1; i++) {
            String triplet = String.valueOf(mirna.charAt(i));
            triplet += mirnaStruc.substring(i - 1, i + 2);
            if (hm.get(triplet) == null) {
                hm.put(triplet, 1);
            } else {
                int tmp = hm.get(triplet);
                tmp++;
                hm.put(triplet, tmp);
            }
        }
        String tripletsvalues = "";
        for (String t : triplets) {
            if (hm.get(t) == null) {
                tripletsvalues += ",0.0";
            } else {
                double a = hm.get(t);
                double l = mirna.length();
                double d = (a / l) * 100;
                tripletsvalues += "," + df.format(d).replace(",", ".");
            }
        }

        return tripletsvalues.substring(1);

    }

    private Character control(char c) {
        if (c == 'A' || c == 'U' || c == 'G' || c == 'C' || c == 'N') {
            return c;
        } else {
            return 'N';
        }
    }

    public Character getNtAtStart0() {
        return control(mirna.charAt(0));
    }

    public Character getNtAtStartPlus1() {
        return control(mirna.charAt(1));
    }

    public Character getNtAtStartPlus2() {
        return control(mirna.charAt(2));
    }

    public Character getNtAtStartPlus3() {
        return control(mirna.charAt(3));
    }

    public Character getNtAtStartMinus1() {
        int mirnaStart = prec.indexOf(mirna);
        try {
            return control(prec.charAt(mirnaStart - 1));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtStartMinus2() {
        int mirnaStart = prec.indexOf(mirna);
        try {
            return control(prec.charAt(mirnaStart - 2));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtStartMinus3() {
        int mirnaStart = prec.indexOf(mirna);
        try {
            return control(prec.charAt(mirnaStart - 3));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtEnd0() {
        return control(mirna.charAt(mirna.length() - 1));
    }

    public Character getNtAtEndPlus1() {
        int mirnaEnd = prec.indexOf(mirna) + mirna.length();
        try {
            return control(prec.charAt(mirnaEnd));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtEndPlus2() {
        int mirnaEnd = prec.indexOf(mirna) + mirna.length() + 1;
        try {
            return control(prec.charAt(mirnaEnd));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtEndPlus3() {
        int mirnaEnd = prec.indexOf(mirna) + mirna.length() + 2;
        try {
            return control(prec.charAt(mirnaEnd));
        } catch (Exception e) {
            return control('?');
        }
    }

    public Character getNtAtEndMinus1() {
        return control(mirna.charAt(mirna.length() - 2));
    }

    public Character getNtAtEndMinus2() {
        return control(mirna.charAt(mirna.length() - 3));
    }

    public Character getNtAtEndMinus3() {
        return control(mirna.charAt(mirna.length() - 4));
    }

    public static String getID() {
        int i = AdaptDataForWeka.id++;
        return "ID" + i;
    }

    /**
     * @return the complementaritySequence
     */
    public String getComplementaritySequence() {
        return complementaritySequence;
    }

    /**
     * @param complementaritySequence the complementaritySequence to set
     */
    public void setComplementaritySequence(String complementaritySequence) {
        this.complementaritySequence = complementaritySequence;
    }

    public String toStringError() {
        //127
        //String features=getID()+",0,0,0,0,0,0,0,false,0,false,0,false,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,false,0,0.0,0.0,0.0,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0.0,0.0,0.0,0.0,0.0,N,N,N,N,N,N,N,N,N,N,N,N,N,N,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,false";
        //100
        String features = getID() + ",0,0,0,0,0,0,0,false,0,false,0,false,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,false,0,0.0,0.0,0.0,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,false,0,0,0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0,0.0,0.0,N,N,N,N,N,N,true";

        return features;
    }

    public String toStringAllAttributes() {
        String features = "";
        try {
            features = ""
                    + getID() + ","
                    + getLength() + ","
                    + getMFE() + ","
                    + df.format(getGCperc()).replace(",", ".") + ","
                    //                    + df.format(getGCpercNormalized()).replace(",", ".") + ","
                    + getMaximumLengthWithoutBulges() + ","
                    + df.format(getMaximumLengthWithoutBulgesPerc()).replace(",", ".") + ","
                    + getStartLengthWithoutBulges() + ","
                    + getBasePairsInDuplex() + ","
                    + getPresenceOfPerfect20MerBasePair() + ","
                    + getStartOfPerfect20MerBasePair() + ","
                    + getPresenceOfPerfect10MerBasePair() + ","
                    + getStartOfPerfect10MerBasePair() + ","
                    + getPresenceOfPerfect5MerBasePair() + ","
                    + getStartOfPerfect5MerBasePair() + ","
                    + df.format(getPercOfA()).replace(",", ".") + ","
                    + df.format(getPercOfU()).replace(",", ".") + ","
                    + df.format(getPercOfG()).replace(",", ".") + ","
                    + df.format(getPercOfC()).replace(",", ".") + ","
                    + df.format(getPercOfAA()).replace(",", ".") + ","
                    + df.format(getPercOfUA()).replace(",", ".") + ","
                    + df.format(getPercOfGA()).replace(",", ".") + ","
                    + df.format(getPercOfCA()).replace(",", ".") + ","
                    + df.format(getPercOfAU()).replace(",", ".") + ","
                    + df.format(getPercOfUU()).replace(",", ".") + ","
                    + df.format(getPercOfGU()).replace(",", ".") + ","
                    + df.format(getPercOfCU()).replace(",", ".") + ","
                    + df.format(getPercOfAG()).replace(",", ".") + ","
                    + df.format(getPercOfUG()).replace(",", ".") + ","
                    + df.format(getPercOfGG()).replace(",", ".") + ","
                    + df.format(getPercOfCG()).replace(",", ".") + ","
                    + df.format(getPercOfAC()).replace(",", ".") + ","
                    + df.format(getPercOfUC()).replace(",", ".") + ","
                    + df.format(getPercOfGC()).replace(",", ".") + ","
                    + df.format(getPercOfCC()).replace(",", ".") + ","
                    + getDistanceFromTerminalLoop() + ","
                    + getDistanceFromHairpinStart() + ","
                    + getMirnaIncludedInLoop() + ","
                    + getLengthOfOverlapInLoop() + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow7()).replace(",", ".") + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow5()).replace(",", ".") + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow3()).replace(",", ".") + ","
                    + getBulgeAtPosition0() + ","
                    + getBulgeAtPositionMinus1() + ","
                    + getBulgeAtPosition1() + ","
                    + getBulgeAtPositionMinus2() + ","
                    + getBulgeAtPosition2() + ","
                    + getBulgeAtPositionMinus3() + ","
                    + getBulgeAtPosition3() + ","
                    + getBulgeAtPositionMinus4() + ","
                    + getBulgeAtEndPosition0() + ","
                    + getBulgeAtEndPositionPlus1() + ","
                    + getBulgeAtEndPositionMinus1() + ","
                    + getBulgeAtEndPositionPlus2() + ","
                    + getBulgeAtEndPositionMinus2() + ","
                    + getBulgeAtEndPositionPlus3() + ","
                    + getBulgeAtEndPositionMinus3() + ","
                    + getBulgeAtEndPositionPlus4() + ","
                    + getNumberOfBulges() + ","
                    + getLenghtOfBiggestBulge() + ","
                    + df.format(getLengthBiggestBulgesPerc()).replace(",", ".") + ","
                    + getAllTriplets() + ","
                    + getPercOfbasepairs() + ","
                    //                    + getPercOfbasepairsNormalized()+","

                    + getNtAtStart0() + ","
                    + getNtAtStartMinus1() + ","
                    //                    + getNtAtStartMinus2()+ ","
                    //                    + getNtAtStartMinus3()+ ","
                    + getNtAtStartPlus1() + ","
                    //                    + getNtAtStartPlus2()+ ","
                    //                    + getNtAtStartPlus3()+ ","
                    + getNtAtEnd0() + ","
                    + getNtAtEndMinus1() + ","
                    //                    + getNtAtEndMinus2()+ ","
                    //                    + getNtAtEndMinus3()+ ","
                    + getNtAtEndPlus1() + ","
                    //                    + getNtAtEndPlus2()+ ","
                    //                    + getNtAtEndPlus3()+ ","
                    //                    + getCoFoldInfos()+","
                    + positive;
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(mirna + " " + prec + " " + precStruc);
        }
        return features.replace("\uFFFD", "?").replace("�", "");
    }

    public String toStringBestAttributes() {
        String features = "";
        try {
            features = ""
                    + getID() + ","
                    //                    + getLength() + ","
                    //                    + getMFE() + ","
                    //                    + df.format(getGCperc()).replace(",", ".") + ","
                    + getMaximumLengthWithoutBulges() + ","
                    + df.format(getMaximumLengthWithoutBulgesPerc()).replace(",", ".") + ","
                    //                    + getStartLengthWithoutBulges() + ","
                    + getBasePairsInDuplex() + ","
                    //                    + getPresenceOfPerfect20MerBasePair() + ","
                    //                    + getStartOfPerfect20MerBasePair() + ","
                    //                    + getPresenceOfPerfect10MerBasePair() + ","
                    //                    + getStartOfPerfect10MerBasePair() + ","
                    //                    + getPresenceOfPerfect5MerBasePair() + ","
                    + getStartOfPerfect5MerBasePair() + ","
                    //
                    //                    + df.format(getPercOfA()).replace(",", ".") + ","
                    //                    + df.format(getPercOfU()).replace(",", ".") + ","
                    //                    + df.format(getPercOfG()).replace(",", ".") + ","
                    //                    + df.format(getPercOfC()).replace(",", ".") + ","
                    //                    + df.format(getPercOfAA()).replace(",", ".") + ","
                    //                    + df.format(getPercOfUA()).replace(",", ".") + ","
                    //                    + df.format(getPercOfGA()).replace(",", ".") + ","
                    //                    + df.format(getPercOfCA()).replace(",", ".") + ","
                    //                    + df.format(getPercOfAU()).replace(",", ".") + ","
                    //                    + df.format(getPercOfUU()).replace(",", ".") + ","
                    //                    + df.format(getPercOfGU()).replace(",", ".") + ","
                    //                    + df.format(getPercOfCU()).replace(",", ".") + ","
                    //                    + df.format(getPercOfAG()).replace(",", ".") + ","
                    //                    + df.format(getPercOfUG()).replace(",", ".") + ","
                    //                    + df.format(getPercOfGG()).replace(",", ".") + ","
                    //                    + df.format(getPercOfCG()).replace(",", ".") + ","
                    //                    + df.format(getPercOfAC()).replace(",", ".") + ","
                    //                    + df.format(getPercOfUC()).replace(",", ".") + ","
                    //                    + df.format(getPercOfGC()).replace(",", ".") + ","
                    //                    + df.format(getPercOfCC()).replace(",", ".") + ","
                    //
                    + getDistanceFromTerminalLoop() + ","
                    + getDistanceFromHairpinStart() + ","
                    + getMirnaIncludedInLoop() + ","
                    + getLengthOfOverlapInLoop() + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow7()).replace(",", ".") + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow5()).replace(",", ".") + ","
                    + df.format(getAverageNumberOfPairedBasesInWindow3()).replace(",", ".") + ","
                    //                    + getBulgeAtPosition0() + ","
                    //                    + getBulgeAtPositionMinus1() + ","
                    //                    + getBulgeAtPosition1() + ","
                    //                    + getBulgeAtPositionMinus2() + ","
                    //                    + getBulgeAtPosition2() + ","
                    //                    + getBulgeAtPositionMinus3() + ","
                    //                    + getBulgeAtPosition3() + ","
                    //                    + getBulgeAtPositionMinus4() + ","
                    //
                    //                    + getBulgeAtEndPosition0() + ","
                    //                    + getBulgeAtEndPositionPlus1() + ","
                    //                    + getBulgeAtEndPositionMinus1() + ","
                    //                    + getBulgeAtEndPositionPlus2() + ","
                    //                    + getBulgeAtEndPositionMinus2() + ","
                    //                    + getBulgeAtEndPositionPlus3() + ","
                    //                    + getBulgeAtEndPositionMinus3() + ","
                    //                    + getBulgeAtEndPositionPlus4() + ","
                    //
                    //                    + getNumberOfBulges() + ","
                    + getLenghtOfBiggestBulge() + ","
                    + df.format(getLengthBiggestBulgesPerc()).replace(",", ".") + ","
                    + getBestTriplets() + ","
                    + getPercOfbasepairs() + ","
                    //
                    //                    + getNtAtStart0()+ ","
                    //                    + getNtAtStartMinus1()+ ","
                    ////                    + getNtAtStartMinus2()+ ","
                    ////                    + getNtAtStartMinus3()+ ","
                    //                    + getNtAtStartPlus1()+ ","
                    ////                    + getNtAtStartPlus2()+ ","
                    ////                    + getNtAtStartPlus3()+ ","
                    //                    + getNtAtEnd0()+ ","
                    //                    + getNtAtEndMinus1()+ ","
                    ////                    + getNtAtEndMinus2()+ ","
                    ////                    + getNtAtEndMinus3()+ ","
                    //                    + getNtAtEndPlus1()+ ","
                    ////                    + getNtAtEndPlus2()+ ","
                    ////                    + getNtAtEndPlus3()+ ","

                    + positive;
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(mirna + " " + prec + " " + precStruc);
        }
        return features;
    }

}
