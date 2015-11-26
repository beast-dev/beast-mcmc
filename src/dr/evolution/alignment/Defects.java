/*
 * Defects.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evolution.alignment;

import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.TaxonList;

import java.io.EOFException;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.TreeSet;

/**
 * @author Alexei Drummond
 *
 * @version $Id: Defects.java,v 1.6 2005/05/23 10:44:07 alexei Exp $
 */
public class Defects {

    private final ArrayList<Defect> defects = new ArrayList<Defect>();
    private final Set<Integer> defectiveSequences = new HashSet<Integer>();
    private final Set<Integer> defectiveSites = new HashSet<Integer>();
    int sequenceCount = 0;
    static final int STOP = -1;
    static final int INDEL = -2;
    int totalReads = 0;

    public Defects(Alignment alignment) {

        Codons codons = Codons.UNIVERSAL;
        Nucleotides nucs = Nucleotides.INSTANCE;

        Alignment codonAlignment = new ConvertAlignment(codons, alignment);

        sequenceCount = codonAlignment.getSequenceCount();

        for (int i = 0; i < codonAlignment.getSequenceCount(); i++) {
            for (int j = 0; j < codonAlignment.getSiteCount(); j++) {
                int state = codonAlignment.getState(i, j);

                if (codons.isStopCodon(state)) {
                    //System.out.print("*");
                    //    System.out.println("Found stop codon: " + triplet + " in sequence " + codonAlignment.getTaxonId(i) + " at site " + j);
                    addDefect(STOP, i, j);
                    totalReads += 1;
                } else if (codons.isGapState(state)) {
                    //System.out.print("?");
                    // go back to original alignment and check
                    int start = j*3;
                    int state1 = alignment.getState(i,start);
                    int state2 = alignment.getState(i,start+1);
                    int state3 = alignment.getState(i,start+2);
                    if (nucs.isGapState(state1) && nucs.isGapState(state2) && nucs.isGapState(state3)) {
                        // real gap
                    } else {
                        // defect
                        //    System.out.println("Found gap in sequence " + codonAlignment.getTaxonId(i) + " at site " + j);
                        addDefect(INDEL, i, j);
                        //totalReads += 1;
                    }

                } else {
                    //System.out.print(" ");
                    totalReads += 1;
                }
            }
            //System.out.println();
        }
    }

    private void addDefect(int code, int sequence, int site) {
        defects.add(new Defect(code, sequence, site));
        defectiveSequences.add(sequence);
        defectiveSites.add(site);
    }

    public int getDefectiveSequenceCount() {
        return defectiveSequences.size();
    }

    public int getDefectiveSiteCount() {
        return defectiveSites.size();
    }

    public int getDefectiveSites(int sequence) {
        int count = 0;
        for( Defect defect : defects ) {
            if( defect.getSequence() == sequence ) {
                count += 1;
            }
        }
        return count;
    }

    public int getStopSites(int sequence) {
        int count = 0;
        for( Defect defect : defects ) {
            if( defect.getSequence() == sequence && defect.isStop() ) {
                count += 1;
            }
        }
        return count;
    }

    public int getDefectiveSequences(int site) {
        int count = 0;
        for(Defect defect : defects) {
            if( defect.getSequence() == site ) {
                count += 1;
            }
        }
        return count;
    }

    /**
     * @return the total number of defects in the alignment
     */
    public int getDefectCount() {
        return defects.size();
    }

    public int getStopCodonCount() {
        int count = 0;
        for( Defect defect : defects ) {
            if( defect.isStop() ) {
                count += 1;
            }
        }
        return count;
    }

    public int getSequenceCount(int defectCount) {
        int count = 0;
        for (int i = 0; i < sequenceCount; i++) {
            if (getDefectiveSites(i) == defectCount) {
                count += 1;
            }
        }
        return count;
    }

    /**
     * @return the maximum number of defects found in any single sequence.
     */
    public int getMaxDefectiveSites() {

        int maxDefects = 0;
        for (int i = 0; i < sequenceCount; i++) {
            int defects = getDefectiveSites(i);
            if (defects > maxDefects) {
                maxDefects = defects;
            }
        }
        return maxDefects;
    }

    public int getMaxStopSites() {

        int maxStops = 0;
        for (int i = 0; i < sequenceCount; i++) {
            int stops = getStopSites(i);
            if (stops > maxStops) {
                maxStops = stops;
            }
        }
        return maxStops;
    }


    /**
     * @param defectCount
     * @return A set of integers representing the sequences with the given defect count.
     */
    public Set<Integer> getSequences(int defectCount) {
        TreeSet<Integer> set = new TreeSet<Integer>();
        for (int i = 0; i < sequenceCount; i++) {
            if (getDefectiveSites(i) == defectCount) {
                set.add(i);
            }
        }
        return set;
    }

    /**
     * @param stopCount
     * @return A set of integers representing the sequences with the given defect count.
     */
    public Set<Integer> getSequencesByStopCount(int stopCount) {
        TreeSet<Integer> set = new TreeSet<Integer>();
        for (int i = 0; i < sequenceCount; i++) {
            if (getStopSites(i) == stopCount) {
                set.add(i);
            }
        }
        return set;
    }


    public int getStopSequenceCount(int stopCount) {
        int count = 0;
        for (int i = 0; i < sequenceCount; i++) {
            if (getStopSites(i) == stopCount) {
                count += 1;
            }
        }
        return count;
    }

    public int getTotalCodonsMinusGaps() { return totalReads; }


    class Defect {

        private int type = 0;
        private int sequence = 0;
        private int site = 0;

        public Defect(int type, int sequence, int site) {
            this.type = type;
            this.sequence = sequence;
            this.site = site;
        }

        public boolean isStop() { return type == STOP; }

        public boolean isIndel() { return type == INDEL; }

        public int getSequence() { return sequence; }
        public int getSite() { return site; }
    }

    /**
     * @param fileName
     * @throws java.io.IOException
     */
    private static Alignment readNexusFile(String fileName) throws java.io.IOException {

        Alignment alignment = null;
        TaxonList taxonList = null;

        try {
            FileReader reader = new FileReader(fileName);

            NexusImporter importer = new NexusImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxonList != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxonList = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = importer.parseDataBlock(taxonList);
                        if (taxonList == null) {
                            taxonList = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        // ignore tree block
                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

        } catch (Importer.ImportException ime) {
            System.err.println("Error reading alignment: " + ime);
        }

        return alignment;
    }

    private static int[][] generateTransversionMutants(int[] triplet) {
        int[][] tvmutants = new int[6][3];

        for (int index = 0; index < 6; index += 1) {
            for (int pos = 0; pos < 3; pos++) {
                tvmutants[index][pos] = triplet[pos];
            }
        }

        int index = 0;
        for (int pos = 0; pos < 3; pos++) {
            if (triplet[pos] % 2 == 0) {
                tvmutants[index++][pos] = 1;
                tvmutants[index++][pos] = 3;
            } else {
                tvmutants[index++][pos] = 0;
                tvmutants[index++][pos] = 2;
            }
        }
        return tvmutants;
    }

    private static int[][] generateTransitionMutants(int[] triplet) {
        int[][] timutants = new int[3][3];

        for (int index = 0; index < 3; index += 1) {
            for (int pos = 0; pos < 3; pos++) {
                timutants[index][pos] = triplet[pos];
            }
        }

        for (int pos = 0; pos < 3; pos++) {
            switch (triplet[pos]) {
                case 0: timutants[pos][pos] = 2; break;
                case 1: timutants[pos][pos] = 3; break;
                case 2: timutants[pos][pos] = 0; break;
                case 3: timutants[pos][pos] = 1; break;
                default: throw new IllegalArgumentException();
            }
        }
        return timutants;
    }

    public static void main(String[] args) throws java.io.IOException {

        String alignmentFileName = null;
        if (args != null && args.length > 0) {
            alignmentFileName = args[0];
        }

        Alignment alignment = readNexusFile(alignmentFileName);

        Defects defects = new Defects(alignment);

        int siteCount = alignment.getSiteCount() / 3;
        int sequenceCount = alignment.getSequenceCount();
        int totalSites = siteCount * sequenceCount;
        int totalReads = defects.getTotalCodonsMinusGaps();

        int defectiveSites = defects.getDefectiveSiteCount();
        int defectiveSequences = defects.getDefectiveSequenceCount();
        int totalDefects = defects.getDefectCount();
        int totalStops = defects.getStopCodonCount();

        double siteP = (double)defectiveSites/(double)siteCount;
        double sequenceP = (double)defectiveSequences/(double)sequenceCount;
        double totalP = (double)totalDefects/(double)totalReads;
        double totalSP = (double)totalStops/(double)totalReads;

        System.out.println("Matrix size=" + totalSites);
        System.out.println("Non-gap codons=" + totalReads);

        System.out.println(defectiveSequences + "/" + sequenceCount + "(" + sequenceP + ") defective sequences.");
        System.out.println(defectiveSites + "/" + siteCount + "(" + siteP + ") defective sites.");
        System.out.println(totalDefects + "/" + totalReads + "(" + totalP + ") defects.");
        System.out.println(totalStops + "/" + totalReads + "(" + totalSP + ") stop codons.");

        double mean = (double)totalDefects/(double)sequenceCount;
        System.out.println(mean + " defects per sequence");

        double meanS = (double)totalStops/(double)sequenceCount;
        System.out.println(meanS + " stops per sequence");

        int maxDefects = defects.getMaxDefectiveSites();
        int maxStops = defects.getMaxStopSites();
        System.out.println("defective sequences:");
        //Set seqs = defects.getSequences(maxDefects);
        //for (Iterator i = seqs.iterator(); i.hasNext();) {
        //    Integer index = (Integer)i.next();
        //    String name = alignment.getTaxonId(index.intValue());
        //    System.out.println("  " + name);
        //}
        for (int d = 1; d <= maxDefects; d++) {
            Set<Integer> seqs = defects.getSequences(d);
            for(Integer seq : seqs) {
                String name = alignment.getTaxonId(seq);
                System.out.println(d + "  " + name);
            }

        }

        System.out.println("Defects\tSequences\texpected");
        double probTerm;
        probTerm = Math.exp(-mean); // probability of 0
        for (int i=0; i<10; i++) {

            System.out.println(i+"\t" + defects.getSequenceCount(i) + "\t" + probTerm*sequenceCount);

            // compute probability of n from prob. of n-1
            probTerm *= mean / (i+1);
        }

        System.out.println("Stops\tSequences\texpected");
        probTerm = Math.exp(-meanS); // probability of 0
        for (int i=0; i<10; i++) {

            System.out.println(i+"\t" + defects.getStopSequenceCount(i) + "\t" + probTerm*sequenceCount);

            // compute probability of n from prob. of n-1
            probTerm *= meanS / (i+1);
        }

        System.out.println("stop-codon sequences:");
        for (int d = 1; d <= maxStops; d++) {
            Set<Integer> seqs = defects.getSequencesByStopCount(d);
            for(Integer index : seqs) {
                String name = alignment.getTaxonId(index);
                System.out.println(d + "  " + name);
            }

        }


        Consensus con = new Consensus("mode", alignment, true);
        Sequence consensus = con.getConsensusSequence();
        SimpleAlignment a = new SimpleAlignment();
        a.addSequence(new Sequence(consensus));
        ConvertAlignment ca = new ConvertAlignment(Codons.UNIVERSAL, a);

        double[] pstop = new double[ca.getSiteCount()];
        int[] counts = new int[10];
        for (double kappa = 2.0; kappa <= 2.0; kappa += 1.0) {
            for (int i = 0; i < ca.getSiteCount(); i++) {
                int state = ca.getState(0,i);
                if (Codons.UNIVERSAL.isStopCodon(state)) {
                    throw new RuntimeException("Consensus has a stop codon in it at position " + i + "!");
                }
                int[] triplet = Codons.UNIVERSAL.getTripletStates(state);
                int[][] tvmutants = generateTransversionMutants(triplet);
                int[][] timutants = generateTransitionMutants(triplet);

                int tvStops = 0;
                for (int j = 0; j < 6; j++) {
                    if (Codons.UNIVERSAL.isStopCodon(Codons.UNIVERSAL.getState(tvmutants[j][0],tvmutants[j][1],tvmutants[j][2]))) {
                        tvStops += 1;
                    }
                }

                int tiStops = 0;
                for (int j = 0; j < 3; j++) {
                    if (Codons.UNIVERSAL.isStopCodon(Codons.UNIVERSAL.getState(timutants[j][0],timutants[j][1],timutants[j][2]))) {
                        tiStops += 1;
                    }
                }

                pstop[i] = (tvStops + kappa*tiStops) / (6.0 + kappa*3.0);
                counts[tvStops+tiStops] += 1;
            }
            System.out.println("kappa = " + kappa + " pstop=" + dr.stats.DiscreteStatistics.mean(pstop));
        }

        System.out.println("stop-mutations\tcodons");
        for (int i = 0; i < 10; i++) {
            System.out.println(i+"\t" + counts[i]);
        }

        /*int count = 0;
        int reps = 100000;
        for (int i = 0; i < reps; i++) {
            SimpleAlignment simpleAlignment = new SimpleAlignment();
            simpleAlignment.addSequence(new Sequence(consensus));

            //pick a random position
            int pos = random.nextInt(simpleAlignment.getSiteCount());
            int state = simpleAlignment.getState(0,pos);
            int newState = random.nextInt(4);
            while (newState == state) {
                newState = random.nextInt(4);
            }
            simpleAlignment.setState(0,pos,newState);
            defects = new Defects(simpleAlignment);
            count += defects.getDefectCount();
        }

        double p = (double)count/(double)reps;
        */

        double p = dr.stats.DiscreteStatistics.mean(pstop);
        double rate = totalSP * (1.0-p) / p;

        System.out.println("Total inferred point-mutation rate = " + rate);
        System.out.println("Total inferred point-mutations = " + (totalStops * (1-p) / p));
        System.out.println("Proportion of point-mutations that produce premature stop codons = " + p);
    }
}
