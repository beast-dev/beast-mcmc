/*
 * SeqGen.java
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

package dr.app.seqgen;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.oldevomodel.substmodel.SubstitutionModel;
import jebl.evolution.alignments.Alignment;
import jebl.evolution.alignments.BasicAlignment;
import jebl.evolution.io.FastaExporter;
import jebl.evolution.sequences.*;
import jebl.evolution.taxa.Taxon;
import jebl.math.Random;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SeqGen {

    public SeqGen(final int length, final double substitutionRate, final FrequencyModel freqModel, final SubstitutionModel substModel, final SiteModel siteModel, final double damageRate) {
        this.length = length;
        this.substitutionRate = substitutionRate;
        this.freqModel = freqModel;
        this.substModel = substModel;
        this.siteModel = siteModel;
        this.damageRate = damageRate;
    }

    public Alignment simulate(Tree tree) {

        int[] initialSequence = new int[length];

        drawSequence(initialSequence, freqModel);

        int[] siteCategories = new int[length];

        drawSiteCategories(siteModel, siteCategories);

        double[] rates = new double[siteModel.getCategoryCount()];
        for (int i = 0; i < rates.length; i++) {
            rates[i] = siteModel.getRateForCategory(i) * substitutionRate;
        }

        for (int i = 0; i < tree.getChildCount(tree.getRoot()); i++) {
            NodeRef child = tree.getChild(tree.getRoot(), i);
            evolveSequences(initialSequence, tree, child, substModel, siteCategories, rates);
        }

        Map<State, State[]> damageMap = new HashMap<State, State[]>();
        damageMap.put(Nucleotides.A_STATE, new State[]{Nucleotides.G_STATE});
        damageMap.put(Nucleotides.C_STATE, new State[]{Nucleotides.T_STATE});
        damageMap.put(Nucleotides.G_STATE, new State[]{Nucleotides.A_STATE});
        damageMap.put(Nucleotides.T_STATE, new State[]{Nucleotides.C_STATE});

        Sequence[] sequences = new Sequence[tree.getExternalNodeCount()];

        List<NucleotideState> nucs = jebl.evolution.sequences.Nucleotides.getCanonicalStates();
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            int[] seq = (int[]) tree.getNodeTaxon(node).getAttribute("seq");
            State[] states = new State[seq.length];
            for (int j = 0; j < states.length; j++) {
                states[j] = nucs.get(seq[j]);
            }

            if (damageRate > 0) {
                damageSequence(states, damageRate, tree.getNodeHeight(node), damageMap);
            }

            sequences[i] = new BasicSequence(SequenceType.NUCLEOTIDE,
                    Taxon.getTaxon(tree.getNodeTaxon(node).getId()),
                    states);
        }

        BasicAlignment alignment = new BasicAlignment(sequences);


        return alignment;
    }

    void drawSiteCategories(SiteModel siteModel, int[] siteCategories) {
        double[] categoryProportions = siteModel.getCategoryProportions();
        double[] cumulativeProportions = new double[categoryProportions.length];
        cumulativeProportions[0] = categoryProportions[0];
        for (int i = 1; i < cumulativeProportions.length; i++) {
            cumulativeProportions[i] = cumulativeProportions[i - 1] + categoryProportions[i];
        }

        for (int i = 0; i < siteCategories.length; i++) {
            siteCategories[i] = draw(cumulativeProportions);
        }
    }

    public void drawSequence(int[] initialSequence, FrequencyModel freqModel) {
        double[] freqs = freqModel.getCumulativeFrequencies();
        for (int i = 0; i < initialSequence.length; i++) {
            initialSequence[i] = draw(freqs);
        }
    }

    void evolveSequences(int[] sequence0,
                         Tree tree, NodeRef node,
                         SubstitutionModel substModel,
                         int[] siteCategories,
                         double[] categoryRates) {
        int stateCount = substModel.getDataType().getStateCount();

        int[] sequence1 = new int[sequence0.length];

        double[][][] transitionProbabilities = new double[siteCategories.length][stateCount][stateCount];

        for (int i = 0; i < categoryRates.length; i++) {
            double branchLength = tree.getBranchLength(node) * categoryRates[i];
            double[] tmp = new double[stateCount * stateCount];
            substModel.getTransitionProbabilities(branchLength, tmp);

            int l = 0;
            for (int j = 0; j < stateCount; j++) {
                transitionProbabilities[i][j][0] = tmp[l];
                l++;
                for (int k = 1; k < stateCount; k++) {
                    transitionProbabilities[i][j][k] = transitionProbabilities[i][j][k - 1] + tmp[l];
                    l++;
                }

            }
        }

        evolveSequence(sequence0, siteCategories, transitionProbabilities, sequence1);

        if (!tree.isExternal(node)) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                evolveSequences(sequence1, tree, child, substModel, siteCategories, categoryRates);
            }
        } else {
            tree.getNodeTaxon(node).setAttribute("seq", sequence1);
        }
    }

    private void evolveSequence(int[] ancestralSequence,
                                int[] siteCategories,
                                double[][][] cumulativeTransitionProbabilities,
                                int[] descendentSequence) {

        for (int i = 0; i < ancestralSequence.length; i++) {
            descendentSequence[i] = draw(cumulativeTransitionProbabilities[siteCategories[i]][ancestralSequence[i]]);
        }
    }

    private void damageSequence(State[] sequence, double rate, double time, Map<State, State[]> damageMap) {
        double pUndamaged = Math.exp(-rate * time);

        for (int i = 0; i < sequence.length; i++) {
            double r = Random.nextDouble();
            if (r >= pUndamaged) {
                State[] states = damageMap.get(sequence[i]);
                if (states.length > 0) {
                    int index = Random.nextInt(states.length);
                    sequence[i] = states[index];
                } else {
                    sequence[i] = states[0];
                }
            }
        }
    }

    /**
     * draws a state from using a set of cumulative frequencies (last value should be 1.0)
     *
     * @param cumulativeFrequencies
     * @return
     */
    private int draw(double[] cumulativeFrequencies) {
        double r = Random.nextDouble();

        // defensive - make sure that it is actually set...
        int state = -1;

        for (int j = 0; j < cumulativeFrequencies.length; j++) {
            if (r < cumulativeFrequencies[j]) {
                state = j;
                break;
            }
        }

        assert (state != -1);

        return state;
    }

    final int length;
    final double substitutionRate;
    final FrequencyModel freqModel;
    final SubstitutionModel substModel;
    final SiteModel siteModel;
    final double damageRate;

    public static void main(String[] argv) {

        String treeFileName = argv[0];
        String outputFileStem = argv[1];

        int length = 500;

        double[] frequencies = new double[]{0.25, 0.25, 0.25, 0.25};
        double kappa = 10.0;
        double alpha = 0.5;
        double substitutionRate = argv.length < 3 ? 1.0E-3 : Double.parseDouble(argv[2]);
        int categoryCount = argv.length < 4 ? 8 : Integer.parseInt(argv[3]);
        double damageRate = argv.length < 5 ? 0 : Double.parseDouble(argv[4]); //1.56E-6;

        System.out.println("substitutionRate = " + substitutionRate  + "; categoryCount = " + categoryCount + "; damageRate = " + damageRate);

        FrequencyModel freqModel = new FrequencyModel(dr.evolution.datatype.Nucleotides.INSTANCE, frequencies);

        HKY hkyModel = new HKY(kappa, freqModel);
        SiteModel siteModel = null;

        if (categoryCount > 1) {
            siteModel = new GammaSiteModel(hkyModel, alpha, categoryCount);
        } else {
            // no rate heterogeneity
            siteModel = new GammaSiteModel(hkyModel);
        }

        List<Tree> trees = new ArrayList<Tree>();

        FileReader reader = null;
        try {
            reader = new FileReader(treeFileName);
//            TreeImporter importer = new NexusImporter(reader);
            TreeImporter importer = new NewickImporter(reader);

            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                trees.add(tree);
                System.out.println("tree height = " + tree.getNodeHeight(tree.getRoot()) + "; leave nodes = " + tree.getExternalNodeCount());
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return;
        } catch (Importer.ImportException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        SeqGen seqGen = new SeqGen(length,
                substitutionRate, freqModel, hkyModel, siteModel,
                damageRate);
        int i = 1;
        for (Tree tree : trees) {
            Alignment alignment = seqGen.simulate(tree);

            FileWriter writer = null;
            try {
//                writer = new FileWriter(outputFileStem + (i < 10 ? "00" : (i < 100 ? "0" : "")) + i + ".nex");
//                NexusExporter exporter = new NexusExporter(writer);
//
//                exporter.exportAlignment(alignment);
//
//                writer.close();
                String outputFileName = outputFileStem + "-" + substitutionRate + ".fasta";

                writer = new FileWriter(outputFileName);
                BufferedWriter bf = new BufferedWriter(writer);
                FastaExporter exporter = new FastaExporter(bf);

                exporter.exportSequences(alignment.getSequenceList());

                bf.close();
                System.out.println("Write " + i + "th sequence file : " + outputFileName);

                i++;
            } catch (IOException e) {
                e.printStackTrace();
                return;
            }
        }

    }
}
