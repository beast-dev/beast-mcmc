/*
 * NormalizedSequenceLikelihoodTest.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package test.dr.evomodel.treelikelihood;

import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class NormalizedSequenceLikelihoodTest extends SequenceLikelihoodTest {

    public NormalizedSequenceLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {

        format.setMaximumFractionDigits(5);

        int numTaxa = TWO_TAXON_SEQUENCE[0].length;

        createAlignmentWithAllUniquePatterns(TWO_TAXON_SEQUENCE, Nucleotides.INSTANCE);

        treeModel = createSillyTreeModel(numTaxa);
    }


    protected TreeModel createSillyTreeModel(int numTaxa) {

        SimpleNode[] nodes = new SimpleNode[2 * numTaxa - 1];
        for (int n = 0; n < 2 * numTaxa - 1; n++) {
            nodes[n] = new SimpleNode();
        }

        nodes[0].setTaxon(taxa[0]);
        nodes[1].setTaxon(taxa[1]);

        nodes[2].setHeight(1);
        nodes[2].addChild(nodes[0]);
        nodes[2].addChild(nodes[1]);

        SimpleNode root = nodes[2];

        if (numTaxa == 3) {
            nodes[3].setTaxon(taxa[2]);
            nodes[4].setHeight(2);
            nodes[4].addChild(nodes[2]);
            nodes[4].addChild(nodes[3]);
            root = nodes[4];
        }

        Tree tree = new SimpleTree(root);
        tree.setUnits(Units.Type.YEARS);

        return new DefaultTreeModel(tree); //treeModel
    }

    public void testAllPossibleAlignments() {
        SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
        tryAllPossibleAlignments(3, patterns);
    }

    private void tryAllPossibleAlignments(int length, SitePatterns patterns) {
        System.out.println("Trying all possible " + length + " site alignments");
        double[] patternLogLikeihoods = computeSitePatternLikelihoods(patterns);

        System.out.println("Site logLikelihoods: " + new Vector(patternLogLikeihoods));

        double total = 0.0;

        List<Double> allAlignmentLogProbabilities = new ArrayList<Double>();
        recursivelyComputeAlignmentLikelihood(allAlignmentLogProbabilities, patternLogLikeihoods, length, 0, 0);

        System.out.println("Total possible alignments: " + allAlignmentLogProbabilities.size());

        for (Double x : allAlignmentLogProbabilities) {
            total += Math.exp(x);
        }
        System.out.println("Total probability = " + total);
        assertEquals("uncorrected", 1.0, total, tolerance);
    }

    private void recursivelyComputeAlignmentLikelihood(List<Double> finalLogProbabilities, double[] patternLogLikelihoods,
                                                       int alignmentLength, int level, double logProbability) {

        if (level < alignmentLength) {
            for (int i = 0; i < patternLogLikelihoods.length; i++) {
                double thisLogProb = logProbability + patternLogLikelihoods[i];
                recursivelyComputeAlignmentLikelihood(finalLogProbabilities, patternLogLikelihoods, alignmentLength,
                        level + 1, thisLogProb);
            }
        } else {
            finalLogProbabilities.add(logProbability);
        }
    }

    protected static final String[][] TWO_TAXON_SEQUENCE = {{"human", "chimp", "gorrila"},
            {"","", ""}};

}
