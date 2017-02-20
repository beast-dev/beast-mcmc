/*
 * IstvanOperator.java
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

package dr.oldevomodel.indel;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import org.w3c.dom.Document;
import org.w3c.dom.Element;


/**
 * Implements branch exchange operations.
 * There is a NARROW and WIDE variety.
 * The narrow exchange is very similar to a rooted-tree
 * nearest-neighbour interchange but with the restriction
 * that node height must remain consistent.
 */
public class IstvanOperator extends SimpleMCMCOperator {

    private double tuning = 0.1;
    private double exponent = 1.5;
    private double gapPenalty = -10;
    private TKF91Likelihood likelihood;
    private IstvansProposal proposal = new IstvansProposal();
    private Alignment alignment;

    int[][] iAlignment;
    double[][][] iProbs;
    double[] iBaseFreqs;
    int[] iParent;
    double[] iTau;

    public IstvanOperator(double iP, double exponent, double gapPenalty, double weight, TKF91Likelihood likelihood) {
        this.tuning = iP;
        this.exponent = exponent;
        this.gapPenalty = gapPenalty;
        this.likelihood = likelihood;
        setWeight(weight);
    }

    public double doOperation() {

        Tree tree = likelihood.getTreeModel();
        alignment = likelihood.getAlignment();

        //System.out.println("Incoming alignment");
        //System.out.println(alignment);
        //System.out.println();

        SubstitutionModel substModel = likelihood.getSiteModel().getSubstitutionModel();

        // initialize the iParent and iTau arrays based on the given tree.
        initTree(tree, likelihood.getSiteModel().getMu());

        int[] treeIndex = new int[tree.getTaxonCount()];
        for (int i = 0; i < treeIndex.length; i++) {
            treeIndex[i] = tree.getTaxonIndex(alignment.getTaxonId(i));
        }

        // initialize the iAlignment array from the given alignment.
        initAlignment(alignment, treeIndex);

        // initialize the iProbs array from the substitution model -- must be called after populating tree!
        initSubstitutionModel(substModel);

        DataType dataType = substModel.getDataType();
        proposal.setGapSymbol(dataType.getGapState());

        int[][] returnedAlignment = new int[iAlignment.length][];

        //System.out.println("Initialization done, starting proposal proper...");


        double logq = proposal.propose(iAlignment, iProbs, iBaseFreqs, iParent, iTau, returnedAlignment, tuning, exponent, gapPenalty);

        //System.out.println("Proposal finished, logq=" + logq);

        //create new alignment object
        SimpleAlignment newAlignment = new SimpleAlignment();
        for (int i = 0; i < alignment.getTaxonCount(); i++) {

            StringBuffer seqBuffer = new StringBuffer();
            for (int j = 0; j < returnedAlignment[i].length; j++) {
                seqBuffer.append(dataType.getChar(returnedAlignment[treeIndex[i]][j]));
            }

            // add sequences in order of tree
            String seqString = seqBuffer.toString();
            Sequence sequence = new Sequence(alignment.getTaxon(i), seqString);
            newAlignment.addSequence(sequence);

            String oldunaligned = alignment.getUnalignedSequenceString(i);
            String unaligned = newAlignment.getUnalignedSequenceString(i);
            if (!unaligned.equals(oldunaligned)) {
                System.err.println("Sequence changed from:");
                System.err.println("old:'" + oldunaligned + "'");
                System.err.println("new:'" + unaligned + "'");
                throw new RuntimeException();
            }
        }
        //System.out.println("Outgoing alignment");
        //System.out.println(newAlignment);
        //System.out.println();


        likelihood.setAlignment(newAlignment);

        return logq;
    }

    // MUST RESET ALIGNMENT IF REJECTED!!
    public void reject() {
        super.reject();
        likelihood.setAlignment(alignment);
    }

    private void initTree(Tree tree, double mutationRate) {
        iParent = new int[tree.getNodeCount()];
        iTau = new double[tree.getNodeCount() - 1];
        populate(tree, tree.getRoot(), new int[]{tree.getExternalNodeCount()}, mutationRate);
        iParent[tree.getNodeCount() - 1] = -1;

    }

    /**
     * initialize the iProbs array from the substitution model -- must be called after populating tree!
     */
    private void initSubstitutionModel(SubstitutionModel model) {

        DataType dataType = model.getDataType();
        int stateCount = dataType.getStateCount();

        iProbs = new double[iTau.length][stateCount][stateCount];

        double[] transProb = new double[stateCount * stateCount];
        int count;
        for (int i = 0; i < iTau.length; i++) {
            model.getTransitionProbabilities(iTau[i], transProb);
            count = 0;
            for (int j = 0; j < stateCount; j++) {
                for (int k = 0; k < stateCount; k++) {
                    iProbs[i][j][k] = transProb[count];
                    count += 1;
                }
            }
        }

        // initialize equlibrium distribution
        iBaseFreqs = new double[stateCount];
        for (int k = 0; k < stateCount; k++) {
            iBaseFreqs[k] = model.getFrequencyModel().getFrequency(k);
        }
    }

    /**
     * Initializes the iAlignment array from the given alignment.
     */
    private void initAlignment(Alignment alignment, int[] treeIndex) {

        int numSeqs = alignment.getSequenceCount();
        int numSites = alignment.getSiteCount();
        DataType dataType = alignment.getDataType();
        int numStates = dataType.getStateCount();

        iAlignment = new int[numSeqs][numSites];

        // populate alignment in order of tree
        for (int i = 0; i < numSeqs; i++) {
            for (int j = 0; j < numSites; j++) {
                iAlignment[treeIndex[i]][j] = alignment.getState(i, j);
            }
        }
    }

    /**
     * Populates the iParent and iTau arrays.
     *
     * @return the node number of the given node.
     */
    private int populate(Tree tree, NodeRef node, int[] current, double mutationRate) {

        int nodeNumber = node.getNumber();

        // if its an external node just return the number
        if (tree.isExternal(node)) {
            iTau[nodeNumber] =
                    (tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node)) * mutationRate;
            return nodeNumber;
        }

        // if internal node, first let your children be assigned numbers
        int[] childNumbers = new int[tree.getChildCount(node)];
        for (int i = 0; i < tree.getChildCount(node); i++) {
            childNumbers[i] = populate(tree, tree.getChild(node, i), current, mutationRate);
        }

        // now, pick the next available number
        nodeNumber = current[0];
        // if you are not the root, then record the branch length above you.
        if (!tree.isRoot(node)) {
            //iTau[nodeNumber] = tree.getBranchLength(node) * mutationRate;
            iTau[nodeNumber] =
                    (tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node)) * mutationRate;
        }
        // increment the next available number
        current[0] += 1;

        // now that you have your number, populate the iParent entries of your children.
        for (int i = 0; i < tree.getChildCount(node); i++) {
            iParent[childNumbers[i]] = nodeNumber;
        }

        // finally return your number so your parent can do the same.
        return nodeNumber;
    }

    public String getOperatorName() {
        return "IstvansOperator";
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.4;
    }

    public String getPerformanceSuggestion() {
        if (MCMCOperator.Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (MCMCOperator.Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public Element createOperatorElement(Document doc) {
        throw new RuntimeException();
    }

}
