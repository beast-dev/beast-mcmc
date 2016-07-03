/*
 * BeagleOperationReport.java
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

package dr.evomodel.treelikelihood;

import beagle.Beagle;
import dr.evomodelxml.treelikelihood.BeagleOperationParser;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tipstatesmodel.TipStatesModel;
import dr.math.matrixAlgebra.Vector;

import java.io.PrintWriter;

/**
 * BeagleTreeLikelihoodModel - implements a Likelihood Function for sequences on a tree.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Marc Suchard
 * @version $Id$
 */

public class BeagleOperationReport extends AbstractSinglePartitionTreeLikelihood {

    /**
     * the patternList
     */
    protected PatternList patternList = null;
    protected DataType dataType = null;
    /**
     * the pattern weights
     */
    protected double[] patternWeights;
    /**
     * the number of patterns
     */
    protected int patternCount;
    /**
     * the number of states in the data
     */
    protected int stateCount;
    /**
     * Flags to specify which patterns are to be updated
     */
    protected boolean[] updatePattern = null;

    public BeagleOperationReport(TreeModel treeModel, PatternList patternList, BranchRateModel branchRateModel, GammaSiteRateModel siteRateModel, Alignment alignment, PrintWriter branch, PrintWriter operation) {
        super(BeagleOperationParser.OPERATION_REPORT, patternList, treeModel);

        boolean useAmbiguities = false;
        this.branchRateModel = branchRateModel;
        this.branchWriter = branch;
        this.operationWriter = operation;
        this.alignment = alignment;
        this.substitutionModel = siteRateModel.getSubstitutionModel();

        try {

            this.tipCount = treeModel.getExternalNodeCount();

            internalNodeCount = nodeCount - tipCount;

            int compactPartialsCount = tipCount;

            // one partials buffer for each tip and two for each internal node (for store restore)
            partialBufferHelper = new BufferIndexHelper(nodeCount, tipCount);

            // two eigen buffers for each decomposition for store and restore.
            eigenBufferHelper = new BufferIndexHelper(eigenCount, 0);

            // two matrices for each node less the root
            matrixBufferHelper = new BufferIndexHelper(nodeCount, 0);

            for (int i = 0; i < tipCount; i++) {
                // Find the id of tip i in the patternList
                String id = treeModel.getTaxonId(i);
                int index = patternList.getTaxonIndex(id);

                if (index == -1) {
                    throw new TaxonList.MissingTaxonException("Taxon, " + id + ", in tree, " + treeModel.getId() +
                            ", is not found in patternList, " + patternList.getId());
                } else {

                    if (useAmbiguities) {
                        setPartials(beagle, patternList, index, i);
                    } else {
                        setStates(beagle, patternList, id, index, i);
                    }

                }
            }


        } catch (TaxonList.MissingTaxonException mte) {
            throw new RuntimeException(mte.toString());
        }
        hasInitialized = true;
    }


    public String toString() {

        calculateLogLikelihood();
        return super.toString();
    }


    public TreeModel getTreeModel() {
        return treeModel;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternList   patternList
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    protected final void setPartials(Beagle beagle,
                                     PatternList patternList,
                                     int sequenceIndex,
                                     int nodeIndex) {
        double[] partials = new double[patternCount * stateCount * categoryCount];

        boolean[] stateSet;

        int v = 0;
        for (int i = 0; i < patternCount; i++) {

            int state = patternList.getPatternState(sequenceIndex, i);
            stateSet = dataType.getStateSet(state);

            for (int j = 0; j < stateCount; j++) {
                if (stateSet[j]) {
                    partials[v] = 1.0;
                } else {
                    partials[v] = 0.0;
                }
                v++;
            }
        }

        // if there is more than one category then replicate the partials for each
        int n = patternCount * stateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        System.err.println("TODO Print partials");
//        beagle.setPartials(nodeIndex, partials);
    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    protected final void setPartials(Beagle beagle,
                                     TipStatesModel tipStatesModel,
                                     int nodeIndex) {
        double[] partials = new double[patternCount * stateCount * categoryCount];

        tipStatesModel.getTipPartials(nodeIndex, partials);

        // if there is more than one category then replicate the partials for each
        int n = patternCount * stateCount;
        int k = n;
        for (int i = 1; i < categoryCount; i++) {
            System.arraycopy(partials, 0, partials, k, n);
            k += n;
        }

        System.err.println("TODO Print partials");
//        beagle.setPartials(nodeIndex, partials);
    }

    public int getPatternCount() {
        return patternCount;
    }

    /**
     * Sets the partials from a sequence in an alignment.
     *
     * @param beagle        beagle
     * @param patternList   patternList
     * @param id
     * @param sequenceIndex sequenceIndex
     * @param nodeIndex     nodeIndex
     */
    protected final void setStates(Beagle beagle,
                                   PatternList patternList,
                                   String id, int sequenceIndex,
                                   int nodeIndex) {
        int i;

        StringBuilder sb = new StringBuilder();
        sb.append("/* ").append(id).append(" */\n\t\tmSeqs[").append(nodeIndex).append("] = \"");
        sb.append(alignment.getAlignedSequenceString(sequenceIndex)).append("\";\n");

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {
            states[i] = patternList.getPatternState(sequenceIndex, i);
        }

        if (alignmentString == null) {
            alignmentString = new StringBuilder();
        }
        alignmentString.append(sb);
    }


    protected double calculateLogLikelihood() {

        if (matrixUpdateIndices == null) {
            matrixUpdateIndices = new int[eigenCount][nodeCount];
            branchLengths = new double[eigenCount][nodeCount];
            branchUpdateCount = new int[eigenCount];
//            scaleBufferIndices = new int[internalNodeCount];
//            storedScaleBufferIndices = new int[internalNodeCount];
        }

        if (operations == null) {
            operations = new int[numRestrictedPartials + 1][internalNodeCount * Beagle.OPERATION_TUPLE_SIZE];
            operationCount = new int[numRestrictedPartials + 1];
        }

        recomputeScaleFactors = false;

        for (int i = 0; i < eigenCount; i++) {
            branchUpdateCount[i] = 0;
        }
        operationListCount = 0;

        if (hasRestrictedPartials) {
            for (int i = 0; i <= numRestrictedPartials; i++) {
                operationCount[i] = 0;
            }
        } else {
            operationCount[0] = 0;
        }

        System.out.println(alignmentString.toString());

        final NodeRef root = treeModel.getRoot();
        traverse(treeModel, root, null, false); // Do not flip buffers

        // Print out eigendecompositions


        for (int i = 0; i < eigenCount; i++) {
            if (branchUpdateCount[i] > 0) {

                if (DEBUG_BEAGLE_OPERATIONS) {
                    StringBuilder sb = new StringBuilder();

                    sb.append("eval = ").append(new Vector(substitutionModel.getEigenDecomposition().getEigenValues())).append("\n");
                    sb.append("evec = ").append(new Vector(substitutionModel.getEigenDecomposition().getEigenVectors())).append("\n");
                    sb.append("ivec = ").append(new Vector(substitutionModel.getEigenDecomposition().getInverseEigenVectors())).append("\n");

                    sb.append("Branch count: ").append(branchUpdateCount[i]);
                    sb.append("\nNode indices:\n");
                    if (SINGLE_LINE) {
                        sb.append("int n[] = {");
                    }
                    for (int k = 0; k < branchUpdateCount[i]; ++k) {
                        if (SINGLE_LINE) {
                            sb.append(" ").append(matrixUpdateIndices[i][k]);
                            if (k < (branchUpdateCount[i] - 1)) {
                                sb.append(",");
                            }
                        } else {
                            sb.append(matrixUpdateIndices[i][k]).append("\n");
                        }
                    }
                    if (SINGLE_LINE) {
                        sb.append(" };\n");
                    }
                    sb.append("\nBranch lengths:\n");
                    if (SINGLE_LINE) {
                        sb.append("double b[] = {");
                    }
                    for (int k = 0; k < branchUpdateCount[i]; ++k) {
                        if (SINGLE_LINE) {
                            sb.append(" ").append(branchLengths[i][k]);
                            if (k < (branchUpdateCount[i] - 1)) {
                                sb.append(",");
                            }
                        } else {
                            sb.append(branchLengths[i][k]).append("\n");
                        }
                    }
                    if (SINGLE_LINE) {
                        sb.append(" };\n");
                    }
                    System.out.println(sb.toString());
                }
            }
        }


        if (DEBUG_BEAGLE_OPERATIONS) {
            StringBuilder sb = new StringBuilder();
            sb.append("Operation count: ").append(operationCount[0]);
            sb.append("\nOperations:\n");
            if (SINGLE_LINE) {
                sb.append("BeagleOperation o[] = {");
            }
            for (int k = 0; k < operationCount[0] * Beagle.OPERATION_TUPLE_SIZE; ++k) {
                if (SINGLE_LINE) {
                    sb.append(" ").append(operations[0][k]);
                    if (k < (operationCount[0] * Beagle.OPERATION_TUPLE_SIZE - 1)) {
                        sb.append(",");
                    }
                } else {
                    sb.append(operations[0][k]).append("\n");
                }
            }
            if (SINGLE_LINE) {
                sb.append(" };\n");
            }
            sb.append("Use scale factors: ").append(useScaleFactors).append("\n");
            System.out.println(sb.toString());

        }


        int rootIndex = partialBufferHelper.getOffsetIndex(root.getNumber());

        System.out.println("Root node: " + rootIndex);

        return 0.0;
    }


    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @param tree           tree
     * @param node           node
     * @param operatorNumber operatorNumber
     * @param flip           flip
     * @return boolean
     */
    private boolean traverse(Tree tree, NodeRef node, int[] operatorNumber, boolean flip) {

        boolean update = false;

        int nodeNum = node.getNumber();

        NodeRef parent = tree.getParent(node);

        if (operatorNumber != null) {
            operatorNumber[0] = -1;
        }

        // First update the transition probability matrix(ices) for this branch
        if (parent != null && updateNode[nodeNum]) {

            final double branchRate = branchRateModel.getBranchRate(tree, node);

            // Get the operational time of the branch
            final double branchTime = branchRate * (tree.getNodeHeight(parent) - tree.getNodeHeight(node));

            if (branchTime < 0.0) {
                throw new RuntimeException("Negative branch length: " + branchTime);
            }

            if (flip) {
                // first flip the matrixBufferHelper
                matrixBufferHelper.flipOffset(nodeNum);
            }

            // then set which matrix to update
            final int eigenIndex = 0; //branchSubstitutionModel.getBranchIndex(tree, node);
            final int updateCount = branchUpdateCount[eigenIndex];
            matrixUpdateIndices[eigenIndex][updateCount] = matrixBufferHelper.getOffsetIndex(nodeNum);

            branchLengths[eigenIndex][updateCount] = branchTime;
            branchUpdateCount[eigenIndex]++;

            update = true;
        }

        // If the node is internal, update the partial likelihoods.
        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            final int[] op1 = {-1};
            final boolean update1 = traverse(tree, child1, op1, flip);

            NodeRef child2 = tree.getChild(node, 1);
            final int[] op2 = {-1};
            final boolean update2 = traverse(tree, child2, op2, flip);

            // If either child node was updated then update this node too
            if (update1 || update2) {

                int x = operationCount[operationListCount] * Beagle.OPERATION_TUPLE_SIZE;

                if (flip) {
                    // first flip the partialBufferHelper
                    partialBufferHelper.flipOffset(nodeNum);
                }

                final int[] operations = this.operations[operationListCount];

                operations[x] = partialBufferHelper.getOffsetIndex(nodeNum);

                if (useScaleFactors) {
                    // get the index of this scaling buffer
                    int n = nodeNum - tipCount;

                    if (recomputeScaleFactors) {
                        // flip the indicator: can take either n or (internalNodeCount + 1) - n
                        scaleBufferHelper.flipOffset(n);

                        // store the index
                        scaleBufferIndices[n] = scaleBufferHelper.getOffsetIndex(n);

                        operations[x + 1] = scaleBufferIndices[n]; // Write new scaleFactor
                        operations[x + 2] = Beagle.NONE;

                    } else {
                        operations[x + 1] = Beagle.NONE;
                        operations[x + 2] = scaleBufferIndices[n]; // Read existing scaleFactor
                    }

                } else {

                    if (useAutoScaling) {
                        scaleBufferIndices[nodeNum - tipCount] = partialBufferHelper.getOffsetIndex(nodeNum);
                    }
                    operations[x + 1] = Beagle.NONE; // Not using scaleFactors
                    operations[x + 2] = Beagle.NONE;
                }

                operations[x + 3] = partialBufferHelper.getOffsetIndex(child1.getNumber()); // source node 1
                operations[x + 4] = matrixBufferHelper.getOffsetIndex(child1.getNumber()); // source matrix 1
                operations[x + 5] = partialBufferHelper.getOffsetIndex(child2.getNumber()); // source node 2
                operations[x + 6] = matrixBufferHelper.getOffsetIndex(child2.getNumber()); // source matrix 2

                operationCount[operationListCount]++;

                update = true;
            }
        }

        return update;

    }

    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    private int eigenCount = 1;
    private int[][] matrixUpdateIndices;
    private double[][] branchLengths;
    private int[] branchUpdateCount;
    private int[] scaleBufferIndices;
    private int[] storedScaleBufferIndices;

    private int[][] operations;
    private int operationListCount;
    private int[] operationCount;
    private static final boolean hasRestrictedPartials = false;

    private final int numRestrictedPartials = 0;

    protected BufferIndexHelper partialBufferHelper;
    private final BufferIndexHelper eigenBufferHelper;
    protected BufferIndexHelper matrixBufferHelper;
    protected BufferIndexHelper scaleBufferHelper;

    protected final int tipCount;
    protected final int internalNodeCount;

    private PartialsRescalingScheme rescalingScheme;
    protected boolean useScaleFactors = false;
    private boolean useAutoScaling = false;
    private boolean recomputeScaleFactors = false;
    private boolean everUnderflowed = false;
    private int rescalingCount = 0;
    private int rescalingCountInner = 0;

    protected final BranchRateModel branchRateModel;
    protected double[] patternLogLikelihoods = null;

    /**
     * the number of rate categories
     */
    protected int categoryCount;

    /**
     * an array used to transfer tip partials
     */
    protected double[] tipPartials;

    /**
     * an array used to transfer tip states
     */
    protected int[] tipStates;

    /**
     * the BEAGLE library instance
     */
    protected Beagle beagle;

    /**
     * Flag to specify that the substitution model has changed
     */
    protected boolean updateSubstitutionModel;

    /**
     * Flag to specify that the site model has changed
     */
    protected boolean updateSiteModel;

    private static final boolean DEBUG_BEAGLE_OPERATIONS = true;
    private static final boolean SINGLE_LINE = true;

    private StringBuilder alignmentString;


    private final PrintWriter branchWriter;
    private final PrintWriter operationWriter;
    private final SubstitutionModel substitutionModel;

    private final Alignment alignment;

    /**
     * Set update flag for a pattern
     */
    protected void updatePattern(int i) {
        if (updatePattern != null) {
            updatePattern[i] = true;
        }
        likelihoodKnown = false;
    }

    /**
     * Set update flag for all patterns
     */
    protected void updateAllPatterns() {
        if (updatePattern != null) {
            for (int i = 0; i < patternCount; i++) {
                updatePattern[i] = true;
            }
        }
        likelihoodKnown = false;
    }

    protected class BufferIndexHelper {
        /**
         * @param maxIndexValue the number of possible input values for the index
         * @param minIndexValue the minimum index value to have the mirrored buffers
         */
        BufferIndexHelper(int maxIndexValue, int minIndexValue) {
            this.maxIndexValue = maxIndexValue;
            this.minIndexValue = minIndexValue;

            offsetCount = maxIndexValue - minIndexValue;
            indexOffsets = new int[offsetCount];
            storedIndexOffsets = new int[offsetCount];
        }

        public int getBufferCount() {
            return 2 * offsetCount + minIndexValue;
        }

        void flipOffset(int i) {
            if (i >= minIndexValue) {
                indexOffsets[i - minIndexValue] = offsetCount - indexOffsets[i - minIndexValue];
            } // else do nothing
        }

        int getOffsetIndex(int i) {
            if (i < minIndexValue) {
                return i;
            }
            return indexOffsets[i - minIndexValue] + i;
        }

        void getIndices(int[] outIndices) {
            for (int i = 0; i < maxIndexValue; i++) {
                outIndices[i] = getOffsetIndex(i);
            }
        }

        void storeState() {
            System.arraycopy(indexOffsets, 0, storedIndexOffsets, 0, indexOffsets.length);

        }

        void restoreState() {
            int[] tmp = storedIndexOffsets;
            storedIndexOffsets = indexOffsets;
            indexOffsets = tmp;
        }

        private final int maxIndexValue;
        private final int minIndexValue;
        private final int offsetCount;

        private int[] indexOffsets;
        private int[] storedIndexOffsets;

    }
}