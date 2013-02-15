/*
 * AbstractARGLikelihood.java
 *
 * (c) 2002-2005 BEAST Development Core Team
 *
 * This package may be distributed under the
 * Lesser Gnu Public Licence (LGPL)
 */

package dr.evomodel.arg.likelihood;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.arg.ARGModel;
import dr.evomodel.treelikelihood.LikelihoodCore;
import dr.inference.model.*;

/**
 * AbstractTreeLikelihood - a base class for likelihood calculators of sites on a tree.
 *
 * @author Andrew Rambaut
 * @version $Id: AbstractARGLikelihood.java,v 1.1 2006/10/10 22:57:55 msuchard Exp $
 */

public abstract class AbstractARGLikelihood extends AbstractModelLikelihood implements ParallelLikelihood {

    public AbstractARGLikelihood(String name, PatternList patternList,
                                 ARGModel treeModel) {

        super(name);

        this.patternList = patternList;
        this.dataType = patternList.getDataType();
        patternCount = patternList.getPatternCount();
        stateCount = dataType.getStateCount();

        patternWeights = patternList.getPatternWeights();

        this.treeModel = treeModel;
        addModel(treeModel);

        nodeCount = treeModel.getNodeCount();

        updateNode = new boolean[nodeCount];
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }

    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    protected final void setStates(LikelihoodCore likelihoodCore, PatternList patternList,
                                   int sequenceIndex, int nodeIndex) {
        int i;

        int[] states = new int[patternCount];

        for (i = 0; i < patternCount; i++) {

            states[i] = patternList.getPatternState(sequenceIndex, i);
            //System.err.print(states[i]+" ");
        }

        likelihoodCore.setNodeStates(nodeIndex, states);
    }

    /**
     * Sets the partials from a sequence in an alignment.
     */
    protected final void setPartials(LikelihoodCore likelihoodCore, PatternList patternList,
                                     int categoryCount,
                                     int sequenceIndex, int nodeIndex) {
        int i, j;

        double[] partials = new double[patternCount * stateCount];

        boolean[] stateSet;

        int v = 0;
        for (i = 0; i < patternCount; i++) {

            int state = patternList.getPatternState(sequenceIndex, i);
            stateSet = dataType.getStateSet(state);

            for (j = 0; j < stateCount; j++) {
                if (stateSet[j]) {
                    partials[v] = 1.0;
                } else {
                    partials[v] = 0.0;
                }
                v++;
            }
        }

        likelihoodCore.setNodePartials(nodeIndex, partials);
    }

    /**
     * Set update flag for a node and its children
     */
    protected void updateNode(NodeRef node) {

        updateNode[node.getNumber()] = true;
        likelihoodKnown = false;
    }

    /**
     * Set update flag for a node and its children
     */
    protected void updateNodeAndChildren(NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNode[child.getNumber()] = true;
        }
        likelihoodKnown = false;
    }   

    /**
     * Set update flag for a node and its children
     */
    protected void updateNodeAndDescendents(NodeRef node) {
        updateNode[node.getNumber()] = true;

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            updateNodeAndDescendents(child);
        }

        likelihoodKnown = false;
    }

    /**
     * Set update flag for all nodes
     */
    protected void updateAllNodes() {
        for (int i = 0; i < nodeCount; i++) {
            updateNode[i] = true;
        }
        likelihoodKnown = false;
    }

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

    public final double[] getPatternWeights() {
        return patternWeights;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    /**
     * Stores the additional state other than model components
     */
    protected void storeState() {

        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restore the additional stored state
     */
    protected void restoreState() {

        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    public void makeDirty() {
        likelihoodKnown = false;
        updateAllNodes();
        updateAllPatterns();
    }

    public void setLikelihood(double likelihood) {
        this.logLikelihood = likelihood;
        likelihoodKnown = true;
    }


    public boolean getLikelihoodKnown() {
        return likelihoodKnown;
    }

    protected abstract double calculateLogLikelihood();

    public String toString() {

        return Double.toString(getLogLikelihood());

    }


    // **************************************************************
    // INSTANCE VARIABLES
    // **************************************************************

    /**
     * the tree
     */
    protected ARGModel treeModel = null;

    /**
     * the partition *
     */
    protected int partition;

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
     * the number of nodes in the tree
     */
    protected int nodeCount;

    /**
     * Flags to specify which patterns are to be updated
     */
    protected boolean[] updatePattern = null;

    /**
     * Flags to specify which nodes are to be updated
     */
    protected boolean[] updateNode;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

}