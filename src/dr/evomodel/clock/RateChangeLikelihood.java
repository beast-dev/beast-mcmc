/*
 * RateChangeLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.clock;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.inference.model.*;

import java.util.logging.Logger;

/**
 * Abstract superclass of likelihoods of rate change through time.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: RateChangeLikelihood.java,v 1.13 2005/07/11 14:06:25 rambaut Exp $
 */
public abstract class RateChangeLikelihood extends AbstractModel implements BranchRateModel, Likelihood {

    public static final int ROOT_RATE_MEAN_OF_CHILDREN = 0;
    public static final int ROOT_RATE_MEAN_OF_ALL = 1;
    public static final int ROOT_RATE_EQUAL_TO_CHILD = 2;
    public static final int ROOT_RATE_IGNORE_ROOT = 3;
    public static final int ROOT_RATE_NONE = 4;

    public static final String RATES = "rates";
    public static final String EPISODIC = "episodic";
    public static final String ROOT_MODEL = "rootModel";
    public static final String MEAN_OF_CHILDREN = "meanOfChildren";
    public static final String MEAN_OF_ALL = "meanOfAll";
    public static final String EQUAL_TO_CHILD = "equalToChild";
    public static final String IGNORE_ROOT = "ignoreRoot";
    public static final String NONE = "none";

    // the index of the root node.
    private int rootNodeNumber;

    private final int rateCount;
    private final double[] rates;

    private boolean ratesKnown = false;

    public RateChangeLikelihood(String name, TreeModel treeModel, Parameter ratesParameter, int rootModel, boolean isEpisodic) {

		super(name);

		this.treeModel = treeModel;
        addModel(treeModel);

        this.ratesParameter = ratesParameter;
        addParameter(ratesParameter);

        rateCount = treeModel.getNodeCount() - 1;
        rates = new double[rateCount];

        if (ratesParameter.getDimension() != rateCount ) {
            throw new IllegalArgumentException("The rate category parameter must be of dimension nodeCount-1");
        }

        ratesKnown = false;

        rootNodeNumber = treeModel.getRoot().getNumber();

        if (rootModel < ROOT_RATE_MEAN_OF_CHILDREN || rootModel > ROOT_RATE_NONE) {
            throw new IllegalArgumentException();
        }

        this.rootModel = rootModel;
        this.isEpisodic = isEpisodic;

        String rootModelName = MEAN_OF_CHILDREN;
        if (rootModel == ROOT_RATE_MEAN_OF_CHILDREN) rootModelName = MEAN_OF_CHILDREN;
        if (rootModel == ROOT_RATE_MEAN_OF_ALL) rootModelName = MEAN_OF_ALL;
        if (rootModel == ROOT_RATE_EQUAL_TO_CHILD) rootModelName = EQUAL_TO_CHILD;
        if (rootModel == ROOT_RATE_IGNORE_ROOT) rootModelName = IGNORE_ROOT;
        if (rootModel == ROOT_RATE_NONE) rootModelName = NONE;
	    Logger.getLogger("dr.evomodel").info("Relaxed clock: " + name + " model, root model = " + rootModelName + (isEpisodic? " (episodic).": "."));

    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    public final void handleModelChangedEvent(Model model, Object object, int index) {

        ratesKnown = false;
        fireModelChanged();
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        ratesKnown = false;
        fireModelChanged();
    }

    protected void storeState() { }

    protected void restoreState() {
        ratesKnown = false;
    }

    protected void acceptState() {}

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Get the model.
     * @return the model.
     */
    public Model getModel() { return this; }

    public final double getLogLikelihood() {
        if (!getLikelihoodKnown()) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Called to decide if the likelihood must be calculated. Can be overridden
     * (for example, to always return false).
     */
    protected boolean getLikelihoodKnown() {
        return likelihoodKnown;
    }

    /**
     * Get the log likelihood of the rate changes in this tree.
     * @return the log likelihood.
     */
    private final double calculateLogLikelihood() {
        NodeRef root = treeModel.getRoot();
        NodeRef node1 = treeModel.getChild(root, 0);
        NodeRef node2 = treeModel.getChild(root, 1);

        if (rootModel == ROOT_RATE_IGNORE_ROOT || rootModel == ROOT_RATE_NONE) {

            double  logL = 0.0;

            if (rootModel == ROOT_RATE_IGNORE_ROOT) {
                logL += branchRateChangeLogLikelihood(getBranchRate(treeModel, node1), getBranchRate(treeModel, node2),
                                                        (treeModel.getBranchLength(node1)) +
                                                        (treeModel.getBranchLength(node2)));
            }

            for (int i = 0; i < treeModel.getChildCount(node1); i++) {
                logL += calculateLogLikelihood(node1, treeModel.getChild(node1, i));
            }
            for (int i = 0; i < treeModel.getChildCount(node2); i++) {
                logL += calculateLogLikelihood(node2, treeModel.getChild(node2, i));
            }
            return logL;

        } else {

            return calculateLogLikelihood(root, node1) + calculateLogLikelihood(root, node2);
        }


    }

    /**
     * Recursively calculate the log likelihood of the rate changes in the given tree.
     * @return the partial log likelihood of the rate changes below the given node plus the
     * branch directly above.
     */
    private final double calculateLogLikelihood(NodeRef parent, NodeRef node) {

        double logL = branchRateChangeLogLikelihood(getBranchRate(treeModel, parent), getBranchRate(treeModel, node),
                                                    treeModel.getBranchLength(node));

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            logL += calculateLogLikelihood(node, treeModel.getChild(node, i));
        }
        return logL;
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[] {
            new LikelihoodColumn(getId())
        };
    }

    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) { super(label); }
        public double getDoubleValue() { return getLogLikelihood(); }
    }

    public double getBranchRate(Tree tree, NodeRef node) {

        if (!ratesKnown) {
            setupRates();
            ratesKnown = true;
        }

        return rates[node.getNumber()];
    }

    public boolean isEpisodic() {
        return isEpisodic;
    }

    /**
     * @return the log likelihood of the rate change from the parent to the given node.
     */
    abstract double branchRateChangeLogLikelihood(double parentRate, double childRate, double time);

    // **************************************************************
    // Private members
    // **************************************************************

    /**
     * Sets the appropriate root rate for this model of rate change.
     */
    private final void setupRates() {

        rootNodeNumber = treeModel.getRoot().getNumber();

        int j = 0;
        for (int i = 0; i < rateCount; i++) {
            if (i != rootNodeNumber) {
                rates[i] = ratesParameter.getParameterValue(j);
                j++;
            }
        }

        double rootRate = -1;

        NodeRef node1 = treeModel.getChild(treeModel.getRoot(), 0);
        NodeRef node2 = treeModel.getChild(treeModel.getRoot(), 1);

        double nodeRate1 = rates[node1.getNumber()];
        double nodeRate2 = rates[node2.getNumber()];

        if (!isEpisodic) {
            // interpolate between the two rates to determine the rate at the root
            // i.e root_rate = (r1*b2 + r2*b1) / (b1+b2)
            double b1 = treeModel.getBranchLength(node1);
            double b2 = treeModel.getBranchLength(node2);
            rootRate = (nodeRate1 * b2 + nodeRate2 * b1) / (b1 + b2);
        } else {
            if (rootModel == ROOT_RATE_MEAN_OF_CHILDREN) {

                rootRate = (nodeRate1 + nodeRate2) / 2.0;

            } else if (rootModel == ROOT_RATE_MEAN_OF_ALL) {

                double sum = 0;
                for (int i = 0; i < rateCount; i++) {
                    if (i != rootNodeNumber) {
                        sum += rates[i];
                    }
                }

                rootRate = sum / (rateCount - 1);

            } else if (rootModel == ROOT_RATE_EQUAL_TO_CHILD) {

                rootRate = nodeRate1;

            } else if (rootModel == ROOT_RATE_IGNORE_ROOT || rootModel == ROOT_RATE_NONE) {

                // The root rate is not used so set it to something silly...
                rootRate = -1;
            }
        }

        rates[rootNodeNumber] = rootRate;


    }

     public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }


    private double logLikelihood;
    private boolean likelihoodKnown = false;

    private final TreeModel treeModel;
    private final Parameter ratesParameter;
    private final int rootModel;
    private final boolean isEpisodic;

}