/*
 * RateEvolutionLikelihood.java
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

package dr.oldevomodel.clock;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.logging.Logger;

/**
 * Abstract superclass of likelihoods of rate evolution through time.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Michael Defoin Platel
 */
public abstract class RateEvolutionLikelihood extends AbstractBranchRateModel {

    public static final String RATES = "rates";
    public static final String EPISODIC = "episodic";
    public static final String LOGSPACE = "logspace";

    public static final String ROOTRATE = "rootRate";

    public RateEvolutionLikelihood(String name, TreeModel treeModel, Parameter ratesParameter, Parameter rootRateParameter, boolean isEpisodic) {

        super(name);

        this.treeModel = treeModel;
        addModel(treeModel);

        this.ratesParameter = new TreeParameterModel(treeModel, ratesParameter, false);
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(Double.MAX_VALUE, 0, ratesParameter.getDimension());
        ratesParameter.addBounds(bound);

        addModel(this.ratesParameter);

        this.rootRateParameter = rootRateParameter;
        rootRateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0, 1));
        addVariable(rootRateParameter);

        if (rootRateParameter.getDimension() != 1) {
            throw new IllegalArgumentException("The root rate parameter must be of dimension 1");
        }

        this.isEpisodic = isEpisodic;

        Logger.getLogger("dr.evomodel").info("AutoCorrelated Relaxed Clock: " + name + (isEpisodic ? " (episodic)." : "."));

    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    public final void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
        if (model == ratesParameter) {
            fireModelChanged(this, index);
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
    }

    protected void storeState() {
    }

    protected void restoreState() {
        likelihoodKnown = false;
    }

    protected void acceptState() {
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Get the model.
     *
     * @return the model.
     */
    public Model getModel() {
        return this;
    }

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
     *
     * @return the likelihood.
     */
    protected boolean getLikelihoodKnown() {
        return likelihoodKnown;
    }

    /**
     * Get the log likelihood of the rate changes in this tree.
     *
     * @return the log likelihood.
     */
    private double calculateLogLikelihood() {
        NodeRef root = treeModel.getRoot();
        NodeRef node1 = treeModel.getChild(root, 0);
        NodeRef node2 = treeModel.getChild(root, 1);

        return calculateLogLikelihood(root, node1) + calculateLogLikelihood(root, node2);
    }

    /**
     * Recursively calculate the log likelihood of the rate changes in the given tree.
     *
     * @return the partial log likelihood of the rate changes below the given node plus the
     *         branch directly above.
     */
    private double calculateLogLikelihood(NodeRef parent, NodeRef node) {

        double logL, length;
        length = treeModel.getBranchLength(node);

        logL = branchRateChangeLogLikelihood(getBranchRate(treeModel, parent), getBranchRate(treeModel, node),
                length);

        //System.out.print(parent.getNumber() + " " + getBranchRate(treeModel, parent)+ " " + node.getNumber() + " " + getBranchRate(treeModel, node) + " " + treeModel.getBranchLength(node) + " " + logL + ", ");

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            logL += calculateLogLikelihood(node, treeModel.getChild(node, i));
        }
        return logL;
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }

    abstract double branchRateSample(double parentRate, double time);

    public void sampleRate(NodeRef node) {

        final NodeRef parent = treeModel.getParent(node);
        final double length = treeModel.getBranchLength(node);
        final double rate = branchRateSample(getBranchRate(treeModel, parent), length);

        treeModel.setNodeRate(node, rate);

    }

    public double getBranchRate(Tree tree, NodeRef node) {

        if (tree.isRoot(node)) return rootRateParameter.getParameterValue(0);
        return ratesParameter.getNodeValue(tree, node);
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

    private double logLikelihood;
    private boolean likelihoodKnown = false;

    private final TreeModel treeModel;
    private final TreeParameterModel ratesParameter;
    protected final Parameter rootRateParameter;
    private final boolean isEpisodic;

}
