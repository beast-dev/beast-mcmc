/*
 * DecayingRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id: DecayingRateModel.java,v 1.3 2006/01/10 16:48:27 rambaut Exp $
 */
public class DecayingRateModel extends AbstractBranchRateModel {

    public static final String DECAYING_RATE_MODEL = "decayingRateModel";
    public static final String MUTATION_RATE = "mutationRate";
    public static final String HALF_LIFE = "halfLife";
    public static final String SUBSTITUTION_RATE = "substitutionRate";
    public static final String PROPORTION = "proportion";
    public static final String AVERAGE = "average";

    private final Parameter mutationRateParameter;
    private final Parameter substitutionRateParameter;
    private final Parameter proportionParameter;
    private final Parameter halfLifeParameter;

    private final boolean useAveraging;

    private final double[] rates;
    private boolean ratesCalculated = false;

    public DecayingRateModel(TreeModel treeModel, Parameter mutationRateParameter, Parameter substitutionRateParameter,
                             Parameter proportionParameter, Parameter halfLifeParameter, boolean useAveraging) {

        super(DECAYING_RATE_MODEL);

        rates = new double[treeModel.getNodeCount()];

        addModel(treeModel);

        this.mutationRateParameter = mutationRateParameter;
        this.substitutionRateParameter = substitutionRateParameter;
        this.proportionParameter = proportionParameter;
        this.halfLifeParameter = halfLifeParameter;

        addVariable(mutationRateParameter);
        if (proportionParameter != null) {
            addVariable(proportionParameter);
        } else {
            addVariable(substitutionRateParameter);
        }
        addVariable(halfLifeParameter);

        this.useAveraging = useAveraging;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // TreeModel has changed...
        ratesCalculated = false;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        // Parameters have changed
        ratesCalculated = false;
        fireModelChanged();
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        ratesCalculated = false;
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        if (!ratesCalculated) {
            double mu = mutationRateParameter.getParameterValue(0);
            double k;
            if (proportionParameter != null) {
                k = mu * proportionParameter.getParameterValue(0);
            } else {
                k = substitutionRateParameter.getParameterValue(0);
            }

            double lambda = Math.log(2) / halfLifeParameter.getParameterValue(0);

            calculateNodeRates(tree, tree.getRoot(), mu, k, lambda);

            ratesCalculated = true;
        }

        return rates[node.getNumber()];
    }

    /**
     * Traverse the tree calculating partial likelihoods.
     *
     * @return whether the partials for this node were recalculated.
     */
    private final double calculateNodeRates(Tree tree, NodeRef node, double mu, double k, double lambda) {

        NodeRef parent = tree.getParent(node);

        double time0 = 0.0;

        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            double t1 = calculateNodeRates(tree, child1, mu, k, lambda);

            NodeRef child2 = tree.getChild(node, 1);
            double t2 = calculateNodeRates(tree, child2, mu, k, lambda);

            if (useAveraging) {
                time0 = (t1 + t2) / 2.0;
            } else {
                // pick larger of the two
                if (t1 > t2) {
                    time0 = t1;
                } else {
                    time0 = t2;
                }
            }
        }

        // don't bother if you are at the root because rate at root is ignored
        if (parent == null) return 0;

        double branchTime = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
        double time1 = time0 + branchTime;

        double branchRate = rateIntegral(time1, mu, k, lambda);

        if (time0 > 0.0) {
            branchRate -= rateIntegral(time0, mu, k, lambda);
        }

        rates[node.getNumber()] = branchRate / branchTime;

        return time1;
    }

    private static double rateIntegral(double time, double mu, double k, double lambda) {
        return ((k * time) + (((mu - k) / lambda) * (1.0 - Math.exp(-lambda * time))));
    }

    private static double rate(double time, double mu, double k, double lambda) {
        return k + ((mu - k) * Math.exp(-lambda * time));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DECAYING_RATE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean useAveraging = xo.getBooleanAttribute(AVERAGE);

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            Parameter mutationRateParameter = (Parameter) xo.getElementFirstChild(MUTATION_RATE);

            Parameter proportionParameter = null;
            Parameter substitutionRateParameter = null;

            if (xo.hasChildNamed(PROPORTION)) {
                proportionParameter = (Parameter) xo.getElementFirstChild(PROPORTION);
            } else {
                substitutionRateParameter = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATE);
            }
            Parameter halfLifeParameter = (Parameter) xo.getElementFirstChild(HALF_LIFE);

            Logger.getLogger("dr.evomodel").info("Using decaying-rate clock model.");

            return new DecayingRateModel(treeModel, mutationRateParameter, substitutionRateParameter, proportionParameter, halfLifeParameter, useAveraging);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element provides a clock model in which the rate of " +
                            "molecular evolution decays to a substitution rate in the past.";
        }

        public Class getReturnType() {
            return DecayingRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(AVERAGE, false),
                new ElementRule(TreeModel.class, "The tree model"),
                new ElementRule(MUTATION_RATE, Parameter.class, "The mutation rate parameter", false),
                new XORRule(
                        new ElementRule(PROPORTION, Parameter.class, "The proportion of neutral mutations", false),
                        new ElementRule(SUBSTITUTION_RATE, Parameter.class, "The long-term substitution", false)
                ),
                new ElementRule(HALF_LIFE, Parameter.class, "The half-life of a deleterious mutation", false)
        };
    };

    public static void main(String[] argv) {
        double mu = 1E-4;
        double k = 1E-5;
        double lambda = Math.log(2) / 1.0;

        double time1 = 0.000001;
        for (int i = 0; i < 100; i++) {
            System.out.println(time1 + "\t" + rate(time1, mu, k, lambda) + "\t" + (rateIntegral(time1, mu, k, lambda)) / time1);
            time1 += 1.0;
        }
    }
}