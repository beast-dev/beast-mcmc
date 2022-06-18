/*
 * RateEpochBranchRateModel.java
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
import dr.evomodelxml.branchratemodel.RateEpochBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * Implements a model where time is broken into 'epochs' each with a different but
 * constant rate. Parameters can be used to sample transition times but it is up
 * to the user to keep them bounded and in strict order...
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class RateEpochBranchRateModel extends AbstractBranchRateModel {

    protected final Parameter[] timeParameters;
    protected final Parameter[] rateParameters;

    /**
     * The constructor. For an N-epoch model, there should be N rate paramters and N-1 transition times.
     *
     * @param timeParameters an array of transition time parameters
     * @param rateParameters an array of rate parameters
     */
    public RateEpochBranchRateModel(Parameter[] timeParameters,
                                    Parameter[] rateParameters) {

        super(RateEpochBranchRateModelParser.RATE_EPOCH_BRANCH_RATES);

        this.timeParameters = timeParameters;
        for (Parameter parameter : timeParameters) {
            addVariable(parameter);
        }

        this.rateParameters = rateParameters;
        for (Parameter parameter : rateParameters) {
            addVariable(parameter);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // nothing to do
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        NodeRef parent = tree.getParent(node);

        if (parent != null) {
            double height0 = tree.getNodeHeight(node);
            double height1 = tree.getNodeHeight(parent);
            int i = 0;

            double rate = 0.0;
            double lastHeight = height0;

            // First find the epoch which contains the node height
            while (i < timeParameters.length && height0 > timeParameters[i].getParameterValue(0)) {
                i++;
            }

            // Now walk up the branch until we reach the last epoch or the height of the parent
            while (i < timeParameters.length && height1 > timeParameters[i].getParameterValue(0)) {
                // add the rate for that epoch multiplied by the time spent at that rate
                rate += rateParameters[i].getParameterValue(0) * (timeParameters[i].getParameterValue(0) - lastHeight);
                lastHeight = timeParameters[i].getParameterValue(0);
                i++;
            }

            // Add that last rate segment
            rate += rateParameters[i].getParameterValue(0) * (height1 - lastHeight);

            // normalize the rate for the branch length
            return normalizeRate(rate / (height1 - height0));
        }
        throw new IllegalArgumentException("root node doesn't have a rate!");
    }
    
    @Override
    public Mapping getBranchRateModelMapping(final Tree tree, final NodeRef node) {
        
        NodeRef parent = tree.getParent(node);
        
        List<Double> weightList = new ArrayList<Double>();
        List<Double> rateList = new ArrayList<Double>();

        if (parent != null) {
            double height0 = tree.getNodeHeight(node);
            double height1 = tree.getNodeHeight(parent);
            int i = 0;

            double rate = 0.0;
            double lastHeight = height0;

            // First find the epoch which contains the node height
            while (i < timeParameters.length && height0 >= timeParameters[i].getParameterValue(0)) {
                i++;
            }

            // Now walk up the branch until we reach the last epoch or the height of the parent
            while (i < timeParameters.length && height1 > timeParameters[i].getParameterValue(0)) {
                // insert each epoch that this branch overlaps with to the list,
                // ordered from parent (rootward) to child (tipward),
                // as the transition probability matrix of a given branch is the product
                // of matrices multiplying from parent to child
                weightList.add( 0, timeParameters[i].getParameterValue(0) - lastHeight );
                rateList.add( 0, rateParameters[i].getParameterValue(0) );
                lastHeight = timeParameters[i].getParameterValue(0);
                i++;
            }

            // Add that last rate segment
            weightList.add( 0, height1 - lastHeight );
            rateList.add( 0, rateParameters[i].getParameterValue(0) );

        } else {
            throw new IllegalArgumentException("root node doesn't have a rate!");
        }
        
        if (weightList.size() == 0) {
            throw new RuntimeException("RateEpochBranchRateModel failed to give a valid mapping");
        }
        
        final double[] rates = new double[rateList.size()];
        final double[] weights = new double[weightList.size()];
        for (int i = 0; i < weightList.size(); i++) {
            rates[i] = rateList.get(i);
            weights[i] = weightList.get(i);
        }
        
        return new Mapping() {
            public double[] getRates() {
                return rates;
            }

            public double[] getWeights() {
                return weights;
            }
        };
    }

    protected double normalizeRate(double rate) {
        return rate;
    }

}
