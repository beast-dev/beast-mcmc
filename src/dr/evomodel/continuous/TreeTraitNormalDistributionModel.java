/*
 * TreeTraitNormalDistributionModel.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.RandomGenerator;

/**
 * A class that acts as a model for multivariate normally distributed data.
 *
 * @author Marc Suchard
 * @author Mandev Gill
 */

public class TreeTraitNormalDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel, RandomGenerator {

    public TreeTraitNormalDistributionModel(FullyConjugateMultivariateTraitLikelihood traitModel, boolean conditionOnRoot) {

        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.traitModel = traitModel;
        this.conditionOnRoot = conditionOnRoot;
        dim = traitModel.getTreeModel().getExternalNodeCount() * traitModel.getDimTrait();
        addModel(traitModel);
        distributionKnown = false;
    }

    // *****************************************************************
    // Interface MultivariateDistribution
    // *****************************************************************

    public double logPdf(double[] x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    public double[][] getScaleMatrix() {
        checkDistribution();
        return distribution.getScaleMatrix();
    }

    public double[] getMean() {
        checkDistribution();
        return distribution.getMean();
    }

    public String getType() {
        return "TreeTraitMVN";
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void handleModelChangedEvent(Model model, Object object, int index) {
        distributionKnown = false;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        distributionKnown = false;
    }

    protected void storeState() {
        storedDistribution = distribution;
        storedDistributionKnown = distributionKnown;
    }

    protected void restoreState() {
        distributionKnown = storedDistributionKnown;
        distribution = storedDistribution;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    // **************************************************************
    // Private instance variables and functions
    // **************************************************************

    private void checkDistribution() {
        if (!distributionKnown) {
            mean = null;
            precision = null;
            distribution = createNewDistribution();
            distributionKnown = true;
        }
    }

    private MultivariateNormalDistribution createNewDistribution() {
        return new MultivariateNormalDistribution(computeMean(), computePrecision());
    }

    private double[] computeMean() {
        return MultivariateTraitUtils.computeTreeTraitMean(traitModel, conditionOnRoot);
    }

    private double[][] computePrecision() {
        return MultivariateTraitUtils.computeTreeTraitPrecision(traitModel, conditionOnRoot);
    }

    private final FullyConjugateMultivariateTraitLikelihood traitModel;

    private double[] mean;
    private double[][] precision;

    private MultivariateNormalDistribution distribution;
    private MultivariateNormalDistribution storedDistribution;

    private boolean distributionKnown;
    private boolean storedDistributionKnown;

    // RandomGenerator interface
    public Object nextRandom() {
        checkDistribution();
        return distribution.nextMultivariateNormal();
    }

    public double logPdf(Object x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    private final boolean conditionOnRoot;
    private double[][] precisionMatrix = null;
    private final int dim;

}
