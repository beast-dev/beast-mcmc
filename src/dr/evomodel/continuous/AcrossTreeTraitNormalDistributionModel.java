/*
 * TreeTraitNormalDistributionModel.java
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

package dr.evomodel.continuous;

import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inferencexml.distribution.MultivariateNormalDistributionModelParser;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.RandomGenerator;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * A class that acts as a model for multivariate normally distributed data.
 *
 * @author Marc Suchard
 * @author Philippe Lemey
 */

public class AcrossTreeTraitNormalDistributionModel extends AbstractModel implements ParametricMultivariateDistributionModel, RandomGenerator {

    private final ContinuousDataLikelihoodDelegate delegate1;
    private final ContinuousDataLikelihoodDelegate delegate2;
    private final Parameter rho;
    private final int dim;


    public AcrossTreeTraitNormalDistributionModel(ContinuousDataLikelihoodDelegate delegate1,
                                                  ContinuousDataLikelihoodDelegate delegate2,
                                                  Parameter rho) {

        super(MultivariateNormalDistributionModelParser.NORMAL_DISTRIBUTION_MODEL);

        this.delegate1 = delegate1;
        this.delegate2 = delegate2;
        this.rho = rho;
        this.dim = delegate1.getTreeTraitPrecision().length;

        if (delegate1.getTraitDim() != delegate2.getTraitDim()) {
            throw new IllegalArgumentException("Unequal trait dimensions");
        }

        if (delegate1.getTraitCount() != 1 || delegate2.getTraitCount() != 1) {
            throw new IllegalArgumentException("Only implemented for single traits");
        }

        addModel(delegate1);
        addModel(delegate2);
        addVariable(rho);

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
        if (model == delegate1 || model == delegate2) {
            distributionKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rho) {
            distributionKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown correlation");
        }
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
            distribution = createNewDistribution();
            distributionKnown = true;
        }
    }

    private MultivariateNormalDistribution createNewDistribution() {
        return new MultivariateNormalDistribution(computeMean(), computePrecision());
    }

    private double[] computeMean() {
        return new double[2 * dim]; // TODO Make work for non-zero means

//        if (traitModel.strengthOfSelection != null) {
//            return MultivariateTraitUtils.computeTreeTraitMeanOU(traitModel, rootValue, conditionOnRoot);
//        } else {
//            return MultivariateTraitUtils.computeTreeTraitMean(traitModel, rootValue, conditionOnRoot);
//        }
    }

    private double[][] variance;

    private double[][] computePrecision() {

        double[][][] var = new double[2][][];
        var[0] = delegate1.getTreeTraitVariance();
        var[1] = delegate2.getTreeTraitVariance();

        if (variance == null) {
            variance = new double[2 * dim][2 * dim];
        }

        for (int tree = 0; tree < 2; ++tree) {

            double[][] v = var[tree];

            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    variance[tree * dim + i][tree * dim + j] = v[i][j];
                }
            }

            double r = rho.getParameterValue(0);
            for (int tree2 = tree + 1; tree2 < 2; ++tree2) {
                for (int i = 0; i < dim; ++i) {

                    double sigma1 = Math.sqrt(variance[tree * dim + i][tree * dim + i]);
                    double sigma2 = Math.sqrt(variance[tree2 * dim + i][tree2 * dim + i]);

                    variance[tree * dim + i][tree2 * dim + i] = r * sigma1 * sigma2;
                    variance[tree2 * dim + i][tree * dim + i] = r * sigma1 + sigma2;
                }
            }
        }

        return new SymmetricMatrix(variance).inverse().toComponents();
    }

    private MultivariateNormalDistribution distribution;
    private MultivariateNormalDistribution storedDistribution;

    private boolean distributionKnown;
    private boolean storedDistributionKnown;

    // RandomGenerator interface
    public double[] nextRandom() {
        checkDistribution();
        return distribution.nextMultivariateNormal();
    }

    @Override
    public int getDimension() {
        return dim;
    }

    public double logPdf(Object x) {
        checkDistribution();
        return distribution.logPdf(x);
    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

    public static XMLObjectParser ACROSS_TREE_TRAIT_MODEL = new AbstractXMLObjectParser() {

        private static final String ACROSS_TREE_TRAIT_NORMAL = "acrossTreeTraitNormalDistribution";

        public String getParserName() {
            return ACROSS_TREE_TRAIT_NORMAL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            List<TreeDataLikelihood> likelihoods = xo.getAllChildren(TreeDataLikelihood.class);
            List<ContinuousDataLikelihoodDelegate> delegates = new ArrayList<>();
            for (TreeDataLikelihood likelihood : likelihoods) {
                DataLikelihoodDelegate delegate = likelihood.getDataLikelihoodDelegate();
                if (delegate instanceof ContinuousDataLikelihoodDelegate) {
                    delegates.add((ContinuousDataLikelihoodDelegate) delegate);
                }
            }

            Parameter rho = (Parameter) xo.getChild(Parameter.class);

            return new AcrossTreeTraitNormalDistributionModel(delegates.get(0), delegates.get(1), rho);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeDataLikelihood.class, 2, 2),
                new ElementRule(Parameter.class),
        };

        public String getParserDescription() {
            return "Parses TreeTraitNormalDistributionModel";
        }

        public Class getReturnType() {
            return TreeTraitNormalDistributionModel.class;
        }
    };
}
