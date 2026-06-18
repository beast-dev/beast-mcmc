/*
 * SimulatedTraitParameterParser.java
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

package dr.evomodelxml.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Materializes one prior-predictive Brownian continuous-trait draw from a TreeDataLikelihood.
 *
 * The continuous trait-data parser uses taxon attributes as the source of truth,
 * so this parser also writes the simulated tip values back onto the tree taxa.
 *
 * @author Filippo Monti
 */
public final class SimulatedTraitParameterParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "simulatedTraitParameter";
    public static final String SEED = "seed";
    public static final String SET_TAXON_ATTRIBUTES = "setTaxonAttributes";

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
        if (xo.hasAttribute(SEED)) {
            MathUtils.setSeed(xo.getIntegerAttribute(SEED));
        }

        final TreeDataLikelihood treeDataLikelihood =
                (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        final Tree tree = treeDataLikelihood.getTree();
        final String traitName = xo.getAttribute(TreeTraitParserUtilities.TRAIT_NAME,
                TreeTraitParserUtilities.DEFAULT_TRAIT_NAME);

        final double[] draw = simulateBrownianTips(treeDataLikelihood);
        final int taxonCount = tree.getExternalNodeCount();
        if (taxonCount == 0 || draw.length % taxonCount != 0) {
            throw new XMLParseException(PARSER_NAME + " draw length " + draw.length +
                    " is not a multiple of external node count " + taxonCount + ".");
        }

        final boolean setTaxonAttributes = xo.getAttribute(SET_TAXON_ATTRIBUTES, true);
        if (setTaxonAttributes) {
            setTaxonAttributes(tree, traitName, draw);
        }

        final Parameter parameter = getTargetParameter(xo, draw);

        Logger.getLogger("dr.evomodel.continuous").info(
                "Simulated " + draw.length + " values for trait '" + traitName +
                        "' from " + treeDataLikelihood.getId() + " by prior-predictive Brownian recursion.");

        return parameter;
    }

    private static double[] simulateBrownianTips(final TreeDataLikelihood treeDataLikelihood) throws XMLParseException {

        if (!(treeDataLikelihood.getDataLikelihoodDelegate() instanceof ContinuousDataLikelihoodDelegate)) {
            throw new XMLParseException(PARSER_NAME + " requires a continuous traitDataLikelihood.");
        }

        final ContinuousDataLikelihoodDelegate delegate =
                (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();
        if (delegate.getTraitCount() != 1) {
            throw new XMLParseException(PARSER_NAME + " currently expects exactly one continuous trait.");
        }

        final int dim = delegate.getTraitDim();
        final Tree tree = treeDataLikelihood.getTree();
        final MultivariateDiffusionModel diffusionModel = delegate.getDiffusionModel();
        final double[][] diffusionVariance = new Matrix(diffusionModel.getPrecisionMatrix()).inverse().toComponents();
        final double[] rootMean = delegate.getRootPrior().getMean();
        final double rootSampleSize = delegate.getRootPrior().getPseudoObservations();

        if (rootSampleSize <= 0.0) {
            throw new XMLParseException(PARSER_NAME + " requires a positive root prior sample size.");
        }

        final double[] nodeValues = new double[tree.getNodeCount() * dim];
        final double[] rootDraw = drawMultivariate(rootMean, diffusionVariance, 1.0 / rootSampleSize);
        copy(rootDraw, nodeValues, tree.getRoot().getNumber() * dim, dim);

        simulateDescendants(tree, tree.getRoot(), treeDataLikelihood.getBranchRateModel(),
                diffusionVariance, nodeValues, dim);

        final int taxonCount = tree.getExternalNodeCount();
        final double[] tipDraw = new double[taxonCount * dim];
        for (int i = 0; i < taxonCount; ++i) {
            final NodeRef node = tree.getExternalNode(i);
            System.arraycopy(nodeValues, node.getNumber() * dim, tipDraw, i * dim, dim);
        }
        return tipDraw;
    }

    private static void simulateDescendants(final Tree tree,
                                            final NodeRef node,
                                            final BranchRateModel branchRateModel,
                                            final double[][] diffusionVariance,
                                            final double[] nodeValues,
                                            final int dim) throws XMLParseException {

        final int parentOffset = node.getNumber() * dim;
        final double[] parentValue = new double[dim];
        System.arraycopy(nodeValues, parentOffset, parentValue, 0, dim);

        for (int i = 0; i < tree.getChildCount(node); ++i) {
            final NodeRef child = tree.getChild(node, i);
            final double scale = getBranchVarianceScale(tree, child, branchRateModel);
            final double[] childValue = drawMultivariate(parentValue, diffusionVariance, scale);
            copy(childValue, nodeValues, child.getNumber() * dim, dim);
            simulateDescendants(tree, child, branchRateModel, diffusionVariance, nodeValues, dim);
        }
    }

    private static double getBranchVarianceScale(final Tree tree,
                                                 final NodeRef node,
                                                 final BranchRateModel branchRateModel) throws XMLParseException {
        final double branchLength = tree.getBranchLength(node);
        final double branchRate = branchRateModel == null ? 1.0 : branchRateModel.getBranchRate(tree, node);
        final double scale = branchLength * branchRate;
        if (scale < 0.0) {
            throw new XMLParseException(PARSER_NAME + " found a negative branch variance scale " + scale +
                    " on node " + node.getNumber() + ".");
        }
        return scale;
    }

    private static double[] drawMultivariate(final double[] mean,
                                             final double[][] variance,
                                             final double scale) {
        if (scale == 0.0) {
            return mean.clone();
        }
        return MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, variance, scale);
    }

    private static void copy(final double[] source,
                             final double[] destination,
                             final int destinationOffset,
                             final int length) {
        System.arraycopy(source, 0, destination, destinationOffset, length);
    }

    private static Parameter getTargetParameter(final XMLObject xo,
                                                final double[] draw) throws XMLParseException {

        Parameter parameter = null;
        if (xo.hasChildNamed(TreeTraitParserUtilities.TRAIT_PARAMETER)) {
            parameter = (Parameter) xo.getChild(TreeTraitParserUtilities.TRAIT_PARAMETER)
                    .getChild(Parameter.class);
            if (parameter.getDimension() != draw.length) {
                try {
                    parameter.setDimension(draw.length);
                } catch (UnsupportedOperationException e) {
                    throw new XMLParseException(PARSER_NAME + " target traitParameter has dimension " +
                            parameter.getDimension() + " but simulated draw has dimension " + draw.length + ".");
                }
            }
        } else {
            final String id = xo.hasId() ? xo.getId() : PARSER_NAME;
            parameter = new Parameter.Default(id, draw);
        }

        for (int i = 0; i < draw.length; ++i) {
            parameter.setParameterValue(i, draw[i]);
        }

        return parameter;
    }

    private static void setTaxonAttributes(final Tree tree,
                                           final String traitName,
                                           final double[] draw) throws XMLParseException {

        final int taxonCount = tree.getExternalNodeCount();
        final int traitDimension = draw.length / taxonCount;

        int offset = 0;
        for (int i = 0; i < taxonCount; ++i) {
            final Taxon taxon = tree.getNodeTaxon(tree.getExternalNode(i));
            if (taxon == null) {
                throw new XMLParseException(PARSER_NAME + " cannot set trait '" +
                        traitName + "' because taxon " + i + " is null.");
            }

            final StringBuilder value = new StringBuilder();
            for (int j = 0; j < traitDimension; ++j) {
                if (j > 0) {
                    value.append(' ');
                }
                value.append(Double.toString(draw[offset + j]));
            }
            taxon.setAttribute(traitName, value.toString());
            offset += traitDimension;
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return RULES;
    }

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true),
            AttributeRule.newIntegerRule(SEED, true),
            AttributeRule.newBooleanRule(SET_TAXON_ATTRIBUTES, true)
    };

    @Override
    public String getParserDescription() {
        return "Draws one continuous-trait realization from a TreeDataLikelihood and returns it as a Parameter.";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }
}
