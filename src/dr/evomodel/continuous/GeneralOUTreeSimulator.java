/*
 * GeneralOUTreeSimulator.java
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

import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.ConjugateRootTraitPrior;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.GaussianBranchTransitionKernel;
import dr.inferencexml.timeseries.OUSelectionChartParserHelper;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Exact simulator for a constant-optimum multivariate OU process on a rooted tree.
 *
 * <p>The branch transition is sampled from the same shared canonical OU kernel
 * used elsewhere in the codebase:
 * <pre>
 *   X_child | X_parent ~ N(F(dt) X_parent + b(dt), V(dt))
 * </pre>
 * with
 * <pre>
 *   F(dt) = exp(-A dt),  b(dt) = (I - F(dt)) mu.
 * </pre>
 *
 * <p>The root state is sampled separately from a user-provided Gaussian root
 * distribution. A convenience constructor maps the standard conjugate root prior
 * used in tree likelihoods to the corresponding Gaussian root covariance
 * {@code Sigma / pseudoObservations}.</p>
 */
public final class GeneralOUTreeSimulator {

    public static final String GENERAL_OU_TREE_SIMULATOR = "generalOuTreeSimulator";
    public static final String TRAIT_NAME = "traitName";
    public static final String CLONE = "clone";
    public static final String ROOT_MEAN = "rootMean";
    public static final String ROOT_COVARIANCE = "rootCovariance";
    public static final String OPTIMUM_MEAN = "optimalTraits";
    public static final String STRENGTH_OF_SELECTION_MATRIX = "strengthOfSelectionMatrix";

    private final Tree tree;
    private final MultivariateElasticModel elasticModel;
    private final MultivariateDiffusionModel diffusionModel;
    private final Parameter optimumMean;
    private final Parameter rootMean;
    private final MatrixParameterInterface explicitRootCovariance;
    private final ConjugateRootTraitPrior conjugateRootPrior;
    private final int dimension;

    private final MatrixParameter diffusionCovarianceSnapshot;
    private final MatrixParameter unusedInitialCovariance;
    private final OUProcessModel processModel;
    private final GaussianBranchTransitionKernel kernel;
    private final DenseMatrix64F precisionMatrix;
    private final DenseMatrix64F covarianceMatrix;

    public GeneralOUTreeSimulator(final Tree tree,
                                  final MultivariateElasticModel elasticModel,
                                  final MultivariateDiffusionModel diffusionModel,
                                  final Parameter optimumMean,
                                  final Parameter rootMean,
                                  final MatrixParameterInterface rootCovariance) {
        if (tree == null) {
            throw new IllegalArgumentException("tree must not be null");
        }
        if (elasticModel == null) {
            throw new IllegalArgumentException("elasticModel must not be null");
        }
        if (diffusionModel == null) {
            throw new IllegalArgumentException("diffusionModel must not be null");
        }
        if (optimumMean == null) {
            throw new IllegalArgumentException("optimumMean must not be null");
        }
        if (rootMean == null) {
            throw new IllegalArgumentException("rootMean must not be null");
        }
        if (rootCovariance == null) {
            throw new IllegalArgumentException("rootCovariance must not be null");
        }
        this.tree = tree;
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        this.optimumMean = optimumMean;
        this.rootMean = rootMean;
        this.explicitRootCovariance = rootCovariance;
        this.conjugateRootPrior = null;
        this.dimension = elasticModel.getStrengthOfSelectionMatrixParameter().getRowDimension();

        validateMeanDimension(optimumMean, "optimumMean");
        validateVectorParameter(rootMean, dimension, "rootMean");
        validateSquareMatrix(rootCovariance, dimension, "rootCovariance");

        this.diffusionCovarianceSnapshot = new MatrixParameter("generalOuTreeSimulator.diffusion", dimension, dimension);
        this.unusedInitialCovariance = new MatrixParameter("generalOuTreeSimulator.initial", dimension, dimension);
        setZero(unusedInitialCovariance);
        this.processModel = new OUProcessModel(
                "generalOuTreeSimulator.process",
                dimension,
                elasticModel.getStrengthOfSelectionMatrixParameter(),
                diffusionCovarianceSnapshot,
                optimumMean,
                unusedInitialCovariance);
        this.kernel = processModel.getRepresentation(GaussianBranchTransitionKernel.class);
        this.precisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.covarianceMatrix = new DenseMatrix64F(dimension, dimension);
    }

    public GeneralOUTreeSimulator(final Tree tree,
                                  final MultivariateElasticModel elasticModel,
                                  final MultivariateDiffusionModel diffusionModel,
                                  final Parameter optimumMean,
                                  final ConjugateRootTraitPrior rootPrior) {
        if (rootPrior == null) {
            throw new IllegalArgumentException("rootPrior must not be null");
        }
        if (tree == null) {
            throw new IllegalArgumentException("tree must not be null");
        }
        if (elasticModel == null) {
            throw new IllegalArgumentException("elasticModel must not be null");
        }
        if (diffusionModel == null) {
            throw new IllegalArgumentException("diffusionModel must not be null");
        }
        if (optimumMean == null) {
            throw new IllegalArgumentException("optimumMean must not be null");
        }
        this.tree = tree;
        this.elasticModel = elasticModel;
        this.diffusionModel = diffusionModel;
        this.optimumMean = optimumMean;
        this.rootMean = rootPrior.getMeanParameter();
        this.explicitRootCovariance = null;
        this.conjugateRootPrior = rootPrior;
        this.dimension = elasticModel.getStrengthOfSelectionMatrixParameter().getRowDimension();

        validateMeanDimension(optimumMean, "optimumMean");
        validateVectorParameter(rootMean, dimension, "rootPrior.mean");

        this.diffusionCovarianceSnapshot = new MatrixParameter("generalOuTreeSimulator.diffusion", dimension, dimension);
        this.unusedInitialCovariance = new MatrixParameter("generalOuTreeSimulator.initial", dimension, dimension);
        setZero(unusedInitialCovariance);
        this.processModel = new OUProcessModel(
                "generalOuTreeSimulator.process",
                dimension,
                elasticModel.getStrengthOfSelectionMatrixParameter(),
                diffusionCovarianceSnapshot,
                optimumMean,
                unusedInitialCovariance);
        this.kernel = processModel.getRepresentation(GaussianBranchTransitionKernel.class);
        this.precisionMatrix = new DenseMatrix64F(dimension, dimension);
        this.covarianceMatrix = new DenseMatrix64F(dimension, dimension);
    }

    public int getDimension() {
        return dimension;
    }

    /**
     * Simulate traits for every node in the tree, indexed by node number.
     */
    public double[][] simulateNodeTraits() {
        refreshDiffusionSnapshot();
        final double[][] nodeTraits = new double[tree.getNodeCount()][];
        simulateFromNode(tree.getRoot(), drawRootState(), nodeTraits);
        return nodeTraits;
    }

    /**
     * Simulate a new tree with the sampled traits written as node attributes.
     * External taxa also receive a space-separated trait string.
     */
    public Tree simulateAnnotatedTree(final String traitName, final boolean cloneTree) {
        if (traitName == null || traitName.length() == 0) {
            throw new IllegalArgumentException("traitName must not be empty");
        }

        final MutableTree mutableTree;
        if (cloneTree) {
            final FlexibleTree cloned = new FlexibleTree(tree);
            cloned.resolveTree();
            mutableTree = cloned;
        } else if (tree instanceof MutableTree) {
            mutableTree = (MutableTree) tree;
        } else {
            throw new IllegalArgumentException("tree must be mutable when cloneTree=false");
        }

        final double[][] nodeTraits = simulateNodeTraits();
        for (int i = 0; i < mutableTree.getNodeCount(); ++i) {
            final NodeRef node = mutableTree.getNode(i);
            final double[] trait = nodeTraits[node.getNumber()].clone();
            mutableTree.setNodeAttribute(node, traitName, trait);
            if (mutableTree.isExternal(node) && mutableTree.getNodeTaxon(node) != null) {
                mutableTree.getNodeTaxon(node).setAttribute(traitName, formatTrait(trait));
            }
        }
        return mutableTree;
    }

    private void simulateFromNode(final NodeRef node,
                                  final double[] state,
                                  final double[][] nodeTraits) {
        nodeTraits[node.getNumber()] = state.clone();

        final int childCount = tree.getChildCount(node);
        if (childCount == 0) {
            return;
        }

        final double[][] transitionMatrix = new double[dimension][dimension];
        final double[] transitionOffset = new double[dimension];
        final double[][] transitionCovariance = new double[dimension][dimension];
        final double[] childMean = new double[dimension];

        for (int i = 0; i < childCount; ++i) {
            final NodeRef child = tree.getChild(node, i);
            final double dt = getBranchLength(child);
            kernel.fillTransitionMatrix(dt, transitionMatrix);
            kernel.fillTransitionOffset(dt, transitionOffset);
            kernel.fillTransitionCovariance(dt, transitionCovariance);

            multiply(transitionMatrix, state, childMean);
            addInPlace(childMean, transitionOffset);

            final double[] childState = drawGaussian(childMean, transitionCovariance);
            simulateFromNode(child, childState, nodeTraits);
        }
    }

    private double[] drawRootState() {
        return drawGaussian(expandVector(rootMean), getRootCovariance());
    }

    private double[][] getRootCovariance() {
        if (explicitRootCovariance != null) {
            return copyMatrix(explicitRootCovariance.getParameterAsMatrix());
        }

        final double[][] diffusionCovariance = copyMatrix(diffusionCovarianceSnapshot.getParameterAsMatrix());
        final double pseudoObservations = conjugateRootPrior.getPseudoObservations();
        if (Double.isInfinite(pseudoObservations)) {
            return new double[dimension][dimension];
        }
        if (pseudoObservations <= 0.0) {
            throw new IllegalArgumentException("root prior pseudoObservations must be positive");
        }
        final double scale = 1.0 / pseudoObservations;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                diffusionCovariance[i][j] *= scale;
            }
        }
        return diffusionCovariance;
    }

    private double[] drawGaussian(final double[] mean, final double[][] covariance) {
        if (isEffectivelyZero(covariance)) {
            return mean.clone();
        }
        return MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, covariance);
    }

    private void refreshDiffusionSnapshot() {
        final double[][] precision = diffusionModel.getPrecisionMatrix();
        final double[] precisionData = precisionMatrix.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                precisionData[i * dimension + j] = precision[i][j];
            }
        }
        covarianceMatrix.set(precisionMatrix);
        CommonOps.invert(covarianceMatrix);
        copyIntoMatrixParameter(covarianceMatrix, diffusionCovarianceSnapshot);
        processModel.fireModelChanged();
    }

    private double getBranchLength(final NodeRef node) {
        final double branchLength = tree.getBranchLength(node);
        if (!Double.isNaN(branchLength)) {
            return branchLength;
        }
        final NodeRef parent = tree.getParent(node);
        return tree.getNodeHeight(parent) - tree.getNodeHeight(node);
    }

    private static double[] expandVector(final Parameter parameter) {
        final int dim = parameter.getDimension();
        if (dim == 1) {
            throw new IllegalArgumentException("Vector parameter must already have expanded dimension here");
        }
        return parameter.getParameterValues();
    }

    private static void validateMeanDimension(final Parameter mean, final String label) {
        if (mean.getDimension() < 1) {
            throw new IllegalArgumentException(label + " must have positive dimension");
        }
    }

    private static void validateVectorParameter(final Parameter parameter,
                                                final int expectedDimension,
                                                final String label) {
        if (parameter.getDimension() != expectedDimension) {
            throw new IllegalArgumentException(label + " dimension must equal " + expectedDimension);
        }
    }

    private static void validateSquareMatrix(final MatrixParameterInterface matrix,
                                             final int expectedDimension,
                                             final String label) {
        if (matrix.getRowDimension() != expectedDimension || matrix.getColumnDimension() != expectedDimension) {
            throw new IllegalArgumentException(label + " must be " + expectedDimension + "x" + expectedDimension);
        }
    }

    private static void multiply(final double[][] matrix, final double[] vector, final double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < matrix.length; ++i) {
            double value = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                value += matrix[i][j] * vector[j];
            }
            out[i] = value;
        }
    }

    private static void addInPlace(final double[] target, final double[] increment) {
        for (int i = 0; i < target.length; ++i) {
            target[i] += increment[i];
        }
    }

    private static boolean isEffectivelyZero(final double[][] matrix) {
        for (double[] row : matrix) {
            for (double value : row) {
                if (Math.abs(value) > 1.0e-15) {
                    return false;
                }
            }
        }
        return true;
    }

    private static double[][] copyMatrix(final double[][] matrix) {
        final double[][] copy = new double[matrix.length][];
        for (int i = 0; i < matrix.length; ++i) {
            copy[i] = matrix[i].clone();
        }
        return copy;
    }

    private static void copyIntoMatrixParameter(final DenseMatrix64F source, final MatrixParameter destination) {
        final int dimension = source.numRows;
        final double[] data = source.data;
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                destination.setParameterValueQuietly(i, j, data[i * dimension + j]);
            }
        }
        destination.fireParameterChangedEvent();
    }

    private static void setZero(final MatrixParameter matrix) {
        for (int i = 0; i < matrix.getRowDimension(); ++i) {
            for (int j = 0; j < matrix.getColumnDimension(); ++j) {
                matrix.setParameterValueQuietly(i, j, 0.0);
            }
        }
        matrix.fireParameterChangedEvent();
    }

    private static String formatTrait(final double[] values) {
        final StringBuilder builder = new StringBuilder();
        for (int i = 0; i < values.length; ++i) {
            if (i > 0) {
                builder.append(' ');
            }
            builder.append(values[i]);
        }
        return builder.toString();
    }

    public static final XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public String getParserName() {
            return GENERAL_OU_TREE_SIMULATOR;
        }

        @Override
        public Object parseXMLObject(final XMLObject xo) throws XMLParseException {
            final String traitName = xo.getStringAttribute(TRAIT_NAME);
            final boolean cloneTree = xo.getAttribute(CLONE, true);

            final Tree tree = (Tree) xo.getChild(Tree.class);
            final MultivariateDiffusionModel diffusionModel =
                    (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);
            final MatrixParameterInterface selectionMatrix =
                    (MatrixParameterInterface) xo.getChild(STRENGTH_OF_SELECTION_MATRIX).getChild(MatrixParameterInterface.class);
            final Parameter optimumMean =
                    (Parameter) xo.getChild(OPTIMUM_MEAN).getChild(Parameter.class);

            OUSelectionChartParserHelper.validateSelectionChart(
                    xo, selectionMatrix, GENERAL_OU_TREE_SIMULATOR);

            final MultivariateElasticModel elasticModel = new MultivariateElasticModel(selectionMatrix);
            final GeneralOUTreeSimulator simulator;

            final XMLObject rootMeanXO = xo.getChild(ROOT_MEAN);
            if (rootMeanXO != null) {
                if (xo.getChild(ROOT_COVARIANCE) == null) {
                    throw new XMLParseException(ROOT_COVARIANCE + " must be provided when " + ROOT_MEAN + " is used");
                }
                if (xo.getChild(ConjugateRootTraitPrior.class) != null) {
                    throw new XMLParseException("Specify either conjugateRootPrior or explicit rootMean/rootCovariance, not both");
                }
                final Parameter rootMean = (Parameter) rootMeanXO.getChild(Parameter.class);
                final MatrixParameterInterface rootCovariance =
                        (MatrixParameterInterface) xo.getChild(ROOT_COVARIANCE).getChild(MatrixParameterInterface.class);
                simulator = new GeneralOUTreeSimulator(
                        tree, elasticModel, diffusionModel, optimumMean, rootMean, rootCovariance);
            } else {
                if (xo.getChild(ConjugateRootTraitPrior.class) == null) {
                    throw new XMLParseException("Either conjugateRootPrior or explicit rootMean/rootCovariance must be provided");
                }
                final ConjugateRootTraitPrior rootPrior = (ConjugateRootTraitPrior) xo.getChild(ConjugateRootTraitPrior.class);
                simulator = new GeneralOUTreeSimulator(
                        tree, elasticModel, diffusionModel, optimumMean, rootPrior);
            }

            return simulator.simulateAnnotatedTree(traitName, cloneTree);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return RULES;
        }

        @Override
        public String getParserDescription() {
            return "Simulates a constant-optimum multivariate OU trait on a tree using exact branch transitions.";
        }

        @Override
        public Class getReturnType() {
            return Tree.class;
        }
    };

    private static final XMLSyntaxRule[] RULES = new XMLSyntaxRule[]{
            new StringAttributeRule(TRAIT_NAME, "The trait attribute name written onto tips and nodes."),
            AttributeRule.newBooleanRule(CLONE, true),
            new StringAttributeRule(OUSelectionChartParserHelper.SELECTION_CHART,
                    "Selection-matrix chart for OU models. Orthogonal block is the default; dense must be explicit.",
                    OUSelectionChartParserHelper.ALLOWED_SELECTION_CHARTS, true),
            new ElementRule(Tree.class),
            new ElementRule(MultivariateDiffusionModel.class),
            new ElementRule(STRENGTH_OF_SELECTION_MATRIX,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(OPTIMUM_MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(ConjugateRootTraitPrior.class, true),
            new ElementRule(ROOT_MEAN,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, true),
            new ElementRule(ROOT_COVARIANCE,
                    new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}, true)
    };
}
