package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.math.matrixAlgebra.SymmetricMatrix;

import java.util.*;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */
public class IntegratedMultivariateTraitLikelihood extends AbstractMultivariateTraitLikelihood {

    public static final double LOG_SQRT_2_PI = 0.5 * Math.log(2 * Math.PI);

    public IntegratedMultivariateTraitLikelihood(String traitName,
                                                 TreeModel treeModel,
                                                 MultivariateDiffusionModel diffusionModel,
                                                 CompoundParameter traitParameter,
                                                 List<Integer> missingIndices,
                                                 boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                                 BranchRateModel rateModel, Model samplingDensity,
                                                 boolean reportAsMultivariate,
                                                 MultivariateNormalDistribution rootPrior) {
        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate);

        dimTrait = diffusionModel.getPrecisionmatrix().length;
        dim = treeModel.getMultivariateNodeTrait(treeModel.getExternalNode(0), traitName).length;


        if (dim % dimTrait != 0)
            throw new RuntimeException("dim is not divisible by dimTrait");

        dimData = dim / dimTrait;

        meanCache = new double[dim * treeModel.getNodeCount()];
        drawnStates = new double[dim * treeModel.getNodeCount()];
        upperPrecisionCache = new double[dim * treeModel.getNodeCount()];
        lowerPrecisionCache = new double[dim * treeModel.getNodeCount()];
        logRemainderDensityCache = new double[dim * treeModel.getNodeCount()];

        // Set tip data values
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            int index = node.getNumber();
            double[] traitValue = treeModel.getMultivariateNodeTrait(node, traitName);
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        }

        setRootPrior(rootPrior);

        StringBuffer sb = new StringBuffer();
        sb.append("\tDiffusion dimension: ").append(dimTrait).append("\n");
        sb.append("\tNumber of observations: ").append(dimData).append("\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }

    private void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();

        try {
            rootPriorPrecisionDeterminant = new Matrix(rootPriorPrecision).determinant();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
    }

    protected String extraInfo() {
        return "\tSample internal node traits: false\n";
    }

    public double getLogDataLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        areStatesRedrawn = false;
        return logLikelihood;
    }

    public double calculateLogLikelihood() {

        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double logDetTreePrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        postOrderTraverse(treeModel, treeModel.getRoot(), treePrecision, logDetTreePrecision);

        // Integrate root trait out against rootPrior
        final int rootIndex = treeModel.getRoot().getNumber();

        if (DEBUG) {
            System.err.println("mean: " + new Vector(meanCache));
            System.err.println("upre: " + new Vector(upperPrecisionCache));
            System.err.println("lpre: " + new Vector(lowerPrecisionCache));
        }

        double logLikelihood = 0;

        for (int datum = 0; datum < dimData; datum++) {

            double[] rootMean = new double[dimTrait]; // TODO Allocate once
            System.arraycopy(meanCache, rootIndex + datum*dimTrait, rootMean, 0, dimTrait);
            double rootPrecision = lowerPrecisionCache[rootIndex];

            double[][] diffusionPrecision = diffusionModel.getPrecisionmatrix();

            if (DEBUG) {
                System.err.println("root mean: " + new Vector(rootMean));
                System.err.println("root prec: " + rootPrecision);
                System.err.println("diffusion prec: " + new Matrix(diffusionPrecision));
            }
                        
            // B = root prior precision
            // z = root prior mean
            // A = likelihood precision
            // y = likelihood mean

            double[] Bz = new double[dimTrait]; // TODO Allocate once
            double[] Ay = new double[dimTrait]; // TODO Allocate once
            double[][] AplusB = new double[dimTrait][dimTrait]; // TODO Allocate once

            double zBz = 0;
            double yAy = 0;
            double detAplusB = 0;
            double square = 0;

            double detDiffusion = diffusionModel.getDeterminantPrecisionMatrix();

            // y'Ay
            for (int i = 0; i < dimTrait; i++) {
                Ay[i] = 0;
                for (int j = 0; j < dimTrait; j++)
                    Ay[i] += diffusionPrecision[i][j] * rootMean[j] * rootPrecision;
                yAy += rootMean[i] * Ay[i];
            }

            logLikelihood += -LOG_SQRT_2_PI * dimTrait
                    + 0.5 * (Math.log(detDiffusion) + dimTrait * Math.log(rootPrecision))
                    - 0.5 * (yAy);

            if (DEBUG) {
                double[][] T = new double[dimTrait][dimTrait];
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++) {
                        T[i][j] = diffusionPrecision[i][j] * rootPrecision;
                    }
                }

                System.err.println("Conditional root MVN precision = \n" + new Matrix(T));

                System.err.println("Conditional root MVN density = " + MultivariateNormalDistribution.logPdf(rootMean, new double[dimTrait], T,
                        Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(T)), 1.0));
            }

            if (integrateRoot) {

                // z'Bz -- TODO only need to calculate when root prior chances (just at initialization once)
                for (int i = 0; i < dimTrait; i++) {
                    Bz[i] = 0;
                    for (int j = 0; j < dimTrait; j++)
                        Bz[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                    zBz += rootPriorMean[i] * Bz[i];
                }

                // (Ay + Bz)' (A+B)^{-1} (Ay + Bz)
                for (int i = 0; i < dimTrait; i++) {
                    Ay[i] += Bz[i];   // Ay is filled with sum
                    for (int j = 0; j < dimTrait; j++) {
                        AplusB[i][j] = diffusionPrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
                    }
                }

                Matrix mat = new Matrix(AplusB);

                try {
                    detAplusB = mat.determinant();

                } catch (IllegalDimension illegalDimension) {
                    illegalDimension.printStackTrace();
                }

                double[][] invAplusB = mat.inverse().toComponents();
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++)
                        square += Ay[i] * invAplusB[i][j] * Ay[j];
                }
                logLikelihood +=
                    + 0.5 * (Math.log(rootPriorPrecisionDeterminant) - Math.log(detAplusB))
                    - 0.5 * (zBz - square);
            }

            if (DEBUG) {
                System.err.println("zBz = " + zBz);
                System.err.println("yAy = " + yAy);
                System.err.println("(Ay+Bz)(A+B)^{-1}(Ay+Bz) = " + square);
            }
            
        }

        if (DEBUG) {
            System.err.println("logLikelihood (before remainders) = "+logLikelihood+ " (should match root MVN density when root not integrated out)");
        }

        double sumLogRemainders = 0;
        for(double r : logRemainderDensityCache) // TODO Skip remainders for tips and root
            sumLogRemainders += r;

        logLikelihood += sumLogRemainders;

        if (DEBUG && dimTrait == 1) { // Root trait is univariate!!!
            System.err.println("logLikelihood (final) = " + logLikelihood);

            // Perform a check based on filling in the tipCount * tipCount variance matrix
            // And then integrating out the root trait value

            // Form \Sigma (variance) and \Sigma^{-1} (precision)
            double[][] treeVarianceMatrix = computeTreeMVNormalVariance();
            double[][] treePrecisionMatrix = new SymmetricMatrix(treeVarianceMatrix).inverse().toComponents();
            double[] tipTraits = fillLeafTraits();

            System.err.println("tipTraits = " + new Vector(tipTraits));
            System.err.println("tipPrecision = \n" + new Matrix(treePrecisionMatrix));

            double checkLogLikelihood = MultivariateNormalDistribution.logPdf(tipTraits, new double[tipTraits.length], treePrecisionMatrix,
                            Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(treePrecisionMatrix)), 1.0);

            System.err.println("tipDensity = " + checkLogLikelihood + " (should match final likelihood when root not integrated out)");

            // Convolve root prior
            if (integrateRoot) {

                final int tipCount = treePrecisionMatrix.length;

                // 1^t\Sigma^{-1} y + Pz
                double mean = rootPriorPrecision[0][0] * rootPriorMean[0];
                for(int i = 0; i < tipCount; i++) {
                    for(int j = 0; j < tipCount; j++) {
                        mean += treePrecisionMatrix[i][j] * tipTraits[j];
                    }
                }

                // 1^t \Sigma^{-1} 1 + P
                double precision = rootPriorPrecision[0][0];
                for(int i = 0; i < tipCount; i++) {
                    for(int j = 0; j < tipCount; j++) {
                        precision += treePrecisionMatrix[i][j];
                    }
                }
                mean /= precision;


                // We know:  y ~ MVN(x, A) and x ~ N(m, B)
                // Therefore p(x | y) = N( (A+B)^{-1}(Ay + Bm), A + B)
                // We want: p( y ) = p( y | x ) p( x ) / p( x | y ) for any value x, say x = 0

                double logMarginalDensity = checkLogLikelihood;
                logMarginalDensity += MultivariateNormalDistribution.logPdf(
                        rootPriorMean, new double[rootPriorMean.length], rootPriorPrecision,
                        Math.log(rootPriorPrecisionDeterminant),1.0
                );
                logMarginalDensity -= NormalDistribution.logPdf(mean, 0.0, 1.0 / Math.sqrt(precision));
                System.err.println("Mean = "+mean);
                System.err.println("Prec = "+precision);
                System.err.println("log density = "+ logMarginalDensity);
            }
        }

        return logLikelihood;
    }

    void postOrderTraverse(TreeModel treeModel, NodeRef node, double[][] precisionMatrix, double logDetPrecisionMatrix) {

        if (treeModel.isExternal(node)) {
            // Fill in precision scalar, traitValues already filled in
            upperPrecisionCache[node.getNumber()] = 1.0 / getRescaledBranchLength(node);
            lowerPrecisionCache[node.getNumber()] = Double.POSITIVE_INFINITY;
            return;
        }

        final NodeRef childNode0 = treeModel.getChild(node, 0);
        final NodeRef childNode1 = treeModel.getChild(node, 1);

        postOrderTraverse(treeModel, childNode0, precisionMatrix, logDetPrecisionMatrix);
        postOrderTraverse(treeModel, childNode1, precisionMatrix, logDetPrecisionMatrix);

        final int childNumber0 = childNode0.getNumber();
        final int childNumber1 = childNode1.getNumber();
        final int thisNumber = node.getNumber();
        final int meanOffset0 = dim * childNumber0;
        final int meanOffset1 = dim * childNumber1;
        final int meanThisOffset = dim * thisNumber;

        final double precision0 = upperPrecisionCache[childNumber0];
        final double precision1 = upperPrecisionCache[childNumber1];
        final double totalPrecision = precision0 + precision1;

        lowerPrecisionCache[thisNumber] = totalPrecision;

        // Multiple child0 and child1 densities
        for (int i = 0; i < dim; i++)
            meanCache[meanThisOffset + i] = (meanCache[meanOffset0 + i] * precision0 +
                    meanCache[meanOffset1 + i] * precision1)
                    / totalPrecision;

        if (!treeModel.isRoot(node)) {

            // Integrate out trait value at this node
            double thisPrecision = 1.0 / getRescaledBranchLength(node);
            upperPrecisionCache[thisNumber] = totalPrecision * thisPrecision / (totalPrecision + thisPrecision);
        }

        // Compute logRemainderDensity

        logRemainderDensityCache[thisNumber] = 0;

        for(int k=0; k<dimData; k++) {

            double childSS0 = 0;
            double childSS1 = 0;
            double crossSS  = 0;

            final double remainderPrecision = precision0 * precision1 / (precision0 + precision1);

            for(int i=0; i<dimTrait; i++) {

                final double wChild0i = meanCache[meanOffset0+i] * precision0;
                final double wChild1i = meanCache[meanOffset1+i] * precision1;

                for(int j=0; j<dimTrait; j++) {

                    final double child0j = meanCache[meanOffset0+j];
                    final double child1j = meanCache[meanOffset1+j];

                    childSS0 += wChild0i * precisionMatrix[i][j] * child0j;
                    childSS1 += wChild1i * precisionMatrix[i][j] * child1j;

                    crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * meanCache[meanThisOffset + j];
                }

                logRemainderDensityCache[thisNumber] +=
                        - dimTrait * LOG_SQRT_2_PI
                        + 0.5 * dimTrait * (Math.log(remainderPrecision)+logDetPrecisionMatrix)
                        - 0.5 * (childSS0 + childSS1 - crossSS);
            }
        }
    }

    class MeanPrecision {
        MeanPrecision(double[] mean, double precisionScalar) {
            this.mean = mean;
            this.precisionScalar = precisionScalar;
        }

        MeanPrecision(double[] mean, double[][] precision) {
            this.mean = mean;
            this.precision = precision;
        }

        double[] mean;
        double precisionScalar;
        double[][] precision;
    }

    protected double[] traitForNode(TreeModel tree, NodeRef node, String traitName) {
        // TODO Must to a pre-order traversal to draw values, see AncestralStateTreeLikelihood
        if (tree != treeModel) {
             throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!areStatesRedrawn)
            redrawAncestralStates();

        int index = node.getNumber();

        return new double[dimTrait];
    }

    protected double[][] computeTreeMVNormalVariance() {

        if (dimTrait != 1) {
            throw new RuntimeException("Currently only computed for 1D traits");
        }
        final int tipCount = treeModel.getExternalNodeCount();
        double[][] variance = new double[dimTrait * tipCount][dim * tipCount];
        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double precision = treePrecision[0][0];

        for (int i = 0; i < tipCount; i++ ) {
            variance[i][i] = getRescaledLengthToRoot(treeModel.getExternalNode(i)) / precision;
            for (int j = i+1; j < tipCount; j++) {
                Set<String> leafNames = new HashSet<String>();
                leafNames.add(treeModel.getTaxonId(i));
                leafNames.add(treeModel.getTaxonId(j));
                NodeRef mrca = Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
                variance[j][i] = variance[i][j] = getRescaledLengthToRoot(mrca) / precision;
            }
        }
        System.err.println("");
        System.err.println("Conditional variance:\n"+new Matrix(variance));
        return variance;
    }

    protected double[] fillLeafTraits() {
        if (dimTrait != 1) {
            throw new RuntimeException("Currently only computed for 1D traits");
        }

        final int tipCount = treeModel.getExternalNodeCount();
        double[] traits = new double[dimTrait * tipCount];
        System.arraycopy(meanCache, 0, traits, 0, dimTrait * tipCount);
        return traits;
    }

    private double getRescaledLengthToRoot(NodeRef node) {
        double length = 0;
        final NodeRef root = treeModel.getRoot();
        while( node != root ) {
            length += getRescaledBranchLength(node);
            node = treeModel.getParent(node);
        }
        return length;
    }

     private boolean areStatesRedrawn = false;

     public void redrawAncestralStates() {
         preOrderTraverseSample(treeModel, treeModel.getRoot(), 0, diffusionModel.getPrecisionmatrix());
         areStatesRedrawn = true;
     }

    void preOrderTraverseSample(TreeModel treeModel, NodeRef node, int parentIndex, double[][] precisionMatrix) {

        final int thisIndex = node.getNumber();

        if (treeModel.isRoot(node)) {
            // draw root
        } else { // draw conditional on parentState

            final int thisOffset = thisIndex * dimTrait;

            if (Double.isInfinite(lowerPrecisionCache[thisIndex]))
                System.arraycopy(meanCache,thisOffset,drawnStates,thisOffset,dimTrait);
            else {
                final int parentOffset = parentIndex * dimTrait;
                // parent trait at drawnStates[parentOffset]
            }
        }

        if (!treeModel.isExternal(node)) {
            preOrderTraverseSample(treeModel, treeModel.getChild(node,0), thisIndex, precisionMatrix);
            preOrderTraverseSample(treeModel, treeModel.getChild(node,1), thisIndex, precisionMatrix);
        }
    }

    private double[] meanCache;
    private double[] upperPrecisionCache;
    private double[] lowerPrecisionCache;
    private double[] logRemainderDensityCache;

    private double[] drawnStates;

    private int dimData;
    private int dimTrait;
    private int dim;

    private double[] rootPriorMean;
    private double[][] rootPriorPrecision;
    private double rootPriorPrecisionDeterminant;

    private final boolean integrateRoot = true;
    private static boolean DEBUG = false;
}
