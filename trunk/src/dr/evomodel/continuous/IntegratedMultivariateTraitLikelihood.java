package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.math.distributions.MultivariateNormalDistribution;
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
        upperPrecisionCache = new double[treeModel.getNodeCount()];
        lowerPrecisionCache = new double[treeModel.getNodeCount()];
        logRemainderDensityCache = new double[treeModel.getNodeCount()];

        // Set up reusable temporary storage
        Bz = new double[dimTrait];
        Ay = new double[dimTrait];
        tmpM = new double[dimTrait][dimTrait];
        tmp2 = new double[dimTrait];

        setRootPrior(rootPrior);

        // Set tip data values
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            int index = node.getNumber();
            double[] traitValue = treeModel.getMultivariateNodeTrait(node, traitName);
            System.arraycopy(traitValue, 0, meanCache, dim * index, dim);
        }

        StringBuffer sb = new StringBuffer();
        sb.append("\tDiffusion dimension: ").append(dimTrait).append("\n");
        sb.append("\tNumber of observations: ").append(dimData).append("\n");
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }

    private void setRootPriorSumOfSquares() {
        if (integrateRoot) {
            // z'Bz -- sum-of-squares root contribution
            for (int i = 0; i < dimTrait; i++) {
                Bz[i] = 0;
                for (int j = 0; j < dimTrait; j++)
                    Bz[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
                zBz += rootPriorMean[i] * Bz[i];
            }
        } else {
            zBz = 0;
        }
    }

    private void setRootPrior(MultivariateNormalDistribution rootPrior) {
        rootPriorMean = rootPrior.getMean();
        rootPriorPrecision = rootPrior.getScaleMatrix();

        try {
            logRootPriorPrecisionDeterminant = Math.log(new Matrix(rootPriorPrecision).determinant());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        setRootPriorSumOfSquares();
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

        double logLikelihood = 0;
        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double logDetTreePrecision = Math.log(diffusionModel.getDeterminantPrecisionMatrix());
        double[] rootMean = tmp2;

        // Use dynamic programming to compute conditional likelihoods at each internal node
        postOrderTraverse(treeModel, treeModel.getRoot(), treePrecision, logDetTreePrecision);

        if (DEBUG) {
            System.err.println("mean: " + new Vector(meanCache));
            System.err.println("upre: " + new Vector(upperPrecisionCache));
            System.err.println("lpre: " + new Vector(lowerPrecisionCache));
            System.err.println("cach: " + new Vector(logRemainderDensityCache));
        }

        // Compute the contribution of each datum at the root
        final int rootIndex = treeModel.getRoot().getNumber();
        double rootPrecision = lowerPrecisionCache[rootIndex];
        for (int datum = 0; datum < dimData; datum++) {

            double thisLogLikelihood = 0;
            System.arraycopy(meanCache, rootIndex * dim + datum * dimTrait, rootMean, 0, dimTrait);

            if (DEBUG) {
                System.err.println("Datum #" + datum);
                System.err.println("root mean: " + new Vector(rootMean));
                System.err.println("root prec: " + rootPrecision);
                System.err.println("diffusion prec: " + new Matrix(treePrecision));
            }

            // B = root prior precision
            // z = root prior mean
            // A = likelihood precision
            // y = likelihood mean

            // y'Ay
            double yAy = determineSumOfSquares(rootMean, Ay, treePrecision, rootPrecision); // Fills in Ay

            thisLogLikelihood += -LOG_SQRT_2_PI * dimTrait
                    + 0.5 * (logDetTreePrecision + dimTrait * Math.log(rootPrecision) - yAy);

            if (DEBUG) {
                double[][] T = new double[dimTrait][dimTrait];
                for (int i = 0; i < dimTrait; i++) {
                    for (int j = 0; j < dimTrait; j++) {
                        T[i][j] = treePrecision[i][j] * rootPrecision;
                    }
                }
                System.err.println("Conditional root MVN precision = \n" + new Matrix(T));
                System.err.println("Conditional root MVN density = " + MultivariateNormalDistribution.logPdf(rootMean, new double[dimTrait], T,
                        Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(T)), 1.0));
            }

            if (integrateRoot) {
                // Integrate root trait out against rootPrior
                thisLogLikelihood += integrateLogLikelihoodAtRoot(Ay, Bz, tmpM, treePrecision, rootPrecision); // Ay is destroyed
            }

            if (DEBUG) {
                System.err.println("zBz = " + zBz);
                System.err.println("yAy = " + yAy);
                System.err.println("logLikelihood (before remainders) = " + thisLogLikelihood + " (should match conditional root MVN density when root not integrated out)");
            }

            logLikelihood += thisLogLikelihood;
        }

        logLikelihood += sumLogRemainders();

        if (DEBUG) { // Root trait is univariate!!!
            System.err.println("logLikelihood (final) = " + logLikelihood);
//            if (dimTrait == 1) {
            checkViaLargeMatrixInversion();
//            }
        }

        return logLikelihood;
    }

    private double determineSumOfSquares(double[] y, double[] Ay, double[][] A, double scale) {
        // returns Ay and yAy
        double yAy = 0;
        for (int i = 0; i < dimTrait; i++) {
            Ay[i] = 0;
            for (int j = 0; j < dimTrait; j++)
                Ay[i] += A[i][j] * y[j] * scale;
            yAy += y[i] * Ay[i];
        }
        return yAy;
    }

    private double sumLogRemainders() {
        double sumLogRemainders = 0;
        for (double r : logRemainderDensityCache)
            sumLogRemainders += r;
        // Could skip leafs
        return sumLogRemainders;
    }


//    private MeanPrecision convolve(double[] Ay, double[] Bz) {
//               for (int i = 0; i < dimTrait; i++) {
//                Ay[i] += Bz[i];   // Ay is filled with sum, and original value is destroyed
//                for (int j = 0; j < dimTrait; j++) {
//                    AplusB[i][j] = treePrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
//                }
//            }
//
//            Matrix mat = new Matrix(AplusB);
//
//            try {
//                detAplusB = mat.determinant();
//
//            } catch (IllegalDimension illegalDimension) {
//                illegalDimension.printStackTrace();
//            }
//
//            double[][] invAplusB = mat.inverse().toComponents();
//    }

    private double integrateLogLikelihoodAtRoot(double[] Ay, double[] Bz,
                                                double[][] AplusB,
                                                double[][] treePrecision, double rootPrecision) {
        double detAplusB = 0;
        double square = 0;

        // square : (Ay + Bz)' (A+B)^{-1} (Ay + Bz)

        if (dimTrait > 1) {
            for (int i = 0; i < dimTrait; i++) {
                Ay[i] += Bz[i];   // Ay is filled with sum, and original value is destroyed
                for (int j = 0; j < dimTrait; j++) {
                    AplusB[i][j] = treePrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
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
        } else {
            // 1D is very simple
            detAplusB = treePrecision[0][0] * rootPrecision + rootPriorPrecision[0][0];
            Ay[0] += Bz[0];
            square = Ay[0] * Ay[0] / detAplusB;
        }

        if (DEBUG) {
            System.err.println("(Ay+Bz)(A+B)^{-1}(Ay+Bz) = " + square);
        }

        return 0.5 * (logRootPriorPrecisionDeterminant - Math.log(detAplusB) - zBz + square);
    }

    private void checkViaLargeMatrixInversion() {

        // Perform a check based on filling in the (dimTrait * tipCount) * (dimTrait * tipCount) variance matrix
        // And then integrating out the root trait value

        // Form \Sigma (variance) and \Sigma^{-1} (precision)
        double[][] treeTraitVarianceMatrix = computeTreeMVNormalVariance();
        double[][] treeTraitPrecisionMatrix = new SymmetricMatrix(treeTraitVarianceMatrix).inverse().toComponents();
        double totalLogDensity = 0;

        for (int datum = 0; datum < dimData; datum++) {

            double[] tipTraits = fillLeafTraits(datum);

            System.err.println("Datum #" + datum);
            System.err.println("tipTraits = " + new Vector(tipTraits));
            System.err.println("tipPrecision = \n" + new Matrix(treeTraitPrecisionMatrix));

            double checkLogLikelihood = MultivariateNormalDistribution.logPdf(tipTraits, new double[tipTraits.length], treeTraitPrecisionMatrix,
                    Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(treeTraitPrecisionMatrix)), 1.0);

            System.err.println("tipDensity = " + checkLogLikelihood + " (should match final likelihood when root not integrated out)");

            // Convolve root prior
            if (integrateRoot) {
                checkLogLikelihood += integrateLogLikelihoodAtRootFromFullTreeMatrix(treeTraitPrecisionMatrix, tipTraits);
            }
            totalLogDensity += checkLogLikelihood;
        }
        System.err.println("Total logLikelihood (via tree) = " + totalLogDensity);
    }

    private double integrateLogLikelihoodAtRootFromFullTreeMatrix(double[][] treeTraitPrecisionMatrix, double[] tipTraits) {

        double logLikelihood = 0;
        final int tipCount = treeModel.getExternalNodeCount();

        // 1^t\Sigma^{-1} y + Pz
        double[] mean = Ay;
        for (int i = 0; i < dimTrait; i++) {
            mean[i] = 0;
            for (int j = 0; j < dimTrait; j++) {
                mean[i] += rootPriorPrecision[i][j] * rootPriorMean[j];
            }

            for (int j = 0; j < tipCount; j++) {
                final int rowOffset = j * dimTrait + i;
                for (int k = 0; k < tipCount * dimTrait; k++) {
//                            System.err.println("rowOffset = "+rowOffset);
//                            System.err.println("k = "+k);
//                            System.err.println("tipC = "+tipCount);
                    mean[i] += treeTraitPrecisionMatrix[rowOffset][k] * tipTraits[k];
                }
            }
        }

        // 1^t \Sigma^{-1} 1 + P
        double[][] precision = tmpM;
        for (int i = 0; i < dimTrait; i++) {
            for (int j = 0; j < dimTrait; j++) {
                precision[i][j] = rootPriorPrecision[i][j];
                for (int k = 0; k < tipCount; k++) {
                    for (int l = 0; l < tipCount; l++) {
                        precision[i][j] += treeTraitPrecisionMatrix[k * dimTrait + i][l * dimTrait + j];
                    }
                }
            }
        }
        double[] normalizedMean = tmp2;
        double[][] variance = new SymmetricMatrix(precision).inverse().toComponents();
        for (int i = 0; i < dimTrait; i++) {
            normalizedMean[i] = 0.0;
            for (int j = 0; j < dimTrait; j++) {
                normalizedMean[i] += variance[i][j] * mean[j];
            }
        }
        mean = normalizedMean;

        // We know:  y ~ MVN(x, A) and x ~ N(m, B)
        // Therefore p(x | y) = N( (A+B)^{-1}(Ay + Bm), A + B)
        // We want: p( y ) = p( y | x ) p( x ) / p( x | y ) for any value x, say x = 0

        logLikelihood += MultivariateNormalDistribution.logPdf(
                rootPriorMean, new double[rootPriorMean.length], rootPriorPrecision,
                logRootPriorPrecisionDeterminant, 1.0
        );

        logLikelihood -= MultivariateNormalDistribution.logPdf(
                mean, new double[mean.length], precision,
                Math.log(MultivariateNormalDistribution.calculatePrecisionMatrixDeterminate(precision)), 1.0
        );
//        logLikelihood -= NormalDistribution.logPdf(mean, 0.0, 1.0 / Math.sqrt(precision));
        System.err.println("Mean = " + new Vector(mean));
        System.err.println("Prec = " + new Matrix(precision));
        System.err.println("log density = " + logLikelihood);
        return logLikelihood;
    }

    public void makeDirty() {
        super.makeDirty();
        areStatesRedrawn = false;
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
        final double remainderPrecision = precision0 * precision1 / (precision0 + precision1);

        for (int k = 0; k < dimData; k++) {

            double childSS0 = 0;
            double childSS1 = 0;
            double crossSS = 0;

            for (int i = 0; i < dimTrait; i++) {

                final double wChild0i = meanCache[meanOffset0 + k * dimTrait + i] * precision0;
                final double wChild1i = meanCache[meanOffset1 + k * dimTrait + i] * precision1;

                for (int j = 0; j < dimTrait; j++) {

                    final double child0j = meanCache[meanOffset0 + k * dimTrait + j];
                    final double child1j = meanCache[meanOffset1 + k * dimTrait + j];

                    childSS0 += wChild0i * precisionMatrix[i][j] * child0j;
                    childSS1 += wChild1i * precisionMatrix[i][j] * child1j;

                    crossSS += (wChild0i + wChild1i) * precisionMatrix[i][j] * meanCache[meanThisOffset + k * dimTrait + j];
                }
            }

            logRemainderDensityCache[thisNumber] +=
                    -dimTrait * LOG_SQRT_2_PI
                            + 0.5 * (dimTrait * Math.log(remainderPrecision) + logDetPrecisionMatrix)
                            - 0.5 * (childSS0 + childSS1 - crossSS);
        }
    }

    protected double[] getRootNodeTrait() {
        return new double[dim];
    }

    protected double[] traitForNode(TreeModel tree, NodeRef node, String traitName) {

        if (tree != treeModel) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!areStatesRedrawn)
            redrawAncestralStates();

        int index = node.getNumber();

        double[] trait = new double[dim];
        System.arraycopy(drawnStates, index * dim, trait, 0, dim);
        return trait;
    }

//    protected double[][] computeTreeMVNormalVariance() {
//
//        final int tipCount = treeModel.getExternalNodeCount();
//        double[][] variance = new double[tipCount][tipCount];
//        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
//        double precision = treePrecision[0][0];
//
//        for (int i = 0; i < tipCount; i++ ) {
//            variance[i][i] = getRescaledLengthToRoot(treeModel.getExternalNode(i)) / precision;
//            for (int j = i+1; j < tipCount; j++) {
//                Set<String> leafNames = new HashSet<String>();
//                leafNames.add(treeModel.getTaxonId(i));
//                leafNames.add(treeModel.getTaxonId(j));
//                NodeRef mrca = Tree.Utils.getCommonAncestorNode(treeModel, leafNames);
//                variance[j][i] = variance[i][j] = getRescaledLengthToRoot(mrca) / precision;
//            }
//        }
//        System.err.println("");
//        System.err.println("Conditional variance:\n"+new Matrix(variance));
//        return variance;
//    }

    protected double[][] computeTreeMVNormalVariance() {

        final int tipCount = treeModel.getExternalNodeCount();
        double[][] variance = new double[dimTrait * tipCount][dimTrait * tipCount];
        double[][] traitPrecision = diffusionModel.getPrecisionmatrix();
        double[][] traitVariance = new SymmetricMatrix(traitPrecision).inverse().toComponents();
        System.err.println("Trait variance:\n" + new Matrix(traitVariance));

        for (int i = 0; i < tipCount; i++) {
            int rowOffset = i * dimTrait;

            // Fill in diagonal block
            double marginalTime = getRescaledLengthToRoot(treeModel.getExternalNode(i));
            for (int k = 0; k < dimTrait; k++) {
                for (int l = k; l < dimTrait; l++) {
                    variance[rowOffset + k][rowOffset + l] = traitVariance[k][l] * marginalTime;
                }
            }

            for (int j = i + 1; j < tipCount; j++) {
                Set<String> leafNames = new HashSet<String>();
                leafNames.add(treeModel.getTaxonId(i));
                leafNames.add(treeModel.getTaxonId(j));
                NodeRef mrca = Tree.Utils.getCommonAncestorNode(treeModel, leafNames);

                // Fill in off-diagonal block
                double sharedTime = getRescaledLengthToRoot(mrca);
                if (sharedTime > 0) {
                    int colOffset = j * dimTrait;
                    for (int k = 0; k < dimTrait; k++) {
                        for (int l = 0; l < dimTrait; l++) {
                            variance[rowOffset + k][colOffset + l] = traitVariance[k][l] * sharedTime;
                        }
                    }
                }
            }
        }

        // Make symmetric
        for (int i = 0; i < dimTrait * tipCount; i++) {
            for (int j = i + 1; j < dimTrait * tipCount; j++) {
                variance[j][i] = variance[i][j];
            }
        }
        System.err.println("");
        System.err.println("Conditional variance:\n" + new Matrix(variance));
        return variance;
    }


    protected double[] fillLeafTraits(int datum) {

        final int tipCount = treeModel.getExternalNodeCount();
        double[] traits = new double[dimTrait * tipCount];
        int index = 0;
        for (int i = 0; i < tipCount; i++) {
            for (int k = 0; k < dimTrait; k++) {
                traits[index++] = meanCache[dim * i + datum * dimTrait + k];
            }
        }
        return traits;
    }

    private double getRescaledLengthToRoot(NodeRef node) {
        double length = 0;
        final NodeRef root = treeModel.getRoot();
        while (node != root) {
            length += getRescaledBranchLength(node);
            node = treeModel.getParent(node);
        }
        return length;
    }

    private boolean areStatesRedrawn = false;

    public void redrawAncestralStates() {

        double[][] treePrecision = diffusionModel.getPrecisionmatrix();
        double[][] treeVariance = new SymmetricMatrix(treePrecision).inverse().toComponents();

        preOrderTraverseSample(treeModel, treeModel.getRoot(), 0, treePrecision, treeVariance);

        if (DEBUG) {
            System.err.println("all draws = " + new Vector(drawnStates));
        }

        areStatesRedrawn = true;
    }

    void preOrderTraverseSample(TreeModel treeModel, NodeRef node, int parentIndex, double[][] treePrecision, double[][] treeVariance) {

        final int thisIndex = node.getNumber();
        double[][] AplusB = tmpM;

        if (treeModel.isRoot(node)) {
            // draw root

            double[] rootMean = new double[dimTrait];
            final int rootIndex = treeModel.getRoot().getNumber();
            double rootPrecision = lowerPrecisionCache[rootIndex];

            for (int datum = 0; datum < dimData; datum++) {
                System.arraycopy(meanCache, thisIndex * dim + datum * dimTrait, rootMean, 0, dimTrait);
                determineSumOfSquares(rootMean, Ay, treePrecision, rootPrecision); // Fills in Ay

                for (int i = 0; i < dimTrait; i++) {
                    Ay[i] += Bz[i];   // Ay is filled with sum, and original value is destroyed
                    for (int j = 0; j < dimTrait; j++) {
                        AplusB[i][j] = treePrecision[i][j] * rootPrecision + rootPriorPrecision[i][j];
                    }
                }
                Matrix mat = new Matrix(AplusB);
                double[][] invAplusB = mat.inverse().toComponents();

                // Expected value: (A + B)^{-1}(Ay + Bz)
                for (int i = 0; i < dimTrait; i++) {
                    rootMean[i] = 0.0;
                    for (int j = 0; j < dimTrait; j++) {
                        rootMean[i] += invAplusB[i][j] * Ay[j];
                    }
                }

                double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(rootMean, invAplusB);
                System.arraycopy(draw, 0, drawnStates, rootIndex * dim + datum * dimTrait, dimTrait);

                if (DEBUG) {
                    System.err.println("Root mean: " + new Vector(rootMean));
                    System.err.println("Root prec: " + new Matrix(AplusB));
                    System.err.println("Root draw: " + new Vector(draw));
                }
            }
        } else { // draw conditional on parentState

            if (Double.isInfinite(lowerPrecisionCache[thisIndex])) {

                System.arraycopy(meanCache, thisIndex * dim, drawnStates, thisIndex * dim, dim);

            } else {
                // parent trait at drawnStates[parentOffset]
                double precisionToParent = 1.0 / getRescaledBranchLength(node);
                double precisionOfNode = lowerPrecisionCache[thisIndex];
                double totalPrecision = precisionOfNode + precisionToParent;

                double[] mean = Ay; // temporary storage
                double[][] var = tmpM; // temporary storage

                for (int datum = 0; datum < dimData; datum++) {

                    int parentOffset = parentIndex * dim + datum * dimTrait;
                    int thisOffset = thisIndex * dim + datum * dimTrait;

                    for (int i = 0; i < dimTrait; i++) {
                        mean[i] = (meanCache[parentOffset + i] * precisionToParent
                                + meanCache[thisIndex + i]) / totalPrecision;
                        for (int j = 0; j < dimTrait; j++) {
                            var[i][j] = treeVariance[i][j] / totalPrecision;
                        }
                    }
                    double[] draw = MultivariateNormalDistribution.nextMultivariateNormalVariance(mean, var);
                    System.arraycopy(draw, 0, drawnStates, thisOffset, dimTrait);

                    if (DEBUG) {
                        System.err.println("Int mean: " + new Vector(mean));
                        System.err.println("Int var : " + new Matrix(var));
                        System.err.println("Int draw: " + new Vector(draw));
                    }
                }
            }
        }

        if (!treeModel.isExternal(node)) {
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 0), thisIndex, treePrecision, treeVariance);
            preOrderTraverseSample(treeModel, treeModel.getChild(node, 1), thisIndex, treePrecision, treeVariance);
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
    private double logRootPriorPrecisionDeterminant;

    private final boolean integrateRoot = true;
    private static boolean DEBUG = false;

    private double zBz; // Prior sum-of-squares contribution

    // Reusable temporary storage
    private double[] Bz;
    private double[] Ay;
    private double[][] tmpM;
    private double[] tmp2;

}
