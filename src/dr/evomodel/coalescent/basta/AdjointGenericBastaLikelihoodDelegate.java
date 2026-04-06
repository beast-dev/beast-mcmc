package dr.evomodel.coalescent.basta;

import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;

import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

public class AdjointGenericBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private static final double EIGEN_TOLERANCE = 1.0e-12;

    private final Map<Integer, double[]> tipPartials = new LinkedHashMap<>();
    private final double[] sizes;
    private final EigenDecomposition[] decompositions;
    private ForwardWorkspace forwardWorkspace;
    private ReverseWorkspace reverseWorkspace;
    private ForwardResult cachedForward;
    private List<BranchIntervalOperation> cachedBranchOperations;
    private List<TransitionMatrixOperation> cachedMatrixOperations;
    private List<Integer> cachedIntervalStarts;

    public AdjointGenericBastaLikelihoodDelegate(String name,
                                                 Tree tree,
                                                 int stateCount,
                                                 boolean transpose) {
        super(name, tree, stateCount, transpose);
        this.sizes = new double[2 * stateCount];
        this.decompositions = new EigenDecomposition[1];
    }

    @Override
    public double calculateLikelihood(List<BranchIntervalOperation> branchOperations,
                                      List<TransitionMatrixOperation> matrixOperations,
                                      List<Integer> intervalStarts,
                                      int rootNodeNumber) {
        ensureSupported(matrixOperations);
        return getOrCreateForward(branchOperations, matrixOperations, intervalStarts).logLikelihood;
    }

    @Override
    public double[][] calculateGradient(List<BranchIntervalOperation> branchOperations,
                                        List<TransitionMatrixOperation> matrixOperations,
                                        List<Integer> intervalStarts,
                                        int rootNodeNumber) {
        ensureSupported(matrixOperations);
        ForwardResult forward = getOrCreateForward(branchOperations, matrixOperations, intervalStarts);
        return reverse(branchOperations, matrixOperations, intervalStarts, forward).rateGradient;
    }

    @Override
    public double[] calculateGradientPopSize(List<BranchIntervalOperation> branchOperations,
                                             List<TransitionMatrixOperation> matrixOperations,
                                             List<Integer> intervalStarts,
                                             int rootNodeNumber) {
        ensureSupported(matrixOperations);
        ForwardResult forward = getOrCreateForward(branchOperations, matrixOperations, intervalStarts);
        return reverse(branchOperations, matrixOperations, intervalStarts, forward).populationSizeGradient;
    }

    private ForwardResult getOrCreateForward(List<BranchIntervalOperation> branchOperations,
                                             List<TransitionMatrixOperation> matrixOperations,
                                             List<Integer> intervalStarts) {
        if (cachedForward != null &&
                cachedBranchOperations == branchOperations &&
                cachedMatrixOperations == matrixOperations &&
                cachedIntervalStarts == intervalStarts) {
            return cachedForward;
        }
        ForwardResult forward = forward(branchOperations, matrixOperations, intervalStarts);
        cachedForward = forward;
        cachedBranchOperations = branchOperations;
        cachedMatrixOperations = matrixOperations;
        cachedIntervalStarts = intervalStarts;
        return forward;
    }

    private void invalidateForwardCache() {
        cachedForward = null;
        cachedBranchOperations = null;
        cachedMatrixOperations = null;
        cachedIntervalStarts = null;
    }

    private void ensureSupported(List<TransitionMatrixOperation> matrixOperations) {
        for (TransitionMatrixOperation operation : matrixOperations) {
            if (operation.decompositionBuffer != 0) {
                throw new IllegalArgumentException("Adjoint route currently supports only decompositionBuffer == 0");
            }
            EigenDecomposition decomposition = decompositions[operation.decompositionBuffer];
            if (decomposition == null) {
                throw new IllegalStateException("Missing eigen decomposition for adjoint route");
            }
            double[] eigenValues = decomposition.getEigenValues();
            if (eigenValues.length > stateCount) {
                for (int i = 0; i < stateCount; ++i) {
                    if (Math.abs(eigenValues[stateCount + i]) > EIGEN_TOLERANCE) {
                        throw new IllegalArgumentException("Adjoint route currently requires real eigenvalues");
                    }
                }
            }
        }
    }

    private ForwardResult forward(List<BranchIntervalOperation> branchOperations,
                                  List<TransitionMatrixOperation> matrixOperations,
                                  List<Integer> intervalStarts) {
        int bufferCount = maxBufferIndex(branchOperations) + 1;
        int matrixCount = maxMatrixIndex(matrixOperations) + 1;
        int intervalCount = intervalStarts.size() - 1;
        int maxIntervalNumber = maxIntervalNumber(branchOperations);
        ForwardWorkspace workspace = ensureForwardWorkspace(bufferCount, matrixCount, intervalCount, maxIntervalNumber);
        workspace.reset(intervalCount, maxIntervalNumber);

        double[][] partials = workspace.partials;
        for (Map.Entry<Integer, double[]> entry : tipPartials.entrySet()) {
            copyInto(partials[entry.getKey()], entry.getValue());
        }

        double[][][] matrices = workspace.matrices;
        double[][][] spectralKernels = workspace.spectralKernels;
        for (TransitionMatrixOperation operation : matrixOperations) {
            computeTransitionProbabilitiesInPlace(operation.time, decompositions[operation.decompositionBuffer],
                    matrices[operation.outputBuffer], workspace.transitionScratch);
            computeLoewnerKernelInPlace(operation.time, decompositions[operation.decompositionBuffer], spectralKernels[operation.outputBuffer]);
        }

        double[] coalescentProb = workspace.coalescentProb;
        boolean[] hasCoalescentProb = workspace.hasCoalescentProb;
        CoalescentDetails[] coalescentDetails = workspace.coalescentDetails;

        for (BranchIntervalOperation operation : branchOperations) {
            double[] left = partials[operation.accBuffer1];
            matVecInPlace(matrices[operation.inputMatrix1], partials[operation.inputBuffer1], left);
            if (operation.inputBuffer2 < 0) {
                copyInto(partials[operation.outputBuffer], left);
                continue;
            }

            double[] right = partials[operation.accBuffer2];
            matVecInPlace(matrices[operation.inputMatrix2], partials[operation.inputBuffer2], right);

            CoalescentDetails details = coalescentDetails[operation.outputBuffer];
            double[] unnormalized = details.w;
            double probability = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                unnormalized[i] = left[i] * right[i] / sizes[i];
                probability += unnormalized[i];
            }

            double[] normalized = partials[operation.outputBuffer];
            for (int i = 0; i < stateCount; ++i) {
                normalized[i] = unnormalized[i] / probability;
            }
            coalescentProb[operation.intervalNumber] = probability;
            hasCoalescentProb[operation.intervalNumber] = true;
            details.leftStart = partials[operation.inputBuffer1];
            details.rightStart = partials[operation.inputBuffer2];
            details.leftEnd = left;
            details.rightEnd = right;
            details.j = probability;
        }

        double[][] intervalE = workspace.intervalE;
        double[][] intervalF = workspace.intervalF;
        double[][] intervalG = workspace.intervalG;
        double[][] intervalH = workspace.intervalH;
        double logLikelihood = 0.0;

        for (int interval = 0; interval < intervalCount; ++interval) {
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);
            double[] e = intervalE[interval];
            double[] f = intervalF[interval];
            double[] g = intervalG[interval];
            double[] h = intervalH[interval];

            for (int idx = start; idx < end; ++idx) {
                BranchIntervalOperation operation = branchOperations.get(idx);
                accumulateSummary(partials[operation.inputBuffer1], e, f);
                if (operation.inputBuffer2 >= 0) {
                    accumulateSummary(partials[operation.inputBuffer2], e, f);
                }
                accumulateSummary(partials[operation.accBuffer1], g, h);
                if (operation.accBuffer2 >= 0) {
                    accumulateSummary(partials[operation.accBuffer2], g, h);
                }
            }

            double hazardSum = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                hazardSum += (e[i] * e[i] - f[i] + g[i] * g[i] - h[i]) / sizes[i];
            }

            BranchIntervalOperation anchor = branchOperations.get(start);
            logLikelihood += -anchor.intervalLength * hazardSum / 4.0;
            if (hasCoalescentProb[anchor.intervalNumber]) {
                logLikelihood += Math.log(coalescentProb[anchor.intervalNumber]);
            }
        }

        return new ForwardResult(logLikelihood, partials, matrices, intervalE, intervalF, intervalG, intervalH,
                coalescentProb, hasCoalescentProb, coalescentDetails, spectralKernels);
    }

    private ReverseResult reverse(List<BranchIntervalOperation> branchOperations,
                                  List<TransitionMatrixOperation> matrixOperations,
                                  List<Integer> intervalStarts,
                                  ForwardResult forward) {
        ReverseWorkspace workspace = ensureReverseWorkspace(forward.partials.length, forward.matrices.length);
        workspace.reset(matrixOperations, forward.partials);

        double[][] partialAdjoint = workspace.partialAdjoint;
        double[][][] matrixAdjoint = workspace.matrixAdjoint;
        double[][] rateGradient = workspace.rateGradient;
        double[] populationSizeGradient = workspace.populationSizeGradient;
        double[] adjE = workspace.adjE;
        double[] adjF = workspace.adjF;
        double[] adjG = workspace.adjG;
        double[] adjH = workspace.adjH;
        double[] zBar = workspace.zBar;
        double[] adjW = workspace.adjW;
        double[] leftEndBar = workspace.leftEndBar;
        double[] rightEndBar = workspace.rightEndBar;

        for (int interval = intervalStarts.size() - 2; interval >= 0; --interval) {
            double[] e = forward.intervalE[interval];
            double[] f = forward.intervalF[interval];
            double[] g = forward.intervalG[interval];
            double[] h = forward.intervalH[interval];

            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);
            double length = branchOperations.get(start).intervalLength;

            for (int i = 0; i < stateCount; ++i) {
                double invSize = 1.0 / sizes[i];
                adjE[i] = -length * e[i] / (2.0 * sizes[i]);
                adjF[i] = length * invSize / 4.0;
                adjG[i] = -length * g[i] / (2.0 * sizes[i]);
                adjH[i] = length * invSize / 4.0;
                populationSizeGradient[i] += length * (e[i] * e[i] - f[i] + g[i] * g[i] - h[i]) / (4.0 * sizes[i] * sizes[i]);
            }

            for (int idx = start; idx < end; ++idx) {
                BranchIntervalOperation operation = branchOperations.get(idx);
                pushReductionAdjoints(forward.partials, partialAdjoint, operation.accBuffer1, adjG, adjH);
                pushReductionAdjoints(forward.partials, partialAdjoint, operation.accBuffer2, adjG, adjH);
                pushReductionAdjoints(forward.partials, partialAdjoint, operation.inputBuffer1, adjE, adjF);
                pushReductionAdjoints(forward.partials, partialAdjoint, operation.inputBuffer2, adjE, adjF);
            }

            BranchIntervalOperation coalescentOperation = firstCoalescentInInterval(branchOperations, start, end);
            if (coalescentOperation != null) {
                CoalescentDetails details = forward.coalescentDetails[coalescentOperation.outputBuffer];
                double[] zAdjoint = partialAdjoint[coalescentOperation.outputBuffer];
                double adjJ = 1.0 / details.j - dot(zAdjoint, details.w) / (details.j * details.j);

                for (int i = 0; i < stateCount; ++i) {
                    zBar[i] = zAdjoint[i];
                    adjW[i] = zBar[i] / details.j + adjJ;
                    leftEndBar[i] = partialAdjoint[coalescentOperation.accBuffer1][i] + adjW[i] * details.rightEnd[i] / sizes[i];
                    rightEndBar[i] = partialAdjoint[coalescentOperation.accBuffer2][i] + adjW[i] * details.leftEnd[i] / sizes[i];
                    populationSizeGradient[i] -= adjW[i] * details.leftEnd[i] * details.rightEnd[i] / (sizes[i] * sizes[i]);
                }

                addMatTVecAndOuterInPlace(partialAdjoint[coalescentOperation.inputBuffer1],
                        matrixAdjoint[coalescentOperation.inputMatrix1],
                        forward.matrices[coalescentOperation.inputMatrix1], leftEndBar, details.leftStart);
                addMatTVecAndOuterInPlace(partialAdjoint[coalescentOperation.inputBuffer2],
                        matrixAdjoint[coalescentOperation.inputMatrix2],
                        forward.matrices[coalescentOperation.inputMatrix2], rightEndBar, details.rightStart);
            }

            for (int idx = end - 1; idx >= start; --idx) {
                BranchIntervalOperation operation = branchOperations.get(idx);
                if (operation.inputBuffer2 >= 0) {
                    continue;
                }
                double[] yBar = partialAdjoint[operation.accBuffer1];
                double[] x = forward.partials[operation.inputBuffer1];
                addMatTVecAndOuterInPlace(partialAdjoint[operation.inputBuffer1],
                        matrixAdjoint[operation.inputMatrix1],
                        forward.matrices[operation.inputMatrix1], yBar, x);
            }
        }

        for (TransitionMatrixOperation operation : matrixOperations) {
            accumulateSpectralExpmAdjointEigenBasis(decompositions[operation.decompositionBuffer],
                    forward.spectralKernels[operation.outputBuffer], matrixAdjoint[operation.outputBuffer],
                    workspace.eigenBasisGradient[operation.decompositionBuffer], workspace);
            workspace.decompositionUsed[operation.decompositionBuffer] = true;
        }

        for (int buffer = 0; buffer < decompositions.length; ++buffer) {
            if (!workspace.decompositionUsed[buffer] || decompositions[buffer] == null) {
                continue;
            }
            accumulateGeneratorGradientFromEigenBasis(decompositions[buffer], workspace.eigenBasisGradient[buffer], rateGradient, workspace);
        }

        return new ReverseResult(copyMatrix(rateGradient), populationSizeGradient.clone(), matrixAdjoint);
    }

    DebugState debugState(List<BranchIntervalOperation> branchOperations,
                          List<TransitionMatrixOperation> matrixOperations,
                          List<Integer> intervalStarts) {
        ensureSupported(matrixOperations);
        ForwardResult forward = getOrCreateForward(branchOperations, matrixOperations, intervalStarts);
        ReverseResult reverse = reverse(branchOperations, matrixOperations, intervalStarts, forward);
        return new DebugState(forward, reverse);
    }

    double[][] debugMatrixAdjoint(List<BranchIntervalOperation> branchOperations,
                                  List<TransitionMatrixOperation> matrixOperations,
                                  List<Integer> intervalStarts,
                                  int matrixId) {
        return debugState(branchOperations, matrixOperations, intervalStarts).reverse.matrixAdjoint[matrixId];
    }

    double[][] debugForwardMatrix(List<BranchIntervalOperation> branchOperations,
                                  List<TransitionMatrixOperation> matrixOperations,
                                  List<Integer> intervalStarts,
                                  int matrixId) {
        return debugState(branchOperations, matrixOperations, intervalStarts).forward.matrices[matrixId];
    }

    double debugLikelihoodWithMatrixOverride(List<BranchIntervalOperation> branchOperations,
                                             List<TransitionMatrixOperation> matrixOperations,
                                             List<Integer> intervalStarts,
                                             int matrixId,
                                             double[][] replacementMatrix) {
        ensureSupported(matrixOperations);

        int bufferCount = maxBufferIndex(branchOperations) + 1;
        int matrixCount = maxMatrixIndex(matrixOperations) + 1;
        int intervalCount = intervalStarts.size() - 1;
        int maxIntervalNumber = maxIntervalNumber(branchOperations);

        double[][] partials = new double[bufferCount][];
        for (Map.Entry<Integer, double[]> entry : tipPartials.entrySet()) {
            partials[entry.getKey()] = entry.getValue().clone();
        }

        double[][][] matrices = new double[matrixCount][][];
        for (TransitionMatrixOperation operation : matrixOperations) {
            if (operation.outputBuffer == matrixId) {
                matrices[operation.outputBuffer] = replacementMatrix;
            } else {
                matrices[operation.outputBuffer] = computeTransitionProbabilities(operation.time, decompositions[operation.decompositionBuffer]);
            }
        }

        double[] coalescentProb = new double[maxIntervalNumber + 1];
        boolean[] hasCoalescentProb = new boolean[maxIntervalNumber + 1];
        double logLikelihood = 0.0;

        for (BranchIntervalOperation operation : branchOperations) {
            double[] left = matVec(matrices[operation.inputMatrix1], partials[operation.inputBuffer1]);
            if (operation.inputBuffer2 < 0) {
                partials[operation.outputBuffer] = left.clone();
                partials[operation.accBuffer1] = left.clone();
                continue;
            }

            double[] right = matVec(matrices[operation.inputMatrix2], partials[operation.inputBuffer2]);
            partials[operation.accBuffer1] = left.clone();
            partials[operation.accBuffer2] = right.clone();

            double[] normalized = new double[stateCount];
            double probability = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                normalized[i] = left[i] * right[i] / sizes[i];
                probability += normalized[i];
            }
            for (int i = 0; i < stateCount; ++i) {
                normalized[i] /= probability;
            }
            partials[operation.outputBuffer] = normalized;
            coalescentProb[operation.intervalNumber] = probability;
            hasCoalescentProb[operation.intervalNumber] = true;
        }

        for (int interval = 0; interval < intervalCount; ++interval) {
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);
            double[] e = new double[stateCount];
            double[] f = new double[stateCount];
            double[] g = new double[stateCount];
            double[] h = new double[stateCount];

            for (int idx = start; idx < end; ++idx) {
                BranchIntervalOperation operation = branchOperations.get(idx);
                accumulateSummary(partials[operation.inputBuffer1], e, f);
                if (operation.inputBuffer2 >= 0) {
                    accumulateSummary(partials[operation.inputBuffer2], e, f);
                }
                accumulateSummary(partials[operation.accBuffer1], g, h);
                if (operation.accBuffer2 >= 0) {
                    accumulateSummary(partials[operation.accBuffer2], g, h);
                }
            }

            double hazardSum = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                hazardSum += (e[i] * e[i] - f[i] + g[i] * g[i] - h[i]) / sizes[i];
            }

            BranchIntervalOperation anchor = branchOperations.get(start);
            logLikelihood += -anchor.intervalLength * hazardSum / 4.0;
            if (hasCoalescentProb[anchor.intervalNumber]) {
                logLikelihood += Math.log(coalescentProb[anchor.intervalNumber]);
            }
        }

        return logLikelihood;
    }

    private void pushReductionAdjoints(double[][] partials,
                                       double[][] partialAdjoint,
                                       int buffer,
                                       double[] adjSum,
                                       double[] adjSq) {
        if (buffer < 0) {
            return;
        }
        double[] vec = partials[buffer];
        double[] adj = partialAdjoint[buffer];
        for (int i = 0; i < stateCount; ++i) {
            adj[i] += adjSum[i] + 2.0 * vec[i] * adjSq[i];
        }
    }

    private BranchIntervalOperation firstCoalescentInInterval(List<BranchIntervalOperation> branchOperations, int start, int end) {
        for (int idx = start; idx < end; ++idx) {
            BranchIntervalOperation operation = branchOperations.get(idx);
            if (operation.inputBuffer2 >= 0) {
                return operation;
            }
        }
        return null;
    }

    private void computeTransitionProbabilitiesInPlace(double distance,
                                                       EigenDecomposition eigen,
                                                       double[][] target,
                                                       double[] iexp) {
        boolean real = eigen.getEigenValues().length == stateCount;
        double[] evec = eigen.getEigenVectors();
        double[] eval = eigen.getEigenValues();
        double[] ievc = eigen.getInverseEigenVectors();

        for (int i = 0; i < stateCount; i++) {
            if (real || eval[stateCount + i] == 0.0) {
                double exp = Math.exp(distance * eval[i]);
                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = ievc[i * stateCount + j] * exp;
                }
            } else {
                int i2 = i + 1;
                double b = eval[stateCount + i];
                double expat = Math.exp(distance * eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);
                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = expatcosbt * ievc[i * stateCount + j] +
                            expatsinbt * ievc[i2 * stateCount + j];
                    iexp[i2 * stateCount + j] = expatcosbt * ievc[i2 * stateCount + j] -
                            expatsinbt * ievc[i * stateCount + j];
                }
                i++;
            }
        }

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    sum += evec[i * stateCount + k] * iexp[k * stateCount + j];
                }
                target[i][j] = sum;
            }
        }
    }

    private void accumulateSpectralExpmAdjointEigenBasis(EigenDecomposition decomposition,
                                                         double[][] spectralKernel,
                                                         double[][] matrixAdjoint,
                                                         double[][] eigenBasisGradient,
                                                         ReverseWorkspace workspace) {
        double[] eigenVectors = decomposition.getEigenVectors();
        double[] inverseEigenVectors = decomposition.getInverseEigenVectors();
        double[] transformed = workspace.transformedFlat;
        double[] temp = workspace.projectedFlat;

        clearFlat(transformed);
        clearFlat(temp);

        for (int i = 0; i < stateCount; ++i) {
            double[] matrixAdjointRow = matrixAdjoint[i];
            int tempOffset = i * stateCount;
            for (int j = 0; j < stateCount; ++j) {
                double value = matrixAdjointRow[j];
                if (value == 0.0) {
                    continue;
                }
                for (int b = 0; b < stateCount; ++b) {
                    temp[tempOffset + b] += value * inverseEigenVectors[b * stateCount + j];
                }
            }
        }

        for (int i = 0; i < stateCount; ++i) {
            int tempOffset = i * stateCount;
            int eigRowOffset = i * stateCount;
            for (int a = 0; a < stateCount; ++a) {
                double viaLeft = eigenVectors[eigRowOffset + a];
                if (viaLeft == 0.0) {
                    continue;
                }
                int transformedOffset = a * stateCount;
                for (int b = 0; b < stateCount; ++b) {
                    transformed[transformedOffset + b] += viaLeft * temp[tempOffset + b];
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
            int offset = a * stateCount;
            for (int b = 0; b < stateCount; ++b) {
                eigenBasisGradient[a][b] += transformed[offset + b] * spectralKernel[a][b];
            }
        }
    }

    private void computeLoewnerKernelInPlace(double distance, EigenDecomposition decomposition, double[][] target) {
        double[] eigenValues = decomposition.getEigenValues();
        for (int a = 0; a < stateCount; ++a) {
            double lambdaA = distance * eigenValues[a];
            double expA = Math.exp(lambdaA);
            for (int b = 0; b < stateCount; ++b) {
                double lambdaB = distance * eigenValues[b];
                double kernel = Math.abs(lambdaA - lambdaB) < EIGEN_TOLERANCE
                        ? expA
                        : (expA - Math.exp(lambdaB)) / (lambdaA - lambdaB);
                target[a][b] = distance * kernel;
            }
        }
    }

    private void accumulateGeneratorGradientFromEigenBasis(EigenDecomposition decomposition,
                                                           double[][] eigenBasisGradient,
                                                           double[][] rateGradient,
                                                           ReverseWorkspace workspace) {
        double[] eigenVectors = decomposition.getEigenVectors();
        double[] inverseEigenVectors = decomposition.getInverseEigenVectors();
        double[] expmAdjoint = workspace.expmAdjointFlat;
        double[] temp = workspace.projectedFlat;

        clearFlat(expmAdjoint);
        clearFlat(temp);

        for (int i = 0; i < stateCount; ++i) {
            int tempOffset = i * stateCount;
            for (int a = 0; a < stateCount; ++a) {
                double value = inverseEigenVectors[a * stateCount + i];
                if (value == 0.0) {
                    continue;
                }
                for (int b = 0; b < stateCount; ++b) {
                    temp[tempOffset + b] += value * eigenBasisGradient[a][b];
                }
            }
        }

        for (int i = 0; i < stateCount; ++i) {
            int tempOffset = i * stateCount;
            int outOffset = i * stateCount;
            for (int b = 0; b < stateCount; ++b) {
                double value = temp[tempOffset + b];
                if (value == 0.0) {
                    continue;
                }
                for (int j = 0; j < stateCount; ++j) {
                    expmAdjoint[outOffset + j] += value * eigenVectors[j * stateCount + b];
                }
            }
        }

        for (int i = 0; i < stateCount; ++i) {
            int offset = i * stateCount;
            for (int j = 0; j < stateCount; ++j) {
                rateGradient[i][j] += expmAdjoint[offset + j];
            }
        }
    }

    private void matVecInPlace(double[][] matrix, double[] vector, double[] out) {
        for (int i = 0; i < matrix.length; ++i) {
            double sum = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                sum += matrix[i][j] * vector[j];
            }
            out[i] = sum;
        }
    }

    private double[][] computeTransitionProbabilities(double distance, EigenDecomposition eigen) {
        double[][] matrix = new double[stateCount][stateCount];
        double[] iexp = new double[stateCount * stateCount];
        computeTransitionProbabilitiesInPlace(distance, eigen, matrix, iexp);
        return matrix;
    }

    private double[] matVec(double[][] matrix, double[] vector) {
        double[] out = new double[matrix.length];
        matVecInPlace(matrix, vector, out);
        return out;
    }

    private void addInPlace(double[] target, double[] source) {
        for (int i = 0; i < target.length; ++i) {
            target[i] += source[i];
        }
    }

    private void addInPlace(double[][] target, double[][] source) {
        for (int i = 0; i < target.length; ++i) {
            for (int j = 0; j < target[i].length; ++j) {
                target[i][j] += source[i][j];
            }
        }
    }

    private void addOuterInPlace(double[][] target, double[] left, double[] right) {
        for (int i = 0; i < left.length; ++i) {
            double li = left[i];
            if (li == 0.0) {
                continue;
            }
            for (int j = 0; j < right.length; ++j) {
                target[i][j] += li * right[j];
            }
        }
    }

    private void addMatTVecInPlace(double[] target, double[][] matrix, double[] vector) {
        for (int i = 0; i < matrix.length; ++i) {
            double value = vector[i];
            if (value == 0.0) {
                continue;
            }
            double[] row = matrix[i];
            for (int j = 0; j < row.length; ++j) {
                target[j] += row[j] * value;
            }
        }
    }

    private void addMatTVecAndOuterInPlace(double[] transposeVecTarget,
                                           double[][] outerTarget,
                                           double[][] matrix,
                                           double[] left,
                                           double[] right) {
        for (int i = 0; i < matrix.length; ++i) {
            double leftValue = left[i];
            if (leftValue == 0.0) {
                continue;
            }
            double[] matrixRow = matrix[i];
            double[] outerRow = outerTarget[i];
            for (int j = 0; j < matrixRow.length; ++j) {
                transposeVecTarget[j] += matrixRow[j] * leftValue;
                outerRow[j] += leftValue * right[j];
            }
        }
    }

    private double dot(double[] left, double[] right) {
        double sum = 0.0;
        for (int i = 0; i < left.length; ++i) {
            sum += left[i] * right[i];
        }
        return sum;
    }

    private void accumulateSummary(double[] vec, double[] sum, double[] sqSum) {
        for (int i = 0; i < vec.length; ++i) {
            sum[i] += vec[i];
            sqSum[i] += vec[i] * vec[i];
        }
    }

    private void copyInto(double[] target, double[] source) {
        System.arraycopy(source, 0, target, 0, stateCount);
    }

    private void clearSquare(double[][] matrix) {
        for (int i = 0; i < matrix.length; ++i) {
            Arrays.fill(matrix[i], 0.0);
        }
    }

    private void clearFlat(double[] vector) {
        Arrays.fill(vector, 0.0);
    }

    private double[][] copyMatrix(double[][] matrix) {
        double[][] copy = new double[matrix.length][matrix[0].length];
        for (int i = 0; i < matrix.length; ++i) {
            System.arraycopy(matrix[i], 0, copy[i], 0, matrix[i].length);
        }
        return copy;
    }

    private int maxBufferIndex(List<BranchIntervalOperation> branchOperations) {
        int max = -1;
        for (Map.Entry<Integer, double[]> entry : tipPartials.entrySet()) {
            max = Math.max(max, entry.getKey());
        }
        for (BranchIntervalOperation operation : branchOperations) {
            max = Math.max(max, operation.outputBuffer);
            max = Math.max(max, operation.inputBuffer1);
            max = Math.max(max, operation.inputBuffer2);
            max = Math.max(max, operation.accBuffer1);
            max = Math.max(max, operation.accBuffer2);
        }
        return max;
    }

    private int maxMatrixIndex(List<TransitionMatrixOperation> matrixOperations) {
        int max = -1;
        for (TransitionMatrixOperation operation : matrixOperations) {
            max = Math.max(max, operation.outputBuffer);
            max = Math.max(max, operation.decompositionBuffer);
        }
        return max;
    }

    private int maxIntervalNumber(List<BranchIntervalOperation> branchOperations) {
        int max = -1;
        for (BranchIntervalOperation operation : branchOperations) {
            max = Math.max(max, operation.intervalNumber);
        }
        return max;
    }

    private ForwardWorkspace ensureForwardWorkspace(int bufferCount, int matrixCount, int intervalCount, int maxIntervalNumber) {
        if (forwardWorkspace == null ||
                forwardWorkspace.partials.length < bufferCount ||
                forwardWorkspace.matrices.length < matrixCount ||
                forwardWorkspace.intervalE.length < intervalCount ||
                forwardWorkspace.coalescentProb.length <= maxIntervalNumber) {
            forwardWorkspace = new ForwardWorkspace(bufferCount, matrixCount, intervalCount, maxIntervalNumber);
        }
        return forwardWorkspace;
    }

    private ReverseWorkspace ensureReverseWorkspace(int bufferCount, int matrixCount) {
        if (reverseWorkspace == null ||
                reverseWorkspace.partialAdjoint.length < bufferCount ||
                reverseWorkspace.matrixAdjoint.length < matrixCount) {
            reverseWorkspace = new ReverseWorkspace(bufferCount, matrixCount, stateCount);
        }
        return reverseWorkspace;
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate overrides calculateLikelihood directly");
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate overrides calculateLikelihood directly");
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate overrides calculateLikelihood directly");
    }

    @Override
    protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts,
                                                       List<TransitionMatrixOperation> matrixOperations,
                                                       List<BranchIntervalOperation> branchIntervalOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate uses reverse mode directly");
    }

    @Override
    protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate uses reverse mode directly");
    }

    @Override
    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts,
                                                                List<BranchIntervalOperation> branchIntervalOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate uses reverse mode directly");
    }

    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts,
                                                                     List<BranchIntervalOperation> branchIntervalOperations) {
        throw new UnsupportedOperationException("AdjointGenericBastaLikelihoodDelegate uses reverse mode directly");
    }

    @Override
    String getStamp() {
        return "adjoint-generic";
    }

    @Override
    public void setPartials(int index, double[] partials) {
        tipPartials.put(index, partials.clone());
        invalidateForwardCache();
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {
        if (index != 0) {
            throw new IllegalArgumentException("Adjoint route currently supports only eigen decomposition index 0");
        }
        if (transpose) {
            decomposition = decomposition.transpose();
        }
        decompositions[index] = decomposition;
        invalidateForwardCache();
    }

    @Override
    public void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        if (index != 0) {
            throw new IllegalArgumentException("Adjoint route currently supports only population-size buffer index 0");
        }
        if (sizes.length != stateCount) {
            throw new IllegalArgumentException("Expected exactly " + stateCount + " population sizes, found " + sizes.length);
        }
        System.arraycopy(sizes, 0, this.sizes, index * stateCount, stateCount);
        invalidateForwardCache();
    }

    private static final class ForwardResult {
        final double logLikelihood;
        final double[][] partials;
        final double[][][] matrices;
        final double[][] intervalE;
        final double[][] intervalF;
        final double[][] intervalG;
        final double[][] intervalH;
        final double[] coalescentProb;
        final boolean[] hasCoalescentProb;
        final CoalescentDetails[] coalescentDetails;
        final double[][][] spectralKernels;

        private ForwardResult(double logLikelihood,
                              double[][] partials,
                              double[][][] matrices,
                              double[][] intervalE,
                              double[][] intervalF,
                              double[][] intervalG,
                              double[][] intervalH,
                              double[] coalescentProb,
                              boolean[] hasCoalescentProb,
                              CoalescentDetails[] coalescentDetails,
                              double[][][] spectralKernels) {
            this.logLikelihood = logLikelihood;
            this.partials = partials;
            this.matrices = matrices;
            this.intervalE = intervalE;
            this.intervalF = intervalF;
            this.intervalG = intervalG;
            this.intervalH = intervalH;
            this.coalescentProb = coalescentProb;
            this.hasCoalescentProb = hasCoalescentProb;
            this.coalescentDetails = coalescentDetails;
            this.spectralKernels = spectralKernels;
        }
    }

    private static final class CoalescentDetails {
        double[] leftStart;
        double[] rightStart;
        double[] leftEnd;
        double[] rightEnd;
        final double[] w;
        double j;

        private CoalescentDetails(double[] w) {
            this.w = w;
        }
    }

    private final class ForwardWorkspace {
        final double[][] partials;
        final double[][][] matrices;
        final double[][][] spectralKernels;
        final double[][] intervalE;
        final double[][] intervalF;
        final double[][] intervalG;
        final double[][] intervalH;
        final double[] coalescentProb;
        final boolean[] hasCoalescentProb;
        final CoalescentDetails[] coalescentDetails;
        final double[] transitionScratch;

        private ForwardWorkspace(int bufferCount, int matrixCount, int intervalCount, int maxIntervalNumber) {
            this.partials = new double[bufferCount][stateCount];
            this.matrices = new double[matrixCount][stateCount][stateCount];
            this.spectralKernels = new double[matrixCount][stateCount][stateCount];
            this.intervalE = new double[intervalCount][stateCount];
            this.intervalF = new double[intervalCount][stateCount];
            this.intervalG = new double[intervalCount][stateCount];
            this.intervalH = new double[intervalCount][stateCount];
            this.coalescentProb = new double[maxIntervalNumber + 1];
            this.hasCoalescentProb = new boolean[maxIntervalNumber + 1];
            this.coalescentDetails = new CoalescentDetails[bufferCount];
            this.transitionScratch = new double[stateCount * stateCount];
            for (int i = 0; i < bufferCount; ++i) {
                coalescentDetails[i] = new CoalescentDetails(new double[stateCount]);
            }
        }

        private void reset(int intervalCount, int maxIntervalNumber) {
            Arrays.fill(coalescentProb, 0.0);
            Arrays.fill(hasCoalescentProb, 0, maxIntervalNumber + 1, false);
            for (int interval = 0; interval < intervalCount; ++interval) {
                Arrays.fill(intervalE[interval], 0.0);
                Arrays.fill(intervalF[interval], 0.0);
                Arrays.fill(intervalG[interval], 0.0);
                Arrays.fill(intervalH[interval], 0.0);
            }
        }
    }

    private final class ReverseWorkspace {
        final double[][] partialAdjoint;
        final double[][][] matrixAdjoint;
        final double[][][] eigenBasisGradient;
        final boolean[] decompositionUsed;
        final double[][] rateGradient;
        final double[] populationSizeGradient;
        final double[] adjE;
        final double[] adjF;
        final double[] adjG;
        final double[] adjH;
        final double[] zBar;
        final double[] adjW;
        final double[] leftEndBar;
        final double[] rightEndBar;
        final double[][] transformed;
        final double[][] projected;
        final double[][] expmAdjoint;
        final double[] transformedFlat;
        final double[] projectedFlat;
        final double[] expmAdjointFlat;

        private ReverseWorkspace(int bufferCount, int matrixCount, int stateCount) {
            this.partialAdjoint = new double[bufferCount][stateCount];
            this.matrixAdjoint = new double[matrixCount][stateCount][stateCount];
            this.eigenBasisGradient = new double[decompositions.length][stateCount][stateCount];
            this.decompositionUsed = new boolean[decompositions.length];
            this.rateGradient = new double[stateCount][stateCount];
            this.populationSizeGradient = new double[stateCount];
            this.adjE = new double[stateCount];
            this.adjF = new double[stateCount];
            this.adjG = new double[stateCount];
            this.adjH = new double[stateCount];
            this.zBar = new double[stateCount];
            this.adjW = new double[stateCount];
            this.leftEndBar = new double[stateCount];
            this.rightEndBar = new double[stateCount];
            this.transformed = new double[stateCount][stateCount];
            this.projected = new double[stateCount][stateCount];
            this.expmAdjoint = new double[stateCount][stateCount];
            this.transformedFlat = new double[stateCount * stateCount];
            this.projectedFlat = new double[stateCount * stateCount];
            this.expmAdjointFlat = new double[stateCount * stateCount];
        }

        private void reset(List<TransitionMatrixOperation> matrixOperations, double[][] partials) {
            for (int i = 0; i < partialAdjoint.length; ++i) {
                if (partials[i] != null) {
                    Arrays.fill(partialAdjoint[i], 0.0);
                }
            }
            for (TransitionMatrixOperation operation : matrixOperations) {
                clearSquare(matrixAdjoint[operation.outputBuffer]);
                clearSquare(eigenBasisGradient[operation.decompositionBuffer]);
                decompositionUsed[operation.decompositionBuffer] = false;
            }
            clearSquare(rateGradient);
            Arrays.fill(populationSizeGradient, 0.0);
        }
    }

    private static final class ReverseResult {
        final double[][] rateGradient;
        final double[] populationSizeGradient;
        final double[][][] matrixAdjoint;

        private ReverseResult(double[][] rateGradient, double[] populationSizeGradient, double[][][] matrixAdjoint) {
            this.rateGradient = rateGradient;
            this.populationSizeGradient = populationSizeGradient;
            this.matrixAdjoint = matrixAdjoint;
        }
    }

    static final class DebugState {
        final ForwardResult forward;
        final ReverseResult reverse;

        private DebugState(ForwardResult forward, ReverseResult reverse) {
            this.forward = forward;
            this.reverse = reverse;
        }
    }
}
