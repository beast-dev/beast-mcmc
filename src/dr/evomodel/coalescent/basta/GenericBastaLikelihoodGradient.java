package dr.evomodel.coalescent.basta;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

import java.util.Arrays;
import java.util.List;

/**
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class GenericBastaLikelihoodGradient
//extends ProcessOnCoalescentIntervalDelegate
        implements ProcessOnCoalescentIntervalDelegate, GradientWrtParameterProvider {

    private final BastaLikelihood likelihood;
    private final GenericBastaLikelihoodDelegate likelihoodDelegate;
    private final CoalescentIntervalTraversal treeTraversalDelegate;
    private final GenericBastaLikelihoodDelegate.InternalStorage likelihoodStorage;
    private final SubstitutionModel substitutionModel; // TODO generify for multiple substitution models (e.g., epoch)

    private final int stateCount;
    private final Tree tree;

    // wrt parameter stuff // TODO delegate to separate classes
    private final Parameter popSizeParameter = null; // TODO

    // Internal storage // TODO greatly simplify
    private final double[][][] partialsGrad;
    private final double[][][] matricesGrad;
    private final double[][][] coalescentGrad;
    private final double[][][] eGrad;
    private final double[][][] fGrad;
    private final double[][][] gGrad;
    private final double[][][] hGrad;
    private final double[][] partialsGradPopSize;
    private final double[][] coalescentGradPopSize;
    private final double[][] eGradPopSize;
    private final double[][] fGradPopSize;
    private final double[][] gGradPopSize;
    private final double[][] hGradPopSize;


    public GenericBastaLikelihoodGradient(BastaLikelihood likelihood) {

        this.likelihood = likelihood;
        this.likelihoodDelegate = (GenericBastaLikelihoodDelegate) likelihood.getLikelihoodDelegate(); // TODO make type-safe
        this.treeTraversalDelegate = likelihood.getTraversalDelegate();

        this.stateCount = likelihoodDelegate.stateCount;
        int maxNumCoalescentIntervals = likelihoodDelegate.maxNumCoalescentIntervals;
        this.tree = likelihoodDelegate.tree;
        this.likelihoodStorage = likelihoodDelegate.getInternalStorage();
        this.substitutionModel = likelihood.getSubstitutionModel();

        this.partialsGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * (tree.getNodeCount() + 1) * stateCount];
        this.matricesGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount * stateCount];
        this.coalescentGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals];

        this.eGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.fGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.gGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
        this.hGrad = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];

        this.partialsGradPopSize = new double[stateCount][maxNumCoalescentIntervals * (tree.getNodeCount() + 1) * stateCount];
        this.coalescentGradPopSize = new double[stateCount][maxNumCoalescentIntervals];

        this.eGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.fGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.gGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        this.hGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
    }


    private void peelPartialsGrad(double[] partials,
                                  double distance, int resultOffset,
                                  int leftPartialOffset, int rightPartialOffset,
                                  double[] matrices,
                                  int leftMatrixOffset, int rightMatrixOffset,
                                  int leftAccOffset, int rightAccOffset,
                                  double[] probability, int probabilityOffset,
                                  double[] sizes, int sizesOffset,
                                  int stateCount) {

        resultOffset *= stateCount;

        // Handle left
        leftPartialOffset *= stateCount;
        leftMatrixOffset *= stateCount * stateCount;
        leftAccOffset *= stateCount;

        // TODO
        boolean transpose = false;
        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double sum = 0.0;
                    if (transpose && i == b) {
                        sum += partials[leftAccOffset + a] * distance;
                    }
                    if (!transpose) {
                        for (int j = 0; j < stateCount; ++j) {
                            sum += matricesGrad[a][b][leftMatrixOffset + i * stateCount + j] * partials[leftPartialOffset + j];
                        }
                    }
                    for (int j = 0; j < stateCount; ++j) {
                        sum += matrices[leftMatrixOffset + i * stateCount + j] * partialsGrad[a][b][leftPartialOffset + j];
                    }
                    partialsGrad[a][b][resultOffset + i] = sum;

                    throw new RuntimeException("Function should not depend on `transpose`");
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
            for (int i = 0; i < stateCount; ++i) {
                double sum = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    sum += matrices[leftMatrixOffset + i * stateCount + j] * partialsGradPopSize[a][leftPartialOffset + j];
                }
                partialsGradPopSize[a][resultOffset + i] = sum;
            }
        }

        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            rightAccOffset *= stateCount;
            // TODO: check bug?
            sizesOffset *= sizesOffset * stateCount;

            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    double J = probability[probabilityOffset];

                    // first half
                    double partial_J_ab = 0.0;
                    for (int i = 0; i < stateCount; ++i) {
                        double rightGrad = 0.0;
                        if (transpose && i == b) {
                            rightGrad += partials[rightAccOffset + a] * distance;
                        }

                        if (!transpose) {
                            for (int j = 0; j < stateCount; ++j) {
                                rightGrad += matricesGrad[a][b][rightMatrixOffset + i * stateCount + j] * partials[rightPartialOffset + j];
                            }
                        }
                        for (int j = 0; j < stateCount; ++j) {
                            rightGrad += matrices[rightMatrixOffset + i * stateCount + j] * partialsGrad[a][b][rightPartialOffset + j];
                        }
                        double leftGrad = partialsGrad[a][b][resultOffset + i];
                        double left = partials[leftAccOffset + i];
                        double right = partials[rightAccOffset + i];

                        double entry = (leftGrad * right + rightGrad * left) / sizes[sizesOffset + i];
                        partial_J_ab += entry;

                        partialsGrad[a][b][resultOffset + i] = entry / J;
                        partialsGrad[a][b][leftAccOffset + i] = leftGrad;
                        partialsGrad[a][b][rightAccOffset + i] = rightGrad;

                        throw new RuntimeException("Function should not depend on `transpose`");
                    }
                    // second half
                    for (int i = 0; i < stateCount; ++i) {
                        double entry = partials[resultOffset + i];
                        partialsGrad[a][b][resultOffset + i] -= partial_J_ab * entry / J;
                    }
                    coalescentGrad[a][b][probabilityOffset] = partial_J_ab;
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                double J = probability[probabilityOffset];
                // first half
                double partial_J_ab_PopSize = 0.0;
                for (int i = 0; i < stateCount; ++i) {
                    double rightGradPopSize = 0.0;
                    for (int j = 0; j < stateCount; ++j) {
                        rightGradPopSize += matrices[rightMatrixOffset + i * stateCount + j] *  partialsGradPopSize[a][rightPartialOffset + j];
                    }
                    double leftGradPopSize = partialsGradPopSize[a][resultOffset + i];
                    double left = partials[leftAccOffset + i];
                    double right = partials[rightAccOffset + i];

                    double entry = (leftGradPopSize * right + rightGradPopSize * left) / sizes[sizesOffset + i];
                    if (i == a){
                        entry += left * right;
                    }
                    partial_J_ab_PopSize += entry;

                    partialsGradPopSize[a][resultOffset + i] = entry / J;
                    partialsGradPopSize[a][leftAccOffset + i] = leftGradPopSize;
                    partialsGradPopSize[a][rightAccOffset + i] = rightGradPopSize;
                }
                // second half
                for (int i = 0; i < stateCount; ++i) {
                    double entry = partials[resultOffset + i];
                    partialsGradPopSize[a][resultOffset + i] -= partial_J_ab_PopSize * entry / J;
                }
                coalescentGradPopSize[a][probabilityOffset] = partial_J_ab_PopSize;
            }
        }
    }

    private void computeTransitionProbabilitiesGrad(double distance, int matrixOffset) {
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                for (int c = 0; c < stateCount; c++) {
                    for (int d = 0; d < stateCount; d++) { // TODO MAS: last loop unnecessary (also S^4 storage is unnecessary)
                        if (d == b) {
                            // TODO MAS: should these be cached at all? why not generate on the fly (t * matrices[])
                            matricesGrad[a][b][matrixOffset + c*stateCount + b] =  distance * likelihoodStorage.matrices[matrixOffset + c*stateCount + a];
                        } else {
                            matricesGrad[a][b][matrixOffset + c*stateCount + d] = 0; // TODO MAS: avoid caching (many) zeros
                        }
                    }
                }
            }
        }
    }

    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts,
                                                                List<BranchIntervalOperation> branchIntervalOperations) {

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinIntervalGrad(likelihoodStorage.partials, partialsGrad,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);

            }
        }

        double[][] grad = new double[stateCount][stateCount];
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            double[][] temp_grad =  reduceAcrossIntervalsGrad(
                    likelihoodStorage.e, likelihoodStorage.f,
                    likelihoodStorage.g, likelihoodStorage.h,
                    operation.intervalNumber, operation.intervalLength,
                    likelihoodStorage.sizes, likelihoodStorage.coalescent, stateCount);

            for (int a = 0; a < stateCount; a++) {
                for (int b = 0; b < stateCount; b++) {
                    grad[a][b] += temp_grad[a][b];
                }
            }
        }
        return grad;
    }

    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.stream(eGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(fGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(gGradPopSize).forEach(a -> Arrays.fill(a, 0));
        Arrays.stream(hGradPopSize).forEach(a -> Arrays.fill(a, 0));

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinIntervalGrad(likelihoodStorage.partials, partialsGrad,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);
            }
        }

        double[] grad = new double[stateCount];
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            double[] temp_grad =  reduceAcrossIntervalsGradPopSize(
                    likelihoodStorage.e, likelihoodStorage.f,
                    likelihoodStorage.g, likelihoodStorage.h,
                    operation.intervalNumber, operation.intervalLength,
                    likelihoodStorage.sizes, likelihoodStorage.coalescent, stateCount);

            for (int a = 0; a < stateCount; a++) {
                    grad[a] += temp_grad[a];
                }
            }
        return grad;
    }

    private double[][] reduceAcrossIntervalsGrad(double[] e, double[] f, double[] g, double[] h,
                                         int interval, double length,
                                         double[] sizes, double[] coalescent,
                                         int stateCount) {

        int offset = interval * stateCount;
        double[][] grad = new double[stateCount][stateCount];

        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * eGrad[a][b][offset + k] - fGrad[a][b][offset + k] +
                            2 * g[offset + k] * gGrad[a][b][offset + k] - hGrad[a][b][offset + k]) / sizes[k];
                }

                grad[a][b] = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    grad[a][b] += coalescentGrad[a][b][interval] / J;
                }
            }
        }
        return grad;
    }

    private double[] reduceAcrossIntervalsGradPopSize(double[] e, double[] f, double[] g, double[] h,
                                                 int interval, double length,
                                                 double[] sizes, double[] coalescent,
                                                 int stateCount) {

        int offset = interval * stateCount;
        double[] grad = new double[stateCount];

        for (int a = 0; a < stateCount; a++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * eGradPopSize[a][offset + k] - fGradPopSize[a][offset + k] +
                            2 * g[offset + k] * gGradPopSize[a][offset + k] - hGradPopSize[a][offset + k]) / sizes[k];
                    if (k == a) {
                        sum += (e[offset + k] * e[offset + k]) - f[offset + k] +
                                (g[offset + k] * g[offset + k]) - h[offset + k];
                    }
                }


                grad[a] = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    grad[a] += coalescentGradPopSize[a][interval] / J;
                }
            }
        return grad;
    }

    private void reduceWithinIntervalGrad(double[] partials, double[][][] partialsGrad,
                                          int startBuffer1, int startBuffer2,
                                          int endBuffer1, int endBuffer2,
                                          int interval, int stateCount)  {
        interval *= stateCount;

        startBuffer1 *= stateCount;
        endBuffer1 *= stateCount;

        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGrad = partialsGrad[a][b][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    eGrad[a][b][interval + i] += startPGrad;
                    fGrad[a][b][interval + i] += 2 * startP * startPGrad;

                    double endPGrad = partialsGrad[a][b][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gGrad[a][b][interval + i] += endPGrad;
                    hGrad[a][b][interval + i] += 2 * endP * endPGrad;
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGradPopSize = partialsGradPopSize[a][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    eGradPopSize[a][interval + i] += startPGradPopSize;
                    fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                    double endPGradPopSize = partialsGradPopSize[a][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gGradPopSize[a][interval + i] += endPGradPopSize;
                    hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                }
            }

        if (startBuffer2 >= 0) {
            startBuffer2 *= stateCount;
            endBuffer2 *= stateCount;
            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGrad = partialsGrad[a][b][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        eGrad[a][b][interval + i] += startPGrad;
                        fGrad[a][b][interval + i] += 2 * startP * startPGrad;

                        double endPGrad = partialsGrad[a][b][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gGrad[a][b][interval + i] += endPGrad;
                        hGrad[a][b][interval + i] += 2 * endP * endPGrad;
                    }
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGradPopSize = partialsGradPopSize[a][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        eGradPopSize[a][interval + i] += startPGradPopSize;
                        fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                        double endPGradPopSize = partialsGradPopSize[a][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gGradPopSize[a][interval + i] += endPGradPopSize;
                        hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                    }

            }
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return null; // TODO delegate depending on parameter
    }

    @Override
    public int getDimension() {
        return 0; // TODO delegate depending on parameter
    }

    @Override
    public double[] getGradientLogDensity() {
        return new double[0]; // TODO delegate depending on parameter
    }

    private double[] getGradientLogDensityWrtSubstitutionModel() {

        assert(substitutionModel instanceof SVSComplexSubstitutionModel);

        SVSComplexSubstitutionModel svsComplexSubstitutionModel = (SVSComplexSubstitutionModel) substitutionModel;
        Parameter parameters = svsComplexSubstitutionModel.getRatesParameter();


        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
//                transitionMatricesKnown ? NO_OPT :
                        treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        final NodeRef root = tree.getRoot();

        likelihood.getLogLikelihood();
        // log likelihood

        double[][] full_gradient = likelihoodDelegate.calculateGradient(branchOperations, matrixOperations, intervalStarts, root.getNumber());
        double[] gradient = new double[stateCount*(stateCount-1)];

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                gradient[k] = (full_gradient[i][j] - full_gradient[i][i]) * substitutionModel.getFrequencyModel().getFrequency(j) ;
                k += 1;
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                gradient[k] =(full_gradient[i][j] - full_gradient[i][i]) * substitutionModel.getFrequencyModel().getFrequency(j);
                k += 1;
            }
        }
        return gradient;
    }

    private double[] getGradientLogDensityWrtPopSizes() {

        Parameter parameters = popSizeParameter;


        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
//                transitionMatricesKnown ? NO_OPT :
                        treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        final NodeRef root = tree.getRoot();

        likelihood.getLogLikelihood();

        double[] full_gradient =  likelihoodDelegate.calculateGradientPopSize(branchOperations, matrixOperations, intervalStarts, root.getNumber());
        double[] gradient = new double[stateCount];
        for (int i = 0; i < stateCount; ++i) {
            gradient[i] = -full_gradient[i]*Math.pow(parameters.getParameterValue(i), -2);
        }
        return gradient;
    }
}

