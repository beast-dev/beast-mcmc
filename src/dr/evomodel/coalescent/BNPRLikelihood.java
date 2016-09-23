package dr.evomodel.coalescent;

import dr.evolution.tree.Tree;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

import java.util.List;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRLikelihood extends GMRFMultilocusSkyrideLikelihood {

    protected double precAlpha;
    protected double precBeta;

    protected double[] samplingBetas;

    protected double logSamplingLikelihood;
    protected double storedLogSamplingLikelihood;

    protected double[] samplingTimes;
    protected boolean samplingTimesKnown;

    public BNPRLikelihood(List<Tree> treeList,
                          Parameter popParameter,
                          Parameter groupParameter,
                          Parameter precParameter,
                          Parameter lambda,
                          Parameter beta,
                          MatrixParameter dMatrix,
                          boolean timeAwareSmoothing,
                          double cutOff,
                          int numGridPoints,
                          Parameter phi,
                          Parameter ploidyFactorsParameter) {
        super(treeList, popParameter, groupParameter, precParameter, lambda, beta, dMatrix, timeAwareSmoothing,
                cutOff, numGridPoints, phi, ploidyFactorsParameter);
    }

    protected double[] getSamplingTimes() {
        if (!samplingTimesKnown) {
            Tree tree = this.getTree(0);
            int n = tree.getExternalNodeCount();
            double[] nodeHeights = new double[n];
            double maxHeight = 0;
            samplingTimes = new double[n];

            for (int i = 0; i < n; i++) {
                nodeHeights[i] = tree.getNodeHeight(tree.getExternalNode(i));
                if (nodeHeights[i] > maxHeight) {
                    maxHeight = nodeHeights[i];
                }
            }

            for (int i = 0; i < n; i++) {
                samplingTimes[i] = maxHeight - nodeHeights[i];
            }
        }

        return samplingTimes;
    }

    protected double calculateLogSamplingLikelihood() {

        if (!intervalsKnown) {
            // intervalsKnown -> false when handleModelChanged event occurs in super.
            wrapSetupIntervals();
            setupSufficientStatistics();
            intervalsKnown = true;
        }

        // Matrix operations taken from block update sampler to calculate data likelihood and field prior

        double currentLike = 0;
        double[] currentGamma = popSizeParameter.getParameterValues();

        for (int i = 0; i < fieldLength; i++) {
            currentLike += -numCoalEvents[i] * currentGamma[i] + ploidySums[i] - sufficientStatistics[i] * Math.exp(-currentGamma[i]);
        }

        return currentLike;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogCoalescentLikelihood();
            logFieldLikelihood = skygridHelper.getLogFieldLikelihood();
            logSamplingLikelihood = calculateLogSamplingLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood + logFieldLikelihood + logSamplingLikelihood;
    }

    public int getNumGridPoints() {
        return this.gridPoints.length;
    }

    public double getPrecAlpha() {
        return precAlpha;
    }

    public void setPrecAlpha(double precAlpha) {
        this.precAlpha = precAlpha;
    }

    public double getPrecBeta() {
        return precBeta;
    }

    public void setPrecBeta(double precBeta) {
        this.precBeta = precBeta;
    }

    public double[] getSamplingBetas() {
        return samplingBetas;
    }

    public void setSamplingBetas(double[] samplingBetas) {
        this.samplingBetas = samplingBetas;
    }
}
