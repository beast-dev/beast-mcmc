package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.BNPRLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import no.uib.cipr.matrix.BandCholesky;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;
import no.uib.cipr.matrix.UpperSPDBandMatrix;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRBlockUpdateOperator extends GMRFMultilocusSkyrideBlockUpdateOperator {

    private int fieldLength;

    private Parameter popSizeParameter;
    private Parameter precisionParameter;

    BNPRLikelihood bnprField;

    public static final double TWO_TIMES_PI =6.283185;

    public BNPRBlockUpdateOperator(BNPRLikelihood bnprLikelihood,
                                   double weight, CoercionMode mode, double scaleFactor,
                                   int maxIterations, double stopValue) {
        super(bnprLikelihood, weight, mode, scaleFactor, maxIterations, stopValue);

        this.bnprField = bnprLikelihood;
    }

    private DenseVector ESS(DenseVector currentGamma, double currentL) {
        // Choose ellipse
        DenseVector nu = new DenseVector(bnprField.getNumGridPoints()); // TODO: Implement!

        // Log-likelihood threshold
        double u = MathUtils.nextDouble();
        double logy = currentL + Math.log(u);

        // Draw a initial proposal, also defining a bracket
        double t = TWO_TIMES_PI + MathUtils.nextDouble();
        double tMin = t - TWO_TIMES_PI;
        double tMax = t;

        DenseVector q = currentGamma.copy();
        q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

        double l = bnprField.getLogLikelihood(); // TODO: How to get likelihood for almost-proposed points?

        while (l < logy) {

            if (t < 0) {
                tMin = t;
            } else {
                tMax = t;
            }

            t = tMin + MathUtils.nextDouble() * (tMax - tMin);
            q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

            l = bnprField.getLogLikelihood(); // TODO: How to get likelihood for almost-proposed points?
        }

        return q;
    }

    private double getNewPrecision(double currentValue) {
        //double length = scaleFactor - 1 / scaleFactor;
        double returnValue;
        double shape, scale;
        DenseVector currentGamma = new DenseVector(bnprField.getPopSizeParameter().getParameterValues());

        shape = bnprField.getPrecAlpha() + (bnprField.getNumGridPoints() - 1) / 2.0; // TODO: Watch for off-by-one
        scale = bnprField.getPrecBeta() + currentGamma.dot(bnprField.getStoredScaledWeightMatrix(currentValue).mult(currentGamma, new DenseVector(currentGamma.size())));

        returnValue = GammaDistribution.nextGamma(shape, scale);

        return returnValue;
    }

    public double doOperation() throws OperatorFailedException {
        double currentPrecision = precisionParameter.getParameterValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision);

        //double currentLambda = this.lambdaParameter.getParameterValue(0);
        //double proposedLambda = this.getNewLambda(currentLambda, lambdaScaleFactor);

        precisionParameter.setParameterValue(0, proposedPrecision);
        //lambdaParameter.setParameterValue(0, proposedLambda);

        DenseVector currentGamma = new DenseVector(gmrfField.getPopSizeParameter().getParameterValues());
        DenseVector proposedGamma;

        SymmTridiagMatrix currentQ = gmrfField.getStoredScaledWeightMatrix(currentPrecision, currentLambda);
        SymmTridiagMatrix proposedQ = gmrfField.getScaledWeightMatrix(proposedPrecision, proposedLambda);

        double[] wNative = gmrfField.getSufficientStatistics();
        double[] numCoalEv = gmrfField.getNumCoalEvents();

        UpperSPDBandMatrix forwardQW = new UpperSPDBandMatrix(proposedQ, 1);
        UpperSPDBandMatrix backwardQW = new UpperSPDBandMatrix(currentQ, 1);

        BandCholesky forwardCholesky = new BandCholesky(wNative.length, 1, true);
        BandCholesky backwardCholesky = new BandCholesky(wNative.length, 1, true);

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector diagonal2 = new DenseVector(fieldLength);
        DenseVector diagonal3 = new DenseVector(fieldLength);
        DenseVector ZBetaVector = getZBeta(covariates, betaParameter);
        DenseVector QZBetaProp = new DenseVector(fieldLength);
        DenseVector QZBetaCurrent = new DenseVector(fieldLength);
        forwardQW.mult(ZBetaVector, QZBetaProp);
        backwardQW.mult(ZBetaVector, QZBetaCurrent);

        DenseVector modeForward = newtonRaphson(numCoalEv, wNative, currentGamma, proposedQ.copy(), ZBetaVector);

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
            diagonal2.set(i, modeForward.get(i) + 1);

            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
            //diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
            diagonal1.set(i, QZBetaProp.get(i) + diagonal1.get(i) * diagonal2.get(i) - numCoalEv[i]);
        }

        forwardCholesky.factor(forwardQW.copy());

        DenseVector forwardMean = getMultiNormalMean(diagonal1, forwardCholesky);

        DenseVector stand_norm = new DenseVector(zeros);

        for (int i = 0; i < zeros.length; i++)
            stand_norm.set(i, MathUtils.nextGaussian());

        proposedGamma = getMultiNormal(stand_norm, forwardMean, forwardCholesky);

        /*
        double hRatio = 0;
        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, proposedGamma.get(i) - forwardMean.get(i));
        }
        diagonal3.zero();
        forwardQW.mult(diagonal1, diagonal3);

        hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * diagonal1.dot(diagonal3);
        */

        for (int i = 0; i < fieldLength; i++)
            popSizeParameter.setParameterValueQuietly(i, proposedGamma.get(i));

        ((Parameter.Abstract) popSizeParameter).fireParameterChangedEvent();


        double hRatio = 0;

        diagonal1.zero();
        diagonal2.zero();
        diagonal3.zero();

        DenseVector modeBackward = newtonRaphson(numCoalEv, wNative, proposedGamma, currentQ.copy(), ZBetaVector);

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
            diagonal2.set(i, modeBackward.get(i) + 1);

            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
            //diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
            diagonal1.set(i, QZBetaCurrent.get(i) + diagonal1.get(i) * diagonal2.get(i) - numCoalEv[i]);
        }

        backwardCholesky.factor(backwardQW.copy());

        DenseVector backwardMean = getMultiNormalMean(diagonal1, backwardCholesky);

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, currentGamma.get(i) - backwardMean.get(i));
        }

        backwardQW.mult(diagonal1, diagonal3);

        hRatio += logGeneralizedDeterminant(backwardCholesky.getU()) - 0.5 * diagonal1.dot(diagonal3);
        hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * stand_norm.dot(stand_norm);

        return hRatio;
    }
}
