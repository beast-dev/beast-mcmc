package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.BNPRLikelihood;
import dr.evomodelxml.coalescent.operators.BNPRBlockUpdateOperatorParser;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import no.uib.cipr.matrix.*;

import java.util.List;

/**
 * Created by mkarcher on 9/15/16.
 */
public class BNPRBlockUpdateOperator extends SimpleMCMCOperator implements GibbsOperator {

    private int fieldLength;

    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private List<Parameter> betaParameter;
    private List<MatrixParameter> covariates;

    BNPRLikelihood bnprField;

    private double[] zeros;

    public static final double TWO_TIMES_PI =6.283185;

    public BNPRBlockUpdateOperator(BNPRLikelihood bnprLikelihood, double weight) {
        bnprField = bnprLikelihood;
        popSizeParameter = bnprLikelihood.getPopSizeParameter();
        precisionParameter = bnprLikelihood.getPrecisionParameter();

        betaParameter = bnprLikelihood.getBetaListParameter();
        covariates = bnprLikelihood.getCovariates();

        fieldLength = popSizeParameter.getDimension();

        setWeight(weight);

        this.bnprField = bnprLikelihood;

        this.zeros = new double[fieldLength];
    }

    private DenseVector ESS(DenseVector currentGamma, double currentLoglik) {
        // Choose ellipse
        DenseVector nu = new DenseVector(bnprField.getNumGridPoints()); // TODO: Implement!

        // Log-likelihood threshold
        double u = MathUtils.nextDouble();
        double logy = currentLoglik + Math.log(u);

        // Draw a initial proposal, also defining a bracket
        double t = TWO_TIMES_PI + MathUtils.nextDouble();
        double tMin = t - TWO_TIMES_PI;
        double tMax = t;

        DenseVector q = currentGamma.copy();
        q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

        double l = bnprField.getLogLikelihoodSubGamma(q.getData());

        while (l < logy) {

            if (t < 0) {
                tMin = t;
            } else {
                tMax = t;
            }

            t = tMin + MathUtils.nextDouble() * (tMax - tMin);
            q = (DenseVector) q.scale(Math.cos(t)).add(Math.sin(t), nu);

            l = bnprField.getLogLikelihoodSubGamma(q.getData());
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

    public DenseVector getMultiNormalMean(DenseVector CanonVector, BandCholesky Cholesky) {

        DenseVector tempValue = new DenseVector(zeros);
        DenseVector Mean = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

        // Assume Cholesky factorization of the precision matrix Q = LL^T

        // 1. Solve L\omega = b

        CholeskyUpper.transSolve(CanonVector, tempValue);

        // 2. Solve L^T \mu = \omega

        CholeskyUpper.solve(tempValue, Mean);

        return Mean;
    }

    public DenseVector getMultiNormal(DenseVector StandNorm, DenseVector Mean, BandCholesky Cholesky) {

        DenseVector returnValue = new DenseVector(zeros);

        UpperTriangBandMatrix CholeskyUpper = Cholesky.getU();

        // 3. Solve L^T v = z

        CholeskyUpper.solve(StandNorm, returnValue);

        // 4. Return x = \mu + v

        returnValue.add(Mean);

        return returnValue;
    }


    public static DenseVector getMultiNormal(DenseVector Mean, UpperSPDDenseMatrix Variance) {
        int length = Mean.size();
        DenseVector tempValue = new DenseVector(length);
        DenseVector returnValue = new DenseVector(length);
        UpperSPDDenseMatrix ab = Variance.copy();

        for (int i = 0; i < returnValue.size(); i++)
            tempValue.set(i, MathUtils.nextGaussian());

        DenseCholesky chol = new DenseCholesky(length, true);
        chol.factor(ab);

        UpperTriangDenseMatrix x = chol.getU();

        x.transMult(tempValue, returnValue);
        returnValue.add(Mean);
        return returnValue;
    }

    public double doOperation() throws OperatorFailedException {
        double currentPrecision = precisionParameter.getParameterValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision);

        //double currentLambda = this.lambdaParameter.getParameterValue(0);
        //double proposedLambda = this.getNewLambda(currentLambda, lambdaScaleFactor);

        precisionParameter.setParameterValue(0, proposedPrecision);
        //lambdaParameter.setParameterValue(0, proposedLambda);

        DenseVector currentGamma = new DenseVector(bnprField.getPopSizeParameter().getParameterValues());
        DenseVector proposedGamma;

        SymmTridiagMatrix currentQ = bnprField.getStoredScaledWeightMatrix(currentPrecision, currentLambda);
        SymmTridiagMatrix proposedQ = bnprField.getScaledWeightMatrix(proposedPrecision, proposedLambda);

        double[] wNative = bnprField.getSufficientStatistics();
        double[] numCoalEv = bnprField.getNumCoalEvents();

        UpperSPDBandMatrix forwardQW = new UpperSPDBandMatrix(proposedQ, 1);
        UpperSPDBandMatrix backwardQW = new UpperSPDBandMatrix(currentQ, 1);

        BandCholesky forwardCholesky = new BandCholesky(wNative.length, 1, true);
        BandCholesky backwardCholesky = new BandCholesky(wNative.length, 1, true);

        DenseVector diagonal1 = new DenseVector(fieldLength);
        DenseVector diagonal2 = new DenseVector(fieldLength);
        DenseVector diagonal3 = new DenseVector(fieldLength);
        //DenseVector ZBetaVector = getZBeta(covariates, betaParameter);
        //DenseVector QZBetaProp = new DenseVector(fieldLength);
        //DenseVector QZBetaCurrent = new DenseVector(fieldLength);
        //forwardQW.mult(ZBetaVector, QZBetaProp);
        //backwardQW.mult(ZBetaVector, QZBetaCurrent);

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

    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        return BNPRBlockUpdateOperatorParser.BNPR_BLOCK_OPERATOR;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {
        return "This operator should not need tuning, and should accept with probability 1.";
    }
}
