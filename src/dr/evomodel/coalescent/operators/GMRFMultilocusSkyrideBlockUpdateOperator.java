package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodelxml.coalescent.operators.GMRFSkyrideBlockUpdateOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import no.uib.cipr.matrix.*;

import java.util.logging.Logger;

/* A Metropolis-Hastings operator to update the log population sizes and precision parameter jointly under a Gaussian Markov random field prior
 *
 * @author Erik Bloomquist
 * @author Marc Suchard
 * @author Mandev Gill
 * @version $Id: GMRFMultilocusSkylineBlockUpdateOperator.java,v 1.5 2007/03/20 11:26:49 msuchard Exp $
 */
public class GMRFMultilocusSkyrideBlockUpdateOperator extends AbstractCoercableOperator {

    private double scaleFactor;
    private double lambdaScaleFactor;
    private int fieldLength;

    private int maxIterations;
    private double stopValue;

    private Parameter popSizeParameter;
    private Parameter precisionParameter;
    private Parameter lambdaParameter;

    GMRFMultilocusSkyrideLikelihood gmrfField;

    private double[] zeros;

    public GMRFMultilocusSkyrideBlockUpdateOperator(GMRFMultilocusSkyrideLikelihood gmrfLikelihood,
                                          double weight, CoercionMode mode, double scaleFactor,
                                          int maxIterations, double stopValue) {
        super(mode);
        gmrfField = gmrfLikelihood;
        popSizeParameter = gmrfLikelihood.getPopSizeParameter();
        precisionParameter = gmrfLikelihood.getPrecisionParameter();
        lambdaParameter = gmrfLikelihood.getLambdaParameter();

        this.scaleFactor = scaleFactor;
        lambdaScaleFactor = 0.0;
        fieldLength = popSizeParameter.getDimension();

        this.maxIterations = maxIterations;
        this.stopValue = stopValue;
        setWeight(weight);

        zeros = new double[fieldLength];
    }

    private double getNewLambda(double currentValue, double lambdaScale) {
        double a = MathUtils.nextDouble() * lambdaScale - lambdaScale / 2;
        double b = currentValue + a;
        if (b > 1)
            b = 2 - b;
        if (b < 0)
            b = -b;

        return b;
    }

    private double getNewPrecision(double currentValue, double scaleFactor) {
        double length = scaleFactor - 1 / scaleFactor;
        double returnValue;


        if (scaleFactor == 1)
            return currentValue;
        if (MathUtils.nextDouble() < length / (length + 2 * Math.log(scaleFactor))) {
            returnValue = (1 / scaleFactor + length * MathUtils.nextDouble()) * currentValue;
        } else {
            returnValue = Math.pow(scaleFactor, 2.0 * MathUtils.nextDouble() - 1) * currentValue;
        }

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


    public static double logGeneralizedDeterminant(UpperTriangBandMatrix Matrix) {
        double returnValue = 0;

        for (int i = 0; i < Matrix.numColumns(); i++) {
            if (Matrix.get(i, i) > 0.0000001) {
                returnValue += Math.log(Matrix.get(i, i));
            }
        }

        return returnValue;
    }
   
    public DenseVector newtonRaphson(double[] data1, double[] data2, DenseVector currentGamma, SymmTridiagMatrix proposedQ) throws OperatorFailedException {
        return newNewtonRaphson(data1, data2, currentGamma, proposedQ, maxIterations, stopValue);
    }

    public static DenseVector newNewtonRaphson(double[] data1, double[] data2, DenseVector currentGamma, SymmTridiagMatrix proposedQ,
                                               int maxIterations, double stopValue) throws OperatorFailedException {

        DenseVector iterateGamma = currentGamma.copy();
        DenseVector tempValue = currentGamma.copy();

        int numberIterations = 0;


        while (gradient(data1, data2, iterateGamma, proposedQ).norm(Vector.Norm.Two) > stopValue) {
           try {
                jacobian(data2, iterateGamma, proposedQ).solve(gradient(data1, data2, iterateGamma, proposedQ), tempValue);
           } catch (no.uib.cipr.matrix.MatrixNotSPDException e) {
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new OperatorFailedException("");
            } catch (no.uib.cipr.matrix.MatrixSingularException e) {
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new OperatorFailedException("");
            }     

            iterateGamma.add(tempValue);
            numberIterations++;

            if (numberIterations > maxIterations) {
                Logger.getLogger("dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator").fine("Newton-Raphson F");
                throw new OperatorFailedException("Newton Raphson algorithm did not converge within " + maxIterations + " step to a norm less than " + stopValue + "\n" +
                        "Try starting BEAST with a more accurate initial tree.");
            }
        }

        Logger.getLogger("dr.evomodel.coalescent.operators.GMRFMultilocusSkyrideBlockUpdateOperator").fine("Newton-Raphson S");
        return iterateGamma;

    }

    private static DenseVector gradient(double[] data1, double[] data2, DenseVector value, SymmTridiagMatrix Q) {

        DenseVector returnValue = new DenseVector(value.size());
        Q.mult(value, returnValue);
        for (int i = 0; i < value.size(); i++) {
            returnValue.set(i, -returnValue.get(i) - data1[i] + data2[i] * Math.exp(-value.get(i)));
        }
        return returnValue;
    }


    private static SPDTridiagMatrix jacobian(double[] data2, DenseVector value, SymmTridiagMatrix Q) {
        SPDTridiagMatrix jacobian = new SPDTridiagMatrix(Q, true);
        for (int i = 0, n = value.size(); i < n; i++) {
            jacobian.set(i, i, jacobian.get(i, i) + Math.exp(-value.get(i)) * data2[i]);
        }
        return jacobian;
    }

    public double doOperation() throws OperatorFailedException {

        double currentPrecision = precisionParameter.getParameterValue(0);
        double proposedPrecision = this.getNewPrecision(currentPrecision, scaleFactor);

        double currentLambda = this.lambdaParameter.getParameterValue(0);
        double proposedLambda = this.getNewLambda(currentLambda, lambdaScaleFactor);

        precisionParameter.setParameterValue(0, proposedPrecision);
        lambdaParameter.setParameterValue(0, proposedLambda);

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

        DenseVector modeForward = newtonRaphson(numCoalEv, wNative, currentGamma, proposedQ.copy());
       
        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeForward.get(i)));
            diagonal2.set(i, modeForward.get(i) + 1);

            forwardQW.set(i, i, diagonal1.get(i) + forwardQW.get(i, i));
            //diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - numCoalEv[i]);
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

        DenseVector modeBackward = newtonRaphson(numCoalEv, wNative, proposedGamma, currentQ.copy());

        for (int i = 0; i < fieldLength; i++) {
            diagonal1.set(i, wNative[i] * Math.exp(-modeBackward.get(i)));
            diagonal2.set(i, modeBackward.get(i) + 1);

            backwardQW.set(i, i, diagonal1.get(i) + backwardQW.get(i, i));
            //diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - 1);
            diagonal1.set(i, diagonal1.get(i) * diagonal2.get(i) - numCoalEv[i]);
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
        return GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
    }

    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1);
    }

    public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
        scaleFactor = 1 + value * value;
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
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

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);

        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


}

