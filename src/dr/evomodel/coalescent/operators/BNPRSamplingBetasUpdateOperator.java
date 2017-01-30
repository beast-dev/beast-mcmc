package dr.evomodel.coalescent.operators;

import dr.evomodel.coalescent.BNPRLikelihood;
import dr.evomodelxml.coalescent.operators.GMRFSkyrideBlockUpdateOperatorParser;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import no.uib.cipr.matrix.BandCholesky;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.SymmTridiagMatrix;
import no.uib.cipr.matrix.UpperSPDBandMatrix;

/**
 * Created by mkarcher on 9/30/16.
 */
public class BNPRSamplingBetasUpdateOperator extends AbstractCoercableOperator {

    private double scaleFactor;
    private int numBetas;

    private Parameter betas;

    BNPRLikelihood bnprField;

    public BNPRSamplingBetasUpdateOperator(BNPRLikelihood bnprLikelihood, double weight, CoercionMode mode,
                                           double scaleFactor) {
        super(mode);
        this.bnprField = bnprLikelihood;
        this.betas = bnprField.getSamplingBetas();
        this.numBetas = betas.getDimension();
    }

    private double getNewBeta(double currentBeta) {
        double newBeta = currentBeta + MathUtils.nextGaussian() * scaleFactor;

        return newBeta;
    }

    public double doOperation() throws OperatorFailedException {
        int betaIndex = MathUtils.nextInt(numBetas);
        double[] currentBetas = betas.getParameterValues();
        double[] newBetas = currentBetas.clone();

        newBetas[betaIndex] = getNewBeta(currentBetas[betaIndex]);

        betas.setParameterValue(betaIndex, newBetas[betaIndex]);

        double hRatio = bnprField.calculateLogSamplingLikelihoodSubBetas(newBetas) -
                bnprField.calculateLogSamplingLikelihoodSubBetas(currentBetas); // TODO: Check if need priors

        return hRatio;
    }


    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        //return GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
        return "bnprSamplingBetasUpdateOperator"; // TODO: Placeholder, redo this correctly
    }

    public double getCoercableParameter() {
//        return Math.log(scaleFactor);
        return Math.sqrt(scaleFactor - 1); // TODO: Ask about this transformation
    }

    public void setCoercableParameter(double value) {
//        scaleFactor = Math.exp(value);
        scaleFactor = 1 + value * value;
    }

    public double getRawParameter() {
        return scaleFactor;
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
