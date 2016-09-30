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

    private Parameter betas;

    BNPRLikelihood bnprField;

    public BNPRSamplingBetasUpdateOperator(BNPRLikelihood bnprLikelihood, double weight, CoercionMode mode,
                                           double scaleFactor) {
        super(mode);
        this.bnprField = bnprLikelihood;
    }

    public double doOperation() throws OperatorFailedException {

        
        double hRatio = 0;

        //hRatio += logGeneralizedDeterminant(backwardCholesky.getU()) - 0.5 * diagonal1.dot(diagonal3);
        //hRatio -= logGeneralizedDeterminant(forwardCholesky.getU() ) - 0.5 * stand_norm.dot(stand_norm);

        return hRatio;
    }


    //MCMCOperator INTERFACE

    public final String getOperatorName() {
        //return GMRFSkyrideBlockUpdateOperatorParser.BLOCK_UPDATE_OPERATOR;
        return "bnprBlockUpdateOperator"; // TODO: Placeholder, redo this correctly
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
