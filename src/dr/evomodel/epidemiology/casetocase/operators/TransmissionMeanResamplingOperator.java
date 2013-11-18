package dr.evomodel.epidemiology.casetocase.operators;

// for convergence, tree-modifying operators are better if the mean infectious period parameter is resampled
// immediately after the move. This class wraps another operator and performs this.

import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.xml.*;

public class TransmissionMeanResamplingOperator extends AbstractCoercableOperator {

    private CaseToCaseTreeLikelihood c2cTL;
    private MCMCOperator operator;
    private Parameter variance;

    public static final String TRANSMISSION_MEAN_RESAMPLING_OPERATOR = "transmissionMeanResamplingOperator";

    public TransmissionMeanResamplingOperator(CaseToCaseTreeLikelihood c2cTL, MCMCOperator operator,
                                              Parameter variance, CoercionMode mode){
        super(mode);
        this.c2cTL=c2cTL;
        this.operator = operator;
        this.variance = variance;
    }

    public TransmissionMeanResamplingOperator(CaseToCaseTreeLikelihood c2cTL, MCMCOperator operator,
                                              Parameter variance){
        this(c2cTL, operator, variance, CoercionMode.COERCION_OFF);
    }


    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return TRANSMISSION_MEAN_RESAMPLING_OPERATOR+" ("+operator.getOperatorName()+")";
    }

    public void reset(){
        operator.reset();
        super.reset();
    }

    public void reject(){
        operator.reject();
        super.reject();
    }

    public void accept(double deviation){
        operator.accept(deviation);
        super.accept(deviation);
    }


    public double doOperation() throws OperatorFailedException {

        double numeratorMod = c2cTL.getMeanSamplingDistributionPDFValue(variance.getParameterValue(0));

        double hr = operator.operate();

        double denominatorMod = c2cTL.resampleInfPeriodDistMean(variance.getParameterValue(0));
        hr = hr + Math.log(numeratorMod) - Math.log(denominatorMod);

        return hr;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CaseToCaseTreeLikelihood c2cL
                    = (CaseToCaseTreeLikelihood)xo.getChild(CaseToCaseTreeLikelihood.class);
            MCMCOperator treeOp = (MCMCOperator)xo.getChild(MCMCOperator.class);
            Parameter variance = (Parameter)xo.getElementFirstChild("variance");

            CoercionMode mode = CoercionMode.COERCION_OFF;

            if(treeOp instanceof AbstractCoercableOperator){
                mode = CoercionMode.COERCION_ON;
            }

            return new TransmissionMeanResamplingOperator(c2cL,treeOp, variance, mode);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        public String getParserDescription() {
            return "Performs a tree move then resamples the parameter for mean of the case infectious periods.";
        }

        public Class getReturnType() {
            return TransmissionTreeOperator.class;
        }

        public String getParserName() {
            return TRANSMISSION_MEAN_RESAMPLING_OPERATOR;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(CaseToCaseTreeLikelihood.class, "The transmission tree likelihood element"),
                new ElementRule("variance", Parameter.class, "The variance of the proposal normal distribution"),
                new ElementRule(MCMCOperator.class, "A phylogenetic tree operator.")
        };
    };

    public double getCoercableParameter() {
        if(operator instanceof CoercableMCMCOperator){
            return ((CoercableMCMCOperator) operator).getCoercableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setCoercableParameter(double value) {
        if(operator instanceof CoercableMCMCOperator){
            ((CoercableMCMCOperator) operator).setCoercableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }

    public double getRawParameter() {
        if(operator instanceof CoercableMCMCOperator){
            return ((CoercableMCMCOperator) operator).getRawParameter();
        }
        throw new IllegalArgumentException();
    }
}


