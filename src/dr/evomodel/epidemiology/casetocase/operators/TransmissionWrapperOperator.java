package dr.evomodel.epidemiology.casetocase.operators;

// for convergence, tree-modifying operators are better if the mean infectious period parameter is resampled
// immediately after the move. This class wraps another operator and performs this.

import dr.evomodel.epidemiology.casetocase.CaseToCaseTreeLikelihood;
import dr.evomodel.epidemiology.casetocase.WithinCaseCategoryOutbreak;
import dr.evomodel.epidemiology.casetocase.WithinCaseCoalescent;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;

public class TransmissionWrapperOperator extends AbstractCoercableOperator {

    private WithinCaseCoalescent c2cTL;
    private MCMCOperator operator;
    private final double variance;

    public static final String TRANSMISSION_WRAPPER_OPERATOR = "transmissionWrapperOperator";

    public TransmissionWrapperOperator(WithinCaseCoalescent c2cTL, MCMCOperator operator,
                                              double variance, CoercionMode mode){
        super(mode);
        this.c2cTL=c2cTL;
        this.operator = operator;
        this.variance = variance;
    }

    public TransmissionWrapperOperator(WithinCaseCoalescent c2cTL, MCMCOperator operator,
                                              double variance){
        this(c2cTL, operator, variance, CoercionMode.COERCION_OFF);
    }


    public String getPerformanceSuggestion() {
        return "Not implemented";
    }

    public String getOperatorName() {
        return TRANSMISSION_WRAPPER_OPERATOR+" ("+operator.getOperatorName()+")";
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

        c2cTL.debugOutputTree("before.nex");

        // todo this currently doesn't cope with categories

        double[] oldInfectiousPeriods = Arrays.copyOf(c2cTL.getInfectiousPeriods(), c2cTL.getOutbreak().size());
        double[] oldLatentPeriods = null;
        if(c2cTL.hasLatentPeriods()){
            oldLatentPeriods = Arrays.copyOf(c2cTL.getLatentPeriods(), c2cTL.getOutbreak().size());
        }

        double hr = operator.operate();

        c2cTL.debugOutputTree("after.nex");

        // todo this should come out eventually

        if(hr==Double.NaN){
            throw new RuntimeException("operator reporting NaN");
        }

        double[] newInfectiousPeriods = Arrays.copyOf(c2cTL.getInfectiousPeriods(), c2cTL.getOutbreak().size());
        double[] newLatentPeriods = null;
        if(c2cTL.hasLatentPeriods()){
            newLatentPeriods = Arrays.copyOf(c2cTL.getLatentPeriods(), c2cTL.getOutbreak().size());;
        }

        for(int i=0 ; i<oldInfectiousPeriods.length; i++){
            oldInfectiousPeriods[i] = Math.log(oldInfectiousPeriods[i]);
            newInfectiousPeriods[i] = Math.log(newInfectiousPeriods[i]);
            if(c2cTL.hasLatentPeriods()){
                oldLatentPeriods[i] = Math.log(oldLatentPeriods[i]);
                newLatentPeriods[i] = Math.log(newLatentPeriods[i]);
            }
        }

        double[] oldInfSS = c2cTL.getSummaryStatistics(oldInfectiousPeriods);
        double[] newInfSS = c2cTL.getSummaryStatistics(newInfectiousPeriods);

        LogNormalDistributionModel infDist
                = ((WithinCaseCategoryOutbreak)c2cTL.getOutbreak()).getInfectiousDist();
        double oldMean = infDist.getM();
        double oldStdev = infDist.getS();

        double proposalDistStdev = Math.sqrt(variance);

        double infNumerator = NormalDistribution.pdf(oldMean, oldInfSS[0], proposalDistStdev) *
                NormalDistribution.pdf(oldStdev, oldInfSS[3], proposalDistStdev);

        double newMean = MathUtils.nextGaussian()*proposalDistStdev + newInfSS[0];
        double newStdev = MathUtils.nextGaussian()*proposalDistStdev + newInfSS[3];

        if(newStdev<0){
            throw new OperatorFailedException("");
        }

        infDist.setM(newMean);
        infDist.setS(newStdev);

        double infDenominator = NormalDistribution.pdf(newMean, newInfSS[0], proposalDistStdev) *
                NormalDistribution.pdf(newStdev, newInfSS[3], proposalDistStdev);

        hr += Math.log(infNumerator) - Math.log(infDenominator);

        if(c2cTL.hasLatentPeriods()){
            double[] oldLatSS = c2cTL.getSummaryStatistics(oldLatentPeriods);
            double[] newLatSS = c2cTL.getSummaryStatistics(newLatentPeriods);

            LogNormalDistributionModel latDist
                    = ((WithinCaseCategoryOutbreak)c2cTL.getOutbreak()).getLatentDist();
            oldMean = latDist.getM();
            oldStdev = latDist.getS();

            double latNumerator = NormalDistribution.pdf(oldMean, oldLatSS[0], proposalDistStdev) *
                    NormalDistribution.pdf(oldStdev, oldLatSS[3], proposalDistStdev);

            newMean = MathUtils.nextGaussian()*proposalDistStdev + newLatSS[0];
            newStdev = MathUtils.nextGaussian()*proposalDistStdev + newLatSS[3];

            if(newStdev<0){
                throw new OperatorFailedException("");
            }

            latDist.setM(newMean);
            latDist.setS(newStdev);

            double latDenominator = NormalDistribution.pdf(newMean, newLatSS[0], proposalDistStdev) *
                    NormalDistribution.pdf(newStdev, newLatSS[3], proposalDistStdev);

            hr += Math.log(latNumerator) - Math.log(latDenominator);

        }

        return hr;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            WithinCaseCoalescent c2cL
                    = (WithinCaseCoalescent)xo.getChild(CaseToCaseTreeLikelihood.class);
            MCMCOperator treeOp = (MCMCOperator)xo.getChild(MCMCOperator.class);
            double variance = (double)xo.getDoubleAttribute("variance");

            CoercionMode mode = CoercionMode.COERCION_OFF;

            if(treeOp instanceof AbstractCoercableOperator){
                mode = CoercionMode.COERCION_ON;
            }

            return new TransmissionWrapperOperator(c2cL, treeOp, variance, mode);

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
            return TRANSMISSION_WRAPPER_OPERATOR;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                AttributeRule.newDoubleRule("variance"),
                new ElementRule(CaseToCaseTreeLikelihood.class, "The transmission tree likelihood element"),
                new ElementRule(MCMCOperator.class, "Another operator.")
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


