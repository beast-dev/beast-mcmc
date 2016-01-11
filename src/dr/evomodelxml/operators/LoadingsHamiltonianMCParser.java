package dr.evomodelxml.operators;

import dr.evomodel.operators.LoadingsHamiltonianMC;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.CoercionMode;
import dr.xml.*;

/**
 * Created by max on 1/11/16.
 */
public class LoadingsHamiltonianMCParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_HAMILTONIAN_MC_PARSER="loadingsHamiltonianMCParser";
    public static final String WEIGHT="weight";
    public static final String STEP_SIZE="stepSize";
    public static final String N_STEPS="nSteps";
    public static final String MOMENTUM_SD="momentumSd";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        LatentFactorModel lfm=(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        MomentDistributionModel prior=(MomentDistributionModel) xo.getChild(MomentDistributionModel.class);
        double weight=xo.getDoubleAttribute(WEIGHT);
        CoercionMode mode=CoercionMode.parseMode(xo);
        int nSteps=xo.getIntegerAttribute(N_STEPS);
        double stepSize=xo.getDoubleAttribute(STEP_SIZE);
        double momentumSd= xo.getDoubleAttribute(MOMENTUM_SD);

        return new LoadingsHamiltonianMC(lfm, prior, weight, mode, stepSize, nSteps, momentumSd);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newIntegerRule(N_STEPS),
            AttributeRule.newDoubleRule(MOMENTUM_SD),
            new ElementRule(LatentFactorModel.class),
            new ElementRule(MomentDistributionModel.class),
    };


    @Override
    public String getParserDescription() {
        return "Hamiltonian Monte Carlo for loadings matrix in a latent factor model";
    }

    @Override
    public Class getReturnType() {
        return LoadingsHamiltonianMC.class;
    }

    @Override
    public String getParserName() {
        return LOADINGS_HAMILTONIAN_MC_PARSER;
    }
}
