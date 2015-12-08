package dr.evomodelxml.operators;

import dr.evolution.tree.MultivariateTraitTree;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.operators.LatentFactorHamiltonianMC;
import dr.inference.model.CompoundParameter;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.CoercionMode;
import dr.xml.*;

/**
 * Created by max on 12/2/15.
 */
public class LatentFactorHamiltonianMCParser extends AbstractXMLObjectParser {
    public static final String LATENT_FACTOR_MODEL_HAMILTONIAN_MC="LatentFactorHamiltonianMC";
    public static final String WEIGHT="weight";
    public static final String N_STEPS="nSteps";
    public static final String STEP_SIZE="stepSize";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LatentFactorModel lfm=(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        AbstractMultivariateTraitLikelihood tree=(AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);
        double weight=xo.getDoubleAttribute(WEIGHT);
        CoercionMode mode=CoercionMode.parseMode(xo);
        int nSteps=xo.getIntegerAttribute(N_STEPS);
        double stepSize=xo.getDoubleAttribute(STEP_SIZE);


        return new LatentFactorHamiltonianMC(lfm, tree, weight, mode, stepSize, nSteps);


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newDoubleRule(STEP_SIZE),
            AttributeRule.newIntegerRule(N_STEPS),
            new ElementRule(LatentFactorModel.class),
            new ElementRule(AbstractMultivariateTraitLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return "Hamiltonian Monte Carlo For factors";
    }

    @Override
    public Class getReturnType() {
        return LatentFactorHamiltonianMC.class;
    }

    @Override
    public String getParserName() {
        return LATENT_FACTOR_MODEL_HAMILTONIAN_MC;
    }
}
