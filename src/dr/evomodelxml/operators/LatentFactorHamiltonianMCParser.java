package dr.evomodelxml.operators;

import dr.evolution.tree.MultivariateTraitTree;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.operators.LatentFactorHamiltonianMC;
import dr.inference.model.CompoundParameter;
import dr.inference.model.LatentFactorModel;
import dr.xml.*;

/**
 * Created by max on 12/2/15.
 */
public class LatentFactorHamiltonianMCParser extends AbstractXMLObjectParser {
    public static final String LATENT_FACTOR_MODEL_HAMILTONIAN_MC="LatentFactorHamiltonianMC";
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        LatentFactorModel lfm=(LatentFactorModel) xo.getChild(LatentFactorModel.class);
        AbstractMultivariateTraitLikelihood tree=(AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);

        return new LatentFactorHamiltonianMC(lfm, tree);


    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private static final XMLSyntaxRule[] rules = {
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
