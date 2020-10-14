package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ReversibleHMCProvider;
import dr.inference.model.Parameter;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.hmc.GeodesicHamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.HamiltonianMonteCarloOperator;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.util.Transform;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class GeodesicHamiltonianMonteCarloOperatorParser extends HamiltonianMonteCarloOperatorParser {
    public final static String OPERATOR_NAME = "geodesicHamiltonianMonteCarloOperator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return super.parseXMLObject(xo);
    }

    @Override
    protected HamiltonianMonteCarloOperator factory(AdaptationMode adaptationMode, double weight, GradientWrtParameterProvider derivative,
                                                    Parameter parameter, Transform transform, Parameter mask,
                                                    HamiltonianMonteCarloOperator.Options runtimeOptions, MassPreconditioner preconditioner,
                                                    ReversibleHMCProvider reversibleHMCprovider) {
        return new GeodesicHamiltonianMonteCarloOperator(adaptationMode, weight, derivative,
                parameter, transform, mask,
                runtimeOptions, preconditioner);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }


    @Override
    public String getParserDescription() {
        return "Returns a geodesic Hamiltonian Monte Carlo transition kernel";
    }

    @Override
    public Class getReturnType() {
        return GeodesicHamiltonianMonteCarloOperator.class;
    }

    @Override
    public String getParserName() {
        return OPERATOR_NAME;
    }
}
