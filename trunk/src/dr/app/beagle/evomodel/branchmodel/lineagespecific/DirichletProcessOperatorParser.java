package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.math.distributions.Distribution;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessOperatorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_OPERATOR = "dpOperator";
	public static final String DATA_LIKELIHOOD = "dataLogLike";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		DirichletProcessPrior dpp = (DirichletProcessPrior) xo.getChild(DirichletProcessPrior.class);

		
		CompoundLikelihood likelihod = (CompoundLikelihood) xo .getElementFirstChild(DATA_LIKELIHOOD);
//		Likelihood likelihod = (Likelihood) xo.getChild(Likelihood.class);
		
		Parameter zParameter = (Parameter) xo.getElementFirstChild(  DirichletProcessPriorParser.CATEGORIES);
		CountableRealizationsParameter countableRealizationsParameter = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);
		
		final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

		return new DirichletProcessOperator(dpp, zParameter, countableRealizationsParameter, likelihod, weight);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

				new ElementRule(DirichletProcessPrior.class, false),
		new ElementRule(Parameter.class, false), //
		AttributeRule.newDoubleRule(MCMCOperator.WEIGHT) //
		};

	}// END: getSyntaxRules

	@Override
	public String getParserName() {
		return DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public String getParserDescription() {
		return DIRICHLET_PROCESS_OPERATOR;
	}

	@Override
	public Class getReturnType() {
		return DirichletProcessOperator.class;
	}

}// END: class
