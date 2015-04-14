package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessOperatorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_OPERATOR = "dpOperator";
	public static final String DATA_LOG_LIKELIHOOD = "dataLogLikelihood";
	public static final String MH_STEPS = "mhSteps";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		DirichletProcessPrior dpp = (DirichletProcessPrior) xo.getChild(DirichletProcessPrior.class);
//		CompoundLikelihood likelihood = (CompoundLikelihood) xo .getElementFirstChild(DATA_LOG_LIKELIHOOD);
		BeagleBranchLikelihood likelihood = (BeagleBranchLikelihood) xo .getElementFirstChild(DATA_LOG_LIKELIHOOD);
		Parameter categoriesParameter = (Parameter) xo.getElementFirstChild(  DirichletProcessPriorParser.CATEGORIES);
//		CountableRealizationsParameter countableRealizationsParameter = (CountableRealizationsParameter) xo.getChild(CountableRealizationsParameter.class);
		
		Parameter uniquelyRealizedParameters = (Parameter) xo.getChild(Parameter.class);
		
		int M = xo.getIntegerAttribute(MH_STEPS);
		final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

		return new DirichletProcessOperator(dpp, categoriesParameter, uniquelyRealizedParameters,
//				countableRealizationsParameter, 
				likelihood, M, weight);
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
