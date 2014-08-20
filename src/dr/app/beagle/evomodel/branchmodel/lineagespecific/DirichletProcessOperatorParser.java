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
	public static final String UNIQUE_REALIZATION_COUNT = "uniqueRealizationCount";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		DirichletProcessPrior dpp = (DirichletProcessPrior) xo.getChild(DirichletProcessPrior.class);
		
		Parameter zParameter = (Parameter) xo.getChild(Parameter.class);
		int uniqueRealizationCount = xo
				.getIntegerAttribute(UNIQUE_REALIZATION_COUNT);

		 final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
		
		return new DirichletProcessOperator(dpp, zParameter,
				uniqueRealizationCount, weight);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

				new ElementRule(DirichletProcessPrior.class, false),
		new ElementRule(Parameter.class, false), //
		AttributeRule.newDoubleRule(MCMCOperator.WEIGHT), //
		AttributeRule.newIntegerRule(UNIQUE_REALIZATION_COUNT)

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
