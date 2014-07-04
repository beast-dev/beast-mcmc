package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DirichletProcessOperatorParser extends AbstractXMLObjectParser {

	public static final String DIRICHLET_PROCESS_OPERATOR = "dpOperator";
	public static final String INTENSITY = "intensity";
	public static final String UNIQUE_REALIZATION_COUNT = "uniqueRealizationCount";
	public static final String CATEGORY_PROBABILITIES = "categoryProbabilities";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		Parameter zParameter = (Parameter) xo.getChild(Parameter.class);
		double intensity = xo.getDoubleAttribute(INTENSITY);
		int uniqueRealizationCount = xo
				.getIntegerAttribute(UNIQUE_REALIZATION_COUNT);

		Parameter categoryProbabilitiesParameter = null;
		if(xo.hasChildNamed(CATEGORY_PROBABILITIES)) {
			
			 categoryProbabilitiesParameter = (Parameter) xo.getElementFirstChild(CATEGORY_PROBABILITIES); 
		}
		
		return new DirichletProcessOperator(zParameter, categoryProbabilitiesParameter, intensity,
				uniqueRealizationCount);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {

				new ElementRule(CATEGORY_PROBABILITIES,
	                    new XMLSyntaxRule[] { new ElementRule(Parameter.class, false) }, true), // 
				
		new ElementRule(Parameter.class, false), //
				AttributeRule.newDoubleRule(INTENSITY), //
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
