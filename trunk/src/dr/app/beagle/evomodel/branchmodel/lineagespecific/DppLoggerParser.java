package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DppLoggerParser extends AbstractXMLObjectParser {

	public static final String DPP_LOGGER = "dppLogger";
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		Parameter categoryProbabilitiesParameter = (Parameter) xo.getElementFirstChild(DirichletProcessOperatorParser.CATEGORY_PROBABILITIES); 
		
		return new DppLogger(categoryProbabilitiesParameter, uniquelyRealizedParameters);
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserName() {
		return DPP_LOGGER;
	}

	@Override
	public String getParserDescription() {
		return DPP_LOGGER;
	}

	@Override
	public Class getReturnType() {
		return DppLogger.class;
	}

}// END: class
