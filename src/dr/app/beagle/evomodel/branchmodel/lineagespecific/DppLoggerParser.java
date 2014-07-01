package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.app.beagle.evomodel.branchmodel.lineagespecific.DppLogger.LogMode;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class DppLoggerParser extends AbstractXMLObjectParser {

	public static final String DPP_LOGGER = "dppLogger";
	public static final String MODE = "mode";
	
	public static final String NEW_VALUE = "newValue";
	public static final String NEW_CATEGORY = "newCategory";
	public static final String CATEGORY_PROBABILITIES = "categoryProbabilities";
	
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		CompoundParameter uniquelyRealizedParameters = (CompoundParameter)xo.getChild(CompoundParameter.class);
		Parameter categoryProbabilitiesParameter = (Parameter) xo.getElementFirstChild(DirichletProcessOperatorParser.CATEGORY_PROBABILITIES); 
		
		LogMode mode = LogMode.NEW_VALUE;
		String modeString = (String) xo.getAttribute(MODE);
	
		if(modeString.equalsIgnoreCase(NEW_VALUE)) {
			
			 mode = LogMode.NEW_VALUE;
			
		} else if(modeString.equalsIgnoreCase(NEW_CATEGORY)) {
			
			mode = LogMode.NEW_CATEGORY;
			
		} else if(modeString.equalsIgnoreCase(CATEGORY_PROBABILITIES)) {
			
			mode = LogMode.CATEGORY_PROBABILITIES;
			
		} else {
			
			//do nothing

		}// END: modeString check
		
		return new DppLogger(categoryProbabilitiesParameter, uniquelyRealizedParameters, mode);
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
