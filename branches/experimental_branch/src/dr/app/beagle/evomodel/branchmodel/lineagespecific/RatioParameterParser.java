package dr.app.beagle.evomodel.branchmodel.lineagespecific;

import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class RatioParameterParser extends AbstractXMLObjectParser {

	public static final String RATIO_PARAMETER = "ratioParameter";
	public static final String PARAMETER1 = "parameter1";
	public static final String PARAMETER2 = "parameter2";
	
	@Override
	public String getParserName() {
		return RATIO_PARAMETER;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		 Parameter parameter1 = (Parameter) xo.getElementFirstChild(PARAMETER1);
		 Parameter parameter2 = (Parameter) xo.getElementFirstChild(PARAMETER2);
		
		 if(parameter1.getDimension() != parameter2.getDimension()){
			 throw new XMLParseException("Parameters in ratio '" + xo.getId() + "' must have the same dimension");
		 }
		
		
		return new RatioParameter(parameter1, parameter2);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		//TODO
		return null;
	}

	@Override
	public String getParserDescription() {
		return RATIO_PARAMETER;
	}

	@Override
	public Class getReturnType() {
		return Parameter.class;
	}

	
	
}//END: class
