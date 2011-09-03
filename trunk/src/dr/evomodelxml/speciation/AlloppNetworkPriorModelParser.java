package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.AlloppNetworkPriorModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppNetworkPriorModelParser extends AbstractXMLObjectParser {

	public static final String ALLOPPNETWORKPRIORMODEL = "alloppNetworkPriorModel";
	public static final String EVENTRATE = "eventRate";
	
	public String getParserName() {
		return ALLOPPNETWORKPRIORMODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);
		final XMLObject erxo = xo.getChild(EVENTRATE);
		final Parameter eventrate = (Parameter) erxo.getChild(Parameter.class);
		return new AlloppNetworkPriorModel(eventrate, units);
	}

	private XMLSyntaxRule[] eventrateRules() {
				return new XMLSyntaxRule[]{
				new ElementRule(Parameter.class)
		};
	}
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				XMLUnits.SYNTAX_RULES[0],
				new ElementRule(EVENTRATE, eventrateRules())
		};
	}

	@Override
	public String getParserDescription() {
		return "Model for speciation, extinction, hybridization in allopolyploid network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppNetworkPriorModel.class;
	}

}
