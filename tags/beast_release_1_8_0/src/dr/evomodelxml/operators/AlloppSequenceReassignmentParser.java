package dr.evomodelxml.operators;


import dr.evomodel.operators.AlloppSequenceReassignment;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.inference.operators.MCMCOperator;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppSequenceReassignmentParser extends AbstractXMLObjectParser {
	public static final String SEQUENCE_REASSIGNMENT = "sequenceReassignment";
	
	
	public String getParserName() {
		return SEQUENCE_REASSIGNMENT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
		AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new AlloppSequenceReassignment(apspnet, apsp, weight);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class)
        };
	}

	@Override
	public String getParserDescription() {
		return "Operator which reassigns sequences within an allopolyploid species.";
	}

	@Override
	public Class getReturnType() {
		return AlloppSequenceReassignment.class;
	}

}
