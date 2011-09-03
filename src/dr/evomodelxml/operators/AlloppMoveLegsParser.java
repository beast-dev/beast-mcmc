package dr.evomodelxml.operators;

import dr.evomodel.operators.AlloppMoveLegs;
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
 *         Date: 31/08/2011
 */

public class AlloppMoveLegsParser extends AbstractXMLObjectParser {
	public static final String MOVE_LEGS = "moveLegs";


	public String getParserName() {
		return MOVE_LEGS;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
		AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

	    final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
	    return new AlloppMoveLegs(apspnet, apsp, weight);
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
		return "Operator which changes the way a tetraploid subtree joins the diploid tree.";

	}

	@Override
	public Class getReturnType() {
		return AlloppMoveLegs.class;
	}

}



