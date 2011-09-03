package dr.evomodelxml.operators;

import dr.evomodel.operators.AlloppNetworkNodeSlide;
import dr.evomodel.operators.TreeNodeSlide;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
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

/*

<networkNodeReHeight weight="94">
<alloppSpecies idref="alloppSpecies"/>
<alloppSpeciesNetwork idref="apspnetwork"/>
</networkNodeReHeight> 

*/

public class AlloppNetworkNodeSlideParser extends AbstractXMLObjectParser {
	public static final String NETWORK_NODE_REHEIGHT = "networkNodeReHeight";


	public String getParserName() {
		return NETWORK_NODE_REHEIGHT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
		AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new AlloppNetworkNodeSlide(apspnet, apsp, weight);
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
		return "Operator for allopolyploid species network: transforms network without breaking embedding of gene trees.";
	}

	@Override
	public Class getReturnType() {
		return AlloppNetworkNodeSlide.class;
	}

}
