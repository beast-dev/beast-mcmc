package dr.evomodelxml.speciation;


import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * 
 * Parses a possibly-allopolyploid species made of individual organisms. 
 * 
 * @author Graham Jones
 *         Date: 18/04/2011
 */



/*
 * 
 * Parses a diploid or allopolyploid species, recording ploidy level
 * and running through all individual organisms belonging to the species. 
 * 
 * Part of parsing a AlloppSpeciesBindings.
 * 
 */


public class AlloppSpeciesBindingsApSpInfoParser extends
		AbstractXMLObjectParser {
	public static final String APSP = "apsp";
	public static final String PLOIDYLEVEL = "ploidylevel";

	public String getParserName() {
		return APSP;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		AlloppSpeciesBindings.Individual[] individuals = 
			new AlloppSpeciesBindings.Individual[xo.getChildCount()];
        for (int ni = 0; ni < individuals.length; ++ni) {
        	individuals[ni] = (AlloppSpeciesBindings.Individual) xo.getChild(ni);
        }
        return new AlloppSpeciesBindings.ApSpInfo(xo.getId(), xo.getIntegerAttribute(PLOIDYLEVEL), individuals);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] {
			AttributeRule.newDoubleRule(PLOIDYLEVEL),
			new ElementRule(AlloppSpeciesBindings.Individual.class, 1, Integer.MAX_VALUE)
			};
	}

	@Override
	public String getParserDescription() {
		return "A diploid or allopolyploid species made of individuals";
	}

	@Override
	public Class getReturnType() {
		return AlloppSpeciesBindings.ApSpInfo.class;
	}

}
