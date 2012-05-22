package dr.evomodelxml.operators;


import dr.evomodel.operators.MulTreeSequenceReassignment;
import dr.evomodel.speciation.MulSpeciesBindings;
import dr.evomodel.speciation.MulSpeciesTreeModel;
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
 *         Date: 20/12/2011
 */


public class MulTreeSequenceReassignmentParser extends AbstractXMLObjectParser {
	public static final String MULTREE_SEQUENCE_REASSIGNMENT = "mulTreeSequenceReassignment";
	
	
	public String getParserName() {
		return MULTREE_SEQUENCE_REASSIGNMENT;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		MulSpeciesBindings mulspb = (MulSpeciesBindings) xo.getChild(MulSpeciesBindings.class);
		MulSpeciesTreeModel multree = (MulSpeciesTreeModel) xo.getChild(MulSpeciesTreeModel.class);

        final double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);
        return new MulTreeSequenceReassignment(multree, mulspb, weight);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(MulSpeciesBindings.class),
                new ElementRule(MulSpeciesTreeModel.class)
        };
	}

	@Override
	public String getParserDescription() {
		return "Operator which reassigns sequences within an allopolyploid species.";
	}

	@Override
	public Class getReturnType() {
		return MulTreeSequenceReassignment.class;
	}

}
