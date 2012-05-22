package dr.evomodelxml.speciation;

import dr.evomodel.speciation.AlloppMSCoalescent;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
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

public class AlloppMSCoalescentParser extends AbstractXMLObjectParser {

    public static final String ALLOPPMSCOALESCENT = "apspCoalescent";


    public String getParserName() {
	    return ALLOPPMSCOALESCENT;
	}



	
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final AlloppSpeciesBindings apsp = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
        final AlloppSpeciesNetworkModel apspnetwork = 
        	       (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);
        return new AlloppMSCoalescent(apsp, apspnetwork);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class),
        };
    }

	

	@Override
	public String getParserDescription() {
		return "Likelihood of a set of gene trees embedded in a allopolyploid species network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppMSCoalescent.class;
	}

}
