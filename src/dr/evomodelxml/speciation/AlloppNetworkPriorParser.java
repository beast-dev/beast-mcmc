package dr.evomodelxml.speciation;

import dr.evolution.util.Units;
import dr.evomodel.speciation.AlloppNetworkPrior;
import dr.evomodel.speciation.AlloppNetworkPriorModel;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.MultiSpeciesCoalescent;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
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


public class AlloppNetworkPriorParser extends AbstractXMLObjectParser {

	public static final String APSPNETWORKPRIOR = "apspNetworkPrior";
	public static final String UNITS = "units";
	public static final String MODEL = "model";
	public static final String APSPNETWORK = "apspNetwork";

	@Override
	public String getParserName() {
		return APSPNETWORKPRIOR;
	}
	
	
/*
	<apspNetworkPrior id="apspnetwork.prior" units="substitutions">
    <eventRate>
      <parameter id="apspnetwork.prior.eventRate" value="0.1" lower="0.0" upper="Infinity"/>
    </eventRate>
    <apspNetwork>
      <apspNetwork idref="apspnetwork"/>
		</apspNetwork>
	</apspNetworkPrior>
*/
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final XMLObject mxo = xo.getChild(MODEL);
		final AlloppNetworkPriorModel anpm = (AlloppNetworkPriorModel) mxo.getChild(AlloppNetworkPriorModel.class);
		final XMLObject asnmxo = xo.getChild(APSPNETWORK);
		final AlloppSpeciesNetworkModel asnm = (AlloppSpeciesNetworkModel) asnmxo.getChild(AlloppSpeciesNetworkModel.class);
		return new AlloppNetworkPrior(anpm, asnm);	
	}
	
	private  XMLSyntaxRule[] modelRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(AlloppNetworkPriorModel.class)
		};

	}

	private  XMLSyntaxRule[] asnmRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(AlloppSpeciesNetworkModel.class)
		};

	}
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(MODEL, modelRules()),
				new ElementRule(APSPNETWORK, asnmRules()),
				
		};
	}

	@Override
	public String getParserDescription() {
		return "Prior for an allopolyploid species network.";
	}

	@Override
	public Class getReturnType() {
		return AlloppNetworkPrior.class;
	}

}
