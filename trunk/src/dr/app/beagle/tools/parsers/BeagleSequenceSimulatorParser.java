package dr.app.beagle.tools.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.evolution.alignment.Alignment;
import dr.evolution.sequence.Sequence;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BeagleSequenceSimulatorParser extends AbstractXMLObjectParser {

	public static final String BEAGLE_SEQUENCE_SIMULATOR = "beagleSequenceSimulator";
	public static final String REPLICATIONS = "replications";

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(TreeModel.class),
			new ElementRule(BranchSubstitutionModel.class),
			new ElementRule(GammaSiteRateModel.class),
			new ElementRule(BranchRateModel.class),
			new ElementRule(FrequencyModel.class),
			new ElementRule(Sequence.class, true),
			AttributeRule.newIntegerRule(REPLICATIONS) };

	@Override
	public String getParserName() {
		return BEAGLE_SEQUENCE_SIMULATOR;
	}

	@Override
	public String getParserDescription() {
		return "Beagle sequence simulator";
	}

	@Override
	public Class<Alignment> getReturnType() {
		return Alignment.class;
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return rules;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
		BranchSubstitutionModel branchSubstitutionModel = (BranchSubstitutionModel) xo.getChild(BranchSubstitutionModel.class);
		GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
		BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
		FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);
		Sequence ancestralSequence = (Sequence) xo.getChild(Sequence.class);
		int sequenceLength = xo.getIntegerAttribute(REPLICATIONS);

		//TODO:
		////////////////////////////////////////////////////////////
		
//        if (xo.hasChildNamed("siteModel")) {
//        	
//            XMLObject cxo = xo.getChild("siteModel");
//    		BranchSubstitutionModel branchSubstitutionModel = (BranchSubstitutionModel) cxo.getChild(BranchSubstitutionModel.class);
//    		siteModel.addModel(branchSubstitutionModel);
//    		
//        } else {
//        	System.err.println("FUBAR");
//        }

		////////////////////////////////////////////////////////////       
        
		if (rateModel == null) {
			rateModel = new DefaultBranchRateModel();
		}
		
		BeagleSequenceSimulator s = new BeagleSequenceSimulator(tree, //
				branchSubstitutionModel, //
				siteModel, //
				rateModel, //
				freqModel, //
				sequenceLength //
		);

		if (ancestralSequence != null) {
			s.setAncestralSequence(ancestralSequence);
		}

		return s.simulate();

	}// END: parseXMLObject

}// END: class
