package dr.app.beagle.tools.parsers;

import dr.app.beagle.evomodel.sitemodel.BranchSubstitutionModel;
import dr.app.beagle.evomodel.sitemodel.GammaSiteRateModel;
import dr.app.beagle.evomodel.sitemodel.HomogenousBranchSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class BeagleSequenceSimulatorParser extends AbstractXMLObjectParser {

	public static final String BEAGLE_SEQUENCE_SIMULATOR = "beagleSequenceSimulator";
	public static final String REPLICATIONS = "replications";
	public static final String SITE_MODEL = "siteModel";
	
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
		
//		return new XMLSyntaxRule[] {
//				new ElementRule(TreeModel.class),
//				new ElementRule(BranchSubstitutionModel.class),
//				new ElementRule(GammaSiteRateModel.class),
////				 new ElementRule(SITE_MODEL, new XMLSyntaxRule[]{new ElementRule(BranchSubstitutionModel.class)}, false),
//				new ElementRule(BranchRateModel.class, true),
//				new ElementRule(FrequencyModel.class),
//				new ElementRule(Sequence.class, true),
//				AttributeRule.newIntegerRule(REPLICATIONS) };
		
		return null;
		
	}//END: getSyntaxRules

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
		
		GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
		FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);
		Sequence ancestralSequence = (Sequence) xo.getChild(Sequence.class);
		int replications = xo.getIntegerAttribute(REPLICATIONS);

		BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
		if (rateModel == null) {
			rateModel = new DefaultBranchRateModel();
		}
		
		// TODO check this for cast errors
		BranchSubstitutionModel branchSubstitutionModel = (BranchSubstitutionModel) xo.getChild(BranchSubstitutionModel.class);
		if (branchSubstitutionModel == null) {
			SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
			branchSubstitutionModel = new HomogenousBranchSubstitutionModel(substitutionModel, freqModel);
		}
		
		BeagleSequenceSimulator s = new BeagleSequenceSimulator(tree, //
				branchSubstitutionModel, //
				siteModel, //
				rateModel, //
				freqModel, //
				replications //
		);

		if (ancestralSequence != null) {

			if (ancestralSequence.getLength() != 3 * replications && freqModel.getDataType() instanceof Codons) {

				throw new RuntimeException("Ancestral codon sequence has "
						+ ancestralSequence.getLength() + " characters "
						+ "expecting " + 3*replications + " characters");

			} else if (ancestralSequence.getLength() != replications && freqModel.getDataType() instanceof Nucleotides) {

				throw new RuntimeException("Ancestral nuleotide sequence has "
						+ ancestralSequence.getLength() + " characters "
						+ "expecting " + replications + " characters");

			} else {

				s.setAncestralSequence(ancestralSequence);

			}// END: dataType check
		}// END: ancestralSequence check

		return s.simulate();

	}// END: parseXMLObject

}// END: class
