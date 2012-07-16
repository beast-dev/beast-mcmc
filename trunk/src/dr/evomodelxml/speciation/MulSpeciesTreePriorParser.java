package dr.evomodelxml.speciation;


import dr.evomodel.speciation.MulSpeciesTreeModel;
import dr.evomodel.speciation.MulSpeciesTreePrior;
import dr.evomodel.speciation.SpeciationModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class MulSpeciesTreePriorParser extends AbstractXMLObjectParser {
	public static final String MUL_SPECIES_TREE_PRIOR = "mulSpeciesTreePrior";
	public static final String MODEL = "model";
	public static final String MUL_SPECIES_TREE = "mulTree";


	public String getParserName() {
		return MUL_SPECIES_TREE_PRIOR;
	}


	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final XMLObject mxo = xo.getChild(MODEL);
		final SpeciationModel sppm = (SpeciationModel) mxo.getChild(SpeciationModel.class);
		final XMLObject mulsptxo = xo.getChild(MUL_SPECIES_TREE);
		final MulSpeciesTreeModel mulspt = (MulSpeciesTreeModel) mulsptxo.getChild(MulSpeciesTreeModel.class);
		return new MulSpeciesTreePrior(sppm, mulspt);	
	}

	private  XMLSyntaxRule[] modelRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciationModel.class)
        };
    }

    private  XMLSyntaxRule[] mulsptRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MulSpeciesTreeModel.class)
        };
    }
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(MODEL, modelRules()),
				new ElementRule(MUL_SPECIES_TREE, mulsptRules()),
				
		};
	}	
	


	@Override
	public String getParserDescription() {
		return "Prior for a multiply-labelled species tree for allopolyploids.";
	}

	@Override
	public Class getReturnType() {
		return MulSpeciesTreePrior.class;
	}

}
