package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.AminoAcidMixture;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.SubstitutionRateMatrixMixture;
import dr.evomodel.substmodel.aminoacid.EmpiricalAminoAcidModel;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

public class SubstitutionRateMatrixMixtureParser extends AbstractXMLObjectParser {
    private static final String MIXTURE_MODEL = "substitutionRateMatrixMixtureModel";

    public String getParserName() {
        return MIXTURE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<SubstitutionModel> modelList = new ArrayList<>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            SubstitutionModel model = (SubstitutionModel) xo.getChild(i);
            // TODO figure out why logic of this class does not apply to EAAM
            if (model instanceof EmpiricalAminoAcidModel) {throw new RuntimeException("For mixtures of empirical amino acid models use aminoAcidMixtureModel");}
            modelList.add(model);
        }

        return new SubstitutionRateMatrixMixture(modelList);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "List of Q-matrices for combining in a GLM.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SubstitutionModel.class, 1, Integer.MAX_VALUE),
    };
}
