package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.*;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a DemographicModel. Recognises
 * ConstantPopulation and ExponentialGrowth.
 */
public class EmpiricalAminoAcidModelParser extends AbstractXMLObjectParser {

    public static final String EMPIRICAL_AMINO_ACID_MODEL = "aminoAcidModel";
    public static final String FREQUENCIES = "frequencies";
    public static final String TYPE = "type";


    public String getParserName() {
        return EMPIRICAL_AMINO_ACID_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        FrequencyModel freqModel = null;

        if (xo.hasAttribute(FREQUENCIES)) {
            XMLObject cxo = xo.getChild(FREQUENCIES);
            freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        }

        EmpiricalRateMatrix rateMatrix = null;

        String type = xo.getStringAttribute(TYPE);

        if (type.equals(AminoAcidModelType.BLOSUM_62.getXMLName())) {
            rateMatrix = Blosum62.INSTANCE;
        } else if (type.equals(AminoAcidModelType.DAYHOFF.getXMLName())) {
            rateMatrix = Dayhoff.INSTANCE;
        } else if (type.equals(AminoAcidModelType.JTT.getXMLName())) {
            rateMatrix = dr.evomodel.substmodel.JTT.INSTANCE;
        } else if (type.equals(AminoAcidModelType.MT_REV_24.getXMLName())) {
            rateMatrix = MTREV.INSTANCE;
        } else if (type.equals(AminoAcidModelType.CP_REV_45.getXMLName())) {
            rateMatrix = CPREV.INSTANCE;
        } else if (type.equals(AminoAcidModelType.WAG.getXMLName())) {
            rateMatrix = dr.evomodel.substmodel.WAG.INSTANCE;
        } else if (type.equals(AminoAcidModelType.FLU.getXMLName())) {
            rateMatrix = dr.evomodel.substmodel.FLU.INSTANCE;
        }

        return new EmpiricalAminoAcidModel(rateMatrix, freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(TYPE, "The type of empirical amino-acid rate matrix", AminoAcidModelType.xmlNames(), false),
            new ElementRule(FREQUENCIES, FrequencyModel.class, "If the frequencies are omitted than the empirical frequencies associated with the selected model are used.", true)
    };

    public String getParserDescription() {
        return "An empirical amino acid substitution model.";
    }

    public Class getReturnType() {
        return EmpiricalAminoAcidModel.class;
    }


}
