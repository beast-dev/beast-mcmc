package dr.evomodelxml.substmodel;

import dr.evolution.datatype.HiddenNucleotides;
import dr.evomodel.substmodel.AbstractCovarionDNAModel;
import dr.evomodel.substmodel.CovarionGTR;
import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a TwoStateCovarionModel
 */
public class CovarionGTRParser extends AbstractXMLObjectParser {
    public static final String GTR_COVARION_MODEL = "gtrCovarionModel";    

    public String getParserName() {
        return GTR_COVARION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(AbstractCovarionDNAModel.FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        HiddenNucleotides dataType = (HiddenNucleotides) freqModel.getDataType();

        Parameter hiddenRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES);
        Parameter switchingRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.SWITCHING_RATES);

        Parameter rateACParameter = null;
        if (xo.hasChildNamed(GTRParser.A_TO_C)) {
            rateACParameter = (Parameter) xo.getElementFirstChild(GTRParser.A_TO_C);
        }
        Parameter rateAGParameter = null;
        if (xo.hasChildNamed(GTRParser.A_TO_G)) {
            rateAGParameter = (Parameter) xo.getElementFirstChild(GTRParser.A_TO_G);
        }
        Parameter rateATParameter = null;
        if (xo.hasChildNamed(GTRParser.A_TO_T)) {
            rateATParameter = (Parameter) xo.getElementFirstChild(GTRParser.A_TO_T);
        }
        Parameter rateCGParameter = null;
        if (xo.hasChildNamed(GTRParser.C_TO_G)) {
            rateCGParameter = (Parameter) xo.getElementFirstChild(GTRParser.C_TO_G);
        }
        Parameter rateCTParameter = null;
        if (xo.hasChildNamed(GTRParser.C_TO_T)) {
            rateCTParameter = (Parameter) xo.getElementFirstChild(GTRParser.C_TO_T);
        }
        Parameter rateGTParameter = null;
        if (xo.hasChildNamed(GTRParser.G_TO_T)) {
            rateGTParameter = (Parameter) xo.getElementFirstChild(GTRParser.G_TO_T);
        }


        if (dataType != freqModel.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
        }

        return new CovarionGTR(
                dataType,
                hiddenRates,
                switchingRates,
                rateACParameter,
                rateAGParameter,
                rateATParameter,
                rateCGParameter,
                rateCTParameter,
                rateGTParameter,
                freqModel);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A covarion substitution model of langauge evolution with binary data and a hidden rate state with two rates.";
    }

    public Class getReturnType() {
        return CovarionGTR.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(GTRParser.A_TO_C, Parameter.class, "relative rate of A<->C substitution", true),
            new ElementRule(GTRParser.A_TO_G, Parameter.class, "relative rate of A<->G substitution", true),
            new ElementRule(GTRParser.A_TO_T, Parameter.class, "relative rate of A<->T substitution", true),
            new ElementRule(GTRParser.C_TO_G, Parameter.class, "relative rate of C<->G substitution", true),
            new ElementRule(GTRParser.C_TO_T, Parameter.class, "relative rate of C<->T substitution", true),
            new ElementRule(GTRParser.G_TO_T, Parameter.class, "relative rate of G<->T substitution", true),
            new ElementRule(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES, Parameter.class),
            new ElementRule(AbstractCovarionDNAModel.SWITCHING_RATES, Parameter.class),
            new ElementRule(AbstractCovarionDNAModel.FREQUENCIES, FrequencyModel.class),
    };



}
