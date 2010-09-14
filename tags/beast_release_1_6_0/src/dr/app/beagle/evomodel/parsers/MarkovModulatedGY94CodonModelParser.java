package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.substmodel.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneticCode;
import dr.evolution.datatype.HiddenCodons;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class MarkovModulatedGY94CodonModelParser extends GY94CodonModelParser {

    public static final String MARKOV_MODULATED_YANG_MODEL = "markovModulatedYangCodonModel";
    public static final String HIDDEN_COUNT = "hiddenCount";
    public static final String SWITCHING_RATES = "switchingRates";
    public static final String DIAGONALIZATION = "diagonalization";

    public String getParserName() {
        return MARKOV_MODULATED_YANG_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = DataTypeUtils.getDataType(xo);
        HiddenCodons codons;
        if (dataType instanceof HiddenCodons)
            codons = (HiddenCodons) dataType;
        else
            throw new XMLParseException("Must construct " + MARKOV_MODULATED_YANG_MODEL + " with hidden codons");

        Parameter omegaParam = (Parameter) xo.getElementFirstChild(OMEGA);
        Parameter kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);
        Parameter switchingParam = (Parameter) xo.getElementFirstChild(SWITCHING_RATES);
        FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);

        EigenSystem eigenSystem;
        if (xo.getAttribute(DIAGONALIZATION,"default").compareToIgnoreCase("colt") == 0)
            eigenSystem = new ColtEigenSystem();
        else
            eigenSystem = new DefaultEigenSystem(dataType.getStateCount());

        return new MarkovModulatedGY94CodonModel(codons, switchingParam, omegaParam, kappaParam, freqModel,eigenSystem);
    }

    public String getParserDescription() {
        return "This element represents the a Markov-modulated Yang model of codon evolution.";
    }

    public Class getReturnType() {
        return MarkovModulatedGY94CodonModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{

            AttributeRule.newStringRule(DataType.DATA_TYPE),
            AttributeRule.newStringRule(GeneticCode.GENETIC_CODE),
            new ElementRule(OMEGA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(KAPPA,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(SWITCHING_RATES,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(FrequencyModel.class),
            AttributeRule.newStringRule(DIAGONALIZATION),
    };

}
