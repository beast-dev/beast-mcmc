package dr.app.beagle.evomodel.parsers;

import dr.xml.*;
import dr.app.beagle.evomodel.substmodel.FrequencyModel;
import dr.app.beagle.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evolution.datatype.*;
import dr.inference.model.Parameter;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class SVSComplexSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String SVS_COMPLEX_SUBSTITUTION_MODEL = "svsComplexSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String RATES = "rates";
    public static final String FREQUENCIES = "frequencies";
    public static final String INDICATOR = "rateIndicator";
    public static final String ROOT_FREQ = "rootFrequencies";

    public String getParserName() {
        return SVS_COMPLEX_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();

        cxo = xo.getChild(RATES);

        int states = dataType.getStateCount();

        Logger.getLogger("dr.app.beagle.evomodel").info("  BSSVS Complex Substitution Model (stateCount=" + states + ")");

        Parameter ratesParameter = (Parameter) cxo.getChild(Parameter.class);

        int rateCount = (dataType.getStateCount() - 1) * dataType.getStateCount();

        if (ratesParameter == null) {

            if (rateCount == 1) {
                // simplest model for binary traits...
            } else {
                throw new XMLParseException("No rates parameter found in " + getParserName());
            }
        } else if (ratesParameter.getDimension() != rateCount) {
            throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + rateCount + " dimensions.");
        }

        cxo = xo.getChild(INDICATOR);

        Parameter indicatorParameter = (Parameter) cxo.getChild(Parameter.class);

        if (indicatorParameter == null || ratesParameter == null || indicatorParameter.getDimension() != ratesParameter.getDimension())
            throw new XMLParseException("Rates and indicator parameters in " + getParserName() + " element must be the same dimension.");

        if (xo.hasChildNamed(ROOT_FREQ)) {

            cxo = xo.getChild(ROOT_FREQ);
            dr.evomodel.substmodel.FrequencyModel rootFreq = (dr.evomodel.substmodel.FrequencyModel) cxo.getChild(dr.evomodel.substmodel.FrequencyModel.class);

            if (dataType != rootFreq.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
            }

//              return new SVSIrreversibleSubstitutionModel(dataType, freqModel, rootFreq, ratesParameter, indicatorParameter);
            throw new RuntimeException("SVSIrreversibleComplexSubstitutionModel is not yet implemented.");

        }

        return new SVSComplexSubstitutionModel(SVS_COMPLEX_SUBSTITUTION_MODEL,dataType, freqModel, ratesParameter, indicatorParameter);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for any data type.";
    }

    public Class getReturnType() {
        return SubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data", new String[]{Nucleotides.DESCRIPTION, AminoAcids.DESCRIPTION, Codons.DESCRIPTION, TwoStates.DESCRIPTION}, false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(FREQUENCIES,FrequencyModel.class),
            new ElementRule(RATES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)}
            ),
            new ElementRule(INDICATOR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }),
            new ElementRule(ROOT_FREQ,
                    new XMLSyntaxRule[]{
                            new ElementRule(FrequencyModel.class)
                    }, 0, 1)
    };
}
