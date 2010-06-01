package dr.evomodelxml.substmodel;

import dr.evolution.datatype.*;
import dr.evomodel.substmodel.*;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Parses a GeneralSubstitutionModel or one of its more specific descendants.
 */
public class GeneralSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String GENERAL_SUBSTITUTION_MODEL = "generalSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String RATES = "rates";
    public static final String RELATIVE_TO = "relativeTo";
    public static final String FREQUENCIES = "frequencies";
    public static final String INDICATOR = "rateIndicator";
    public static final String ROOT_FREQ = "rootFrequencies";

    public String getParserName() {
        return GENERAL_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter;
        Parameter indicatorParameter;

        XMLObject cxo = xo.getChild(FREQUENCIES);
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = DataTypeUtils.getDataType(xo);

        if (dataType == null) dataType = (DataType) xo.getChild(DataType.class);

//        if (xo.hasAttribute(DataType.DATA_TYPE)) {
//            String dataTypeStr = xo.getStringAttribute(DataType.DATA_TYPE);
//            if (dataTypeStr.equals(Nucleotides.DESCRIPTION)) {
//                dataType = Nucleotides.INSTANCE;
//            } else if (dataTypeStr.equals(AminoAcids.DESCRIPTION)) {
//                dataType = AminoAcids.INSTANCE;
//            } else if (dataTypeStr.equals(Codons.DESCRIPTION)) {
//                dataType = Codons.UNIVERSAL;
//            } else if (dataTypeStr.equals(TwoStates.DESCRIPTION)) {
//                dataType = TwoStates.INSTANCE;
//            }
//        }

//        if (dataType == null) dataType = freqModel.getDataType();

        if (dataType != freqModel.getDataType()) {
            throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
        }

        cxo = xo.getChild(RATES);
        ratesParameter = (Parameter) cxo.getChild(Parameter.class);

        int states = dataType.getStateCount();
        Logger.getLogger("dr.evomodel").info("  General Substitution Model (stateCount=" + states + ")");
        int rateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount()) / 2;

        if (xo.hasChildNamed(INDICATOR)) {// has indicator
            if (ratesParameter.getDimension() != rateCount) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (rateCount)
                        + " dimensions.  However parameter dimension is " + ratesParameter.getDimension());
            }

            cxo = xo.getChild(INDICATOR);
            indicatorParameter = (Parameter) cxo.getChild(Parameter.class);

            if (indicatorParameter == null || ratesParameter == null || indicatorParameter.getDimension() != ratesParameter.getDimension())
                throw new XMLParseException("Rates and indicator parameters in " + getParserName() + " element must be the same dimension.");

            if (xo.hasChildNamed(ROOT_FREQ)) {
                cxo = xo.getChild(ROOT_FREQ);
                FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);

                if (dataType != rootFreq.getDataType()) {
                    throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
                }
                Logger.getLogger("dr.evomodel").info("  SVS Irreversible Substitution Model is applied with "
                        + INDICATOR + " " + indicatorParameter.getParameterName() + " and " + ROOT_FREQ + " " + rootFreq.getId());
                return new SVSIrreversibleSubstitutionModel(dataType, freqModel, rootFreq, ratesParameter, indicatorParameter);

            }
            Logger.getLogger("dr.evomodel").info("  SVS General Substitution Model is applied with "
                    + INDICATOR + " " + indicatorParameter.getParameterName());
            return new SVSGeneralSubstitutionModel(dataType, freqModel, ratesParameter, indicatorParameter);

        } else {// no indicator
            if (!cxo.hasAttribute(RELATIVE_TO)) {
                throw new XMLParseException("The index of the implicit rate (value 1.0) that all other rates are relative to."
                        + " In DNA this is usually G<->T (6)");
            }
            int relativeTo = cxo.getIntegerAttribute(RELATIVE_TO) - 1;
            if (relativeTo < 0) {
                throw new XMLParseException(RELATIVE_TO + " must be 1 or greater");
            } else {
                int t = relativeTo;
                int s = states - 1;
                int row = 0;
                while (t >= s) {
                    t -= s;
                    s -= 1;
                    row += 1;
                }
                int col = t + row + 1;

                Logger.getLogger("dr.evomodel").info("  Rates relative to "
                        + dataType.getCode(row) + "<->" + dataType.getCode(col));
            }

            if (ratesParameter == null) {
                if (rateCount == 1) {
                    // simplest model for binary traits...
                } else {
                    throw new XMLParseException("No rates parameter found in " + getParserName());
                }
            } else if (ratesParameter.getDimension() != (rateCount - 1)) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have ("
                        + rateCount  + "- 1) dimensions because one of dimensions is fixed.");
            }

            return new GeneralSubstitutionModel(dataType, freqModel, ratesParameter, relativeTo);
        }
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general reversible model of sequence substitution for any data type.";
    }

    public Class getReturnType() {
        return GeneralSubstitutionModelParser.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class)
            ),
            new ElementRule(FREQUENCIES, FrequencyModel.class),
            new ElementRule(RATES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)}
            ),
            new ElementRule(INDICATOR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            new ElementRule(ROOT_FREQ,
                    new XMLSyntaxRule[]{
                            new ElementRule(FrequencyModel.class)
                    }, 0, 1)
    };
}
