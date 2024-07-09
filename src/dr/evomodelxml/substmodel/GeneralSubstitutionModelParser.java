/*
 * GeneralSubstitutionModelParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.evolution.datatype.DataType;
import dr.oldevomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
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
    public static final String NORMALIZED = "normalized";
    public static final String COMPUTE_STATIONARY = "computeStationary";

    public String getParserName() {
        return GENERAL_SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter = null;
        FrequencyModel freqModel = null;

        if (xo.hasChildNamed(FREQUENCIES)) {
            XMLObject cxo = xo.getChild(FREQUENCIES);
            freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);
        } else if (!xo.getAttribute(COMPUTE_STATIONARY, false)) {
            throw new XMLParseException("No frequency model found in " + getParserName());
        }

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

        if (dataType == null && freqModel != null) {
            dataType = freqModel.getDataType(); 
        } else if (freqModel != null) {
           if (dataType != freqModel.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
            }
        }

        if (dataType == null) {
            throw new XMLParseException("Data type of " + getParserName() + " element can not be found.");
        }

        XMLObject cxo = xo.getChild(RATES);
        ratesParameter = (Parameter) cxo.getChild(Parameter.class);

        int states = dataType.getStateCount();
        Logger.getLogger("dr.evomodel").info("  General Substitution Model (stateCount=" + states + ")");

        boolean hasRelativeRates = cxo.hasChildNamed(RELATIVE_TO) || (cxo.hasAttribute(RELATIVE_TO) && cxo.getIntegerAttribute(RELATIVE_TO) > 0);

        int nonReversibleRateCount = ((dataType.getStateCount() - 1) * dataType.getStateCount());
        int reversibleRateCount = (nonReversibleRateCount / 2);

        boolean isNonReversible = ratesParameter.getDimension() == nonReversibleRateCount;
        boolean hasIndicator = xo.hasChildNamed(INDICATOR);

        GeneralSubstitutionModel model;

        if (!hasRelativeRates) {
            Parameter indicatorParameter = null;

            if (ratesParameter.getDimension() != reversibleRateCount && ratesParameter.getDimension() != nonReversibleRateCount) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (reversibleRateCount)
                        + " dimensions for reversible model or " + nonReversibleRateCount + " dimensions for non-reversible. " +
                        "However parameter dimension is " + ratesParameter.getDimension());
            }

            if (hasIndicator) { // this is using BSSVS
                cxo = xo.getChild(INDICATOR);
                indicatorParameter = (Parameter) cxo.getChild(Parameter.class);

                if (indicatorParameter.getDimension() != ratesParameter.getDimension()) {
                    throw new XMLParseException("Rates and indicator parameters in " + getParserName() + " element must be the same dimension.");
                }

                boolean randomize = xo.getAttribute(ComplexSubstitutionModelParser.RANDOMIZE, false);
                if (randomize) {
                    BayesianStochasticSearchVariableSelection.Utils.randomize(indicatorParameter,
                        dataType.getStateCount(), !isNonReversible);
                }

            }

            if (isNonReversible) {
//                if (xo.hasChildNamed(ROOT_FREQ)) {
//                    cxo = xo.getChild(ROOT_FREQ);
//                    FrequencyModel rootFreq = (FrequencyModel) cxo.getChild(FrequencyModel.class);
//
//                    if (dataType != rootFreq.getDataType()) {
//                        throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its rootFrequencyModel.");
//                    }
//
//                    Logger.getLogger("dr.evomodel").info("  Using BSSVS Complex Substitution Model");
//                    return new SVSComplexSubstitutionModel(getParserName(), dataType, freqModel, ratesParameter, indicatorParameter);
//
//                } else {
//                    throw new XMLParseException("Non-reversible model missing " + ROOT_FREQ + " element");
//                }
                Logger.getLogger("dr.evomodel").info("  Using BSSVS Complex Substitution Model");
                model = new SVSComplexSubstitutionModel(getParserName(), dataType, freqModel, ratesParameter, indicatorParameter);
            } else {
                Logger.getLogger("dr.evomodel").info("  Using BSSVS General Substitution Model");
                model = new SVSGeneralSubstitutionModel(getParserName(), dataType, freqModel, ratesParameter, indicatorParameter);
            }


        } else {
            // if we have relativeTo attribute then we use the old GeneralSubstitutionModel

            if (ratesParameter.getDimension() != reversibleRateCount - 1) {
                throw new XMLParseException("Rates parameter in " + getParserName() + " element should have " + (reversibleRateCount - 1)
                        + " dimensions. However parameter dimension is " + ratesParameter.getDimension());
            }

            int relativeTo = 0;
            if (hasRelativeRates) {
                relativeTo = cxo.getIntegerAttribute(RELATIVE_TO) - 1;
            }

            if (relativeTo < 0 || relativeTo >= reversibleRateCount) {
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
                if (reversibleRateCount == 1) {
                    // simplest model for binary traits...
                } else {
                    throw new XMLParseException("No rates parameter found in " + getParserName());
                }
            }

            model = new GeneralSubstitutionModel(getParserName(), dataType, freqModel, ratesParameter, relativeTo);
        }

        if (!xo.getAttribute(NORMALIZED, true)) {
            model.setNormalization(false);
            Logger.getLogger("dr.app.beagle.evomodel").info("\tNormalization: false");
        }

        return model;
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
                    , true),
            new ElementRule(FREQUENCIES, 
                    new XMLSyntaxRule[]{
                            new ElementRule(FrequencyModel.class)}, true),
            new ElementRule(RATES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)}
            ),
            new ElementRule(INDICATOR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class),
                    }, true),
            AttributeRule.newBooleanRule(ComplexSubstitutionModelParser.RANDOMIZE,true),
            AttributeRule.newBooleanRule(NORMALIZED, true),
            AttributeRule.newBooleanRule(COMPUTE_STATIONARY, true),
    };
}
