/*
 * ComplexSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodelxml.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.*;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class ComplexSubstitutionModelParser extends AbstractXMLObjectParser {

    public static final String COMPLEX_SUBSTITUTION_MODEL = "complexSubstitutionModel";
    public static final String SVS_COMPLEX_SUBSTITUTION_MODEL = "svsComplexSubstitutionModel";
    public static final String DATA_TYPE = "dataType";
    public static final String RATES = "rates";
    public static final String FREQUENCIES = "frequencies";
    public static final String ROOT_FREQUENCIES = "rootFrequencies";
    public static final String RANDOMIZE = "randomizeIndicator";
    public static final String INDICATOR = "rateIndicator";
    public static final String BSSVS_TOLERANCE = "bssvsTolerance";
    public static final String BSSVS_SCALAR = "bssvsScalar";
    public static final String CHECK_CONDITIONING = "checkConditioning";
    public static final String NORMALIZED = "normalized";

    public static final int maxRandomizationTries = 100;

    public String getParserName() {
        return COMPLEX_SUBSTITUTION_MODEL;
    }

    public String[] getParserNames() {
        return new String[]{COMPLEX_SUBSTITUTION_MODEL, SVS_COMPLEX_SUBSTITUTION_MODEL};
    }


    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter ratesParameter;

        XMLObject cxo;
        if (xo.hasChildNamed(FREQUENCIES)) {
            cxo = xo.getChild(FREQUENCIES);
        } else {
            cxo = xo.getChild(ROOT_FREQUENCIES);
        }
        FrequencyModel freqModel = (FrequencyModel) cxo.getChild(FrequencyModel.class);

        DataType dataType = freqModel.getDataType();

        cxo = xo.getChild(RATES);

        int states = dataType.getStateCount();

        Logger.getLogger("dr.app.beagle.evomodel").info("  Complex Substitution Model (stateCount=" + states + ")");

        ratesParameter = (Parameter) cxo.getChild(Parameter.class);

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

        boolean checkConditioning = xo.getAttribute(CHECK_CONDITIONING, true);

        if (!xo.hasChildNamed(INDICATOR)) {
            if (!checkConditioning) {
                return new ComplexSubstitutionModel(COMPLEX_SUBSTITUTION_MODEL, dataType, freqModel, ratesParameter) {
                    protected EigenSystem getDefaultEigenSystem(int stateCount) {
                        return new ComplexColtEigenSystem(stateCount, false, ColtEigenSystem.defaultMaxConditionNumber, ColtEigenSystem.defaultMaxIterations);
                    }
                };
            } else {
                return new ComplexSubstitutionModel(COMPLEX_SUBSTITUTION_MODEL, dataType, freqModel, ratesParameter);
            }
        }

        cxo = xo.getChild(INDICATOR);

        Parameter indicatorParameter = (Parameter) cxo.getChild(Parameter.class);
        if (indicatorParameter == null || ratesParameter == null || indicatorParameter.getDimension() != ratesParameter.getDimension())
            throw new XMLParseException("Rates and indicator parameters in " + getParserName() + " element must be the same dimension.");

        if (xo.hasAttribute(BSSVS_TOLERANCE)) {
            double tolerance = xo.getAttribute(BSSVS_TOLERANCE,
                    BayesianStochasticSearchVariableSelection.Utils.getTolerance());
            if (tolerance > BayesianStochasticSearchVariableSelection.Utils.getTolerance()) {
                // Only increase smallest allowed tolerance
                BayesianStochasticSearchVariableSelection.Utils.setTolerance(tolerance);
                Logger.getLogger("dr.app.beagle.evomodel").info("\tIncreasing BSSVS tolerance to " + tolerance);
            }
        }

        if (xo.hasAttribute(BSSVS_SCALAR)) {
            double scalar = xo.getAttribute(BSSVS_SCALAR,
                    BayesianStochasticSearchVariableSelection.Utils.getScalar());
            if (scalar < BayesianStochasticSearchVariableSelection.Utils.getScalar()) {
                BayesianStochasticSearchVariableSelection.Utils.setScalar(scalar);
                Logger.getLogger("dr.app.beagle.evomodel").info("\tDecreasing BSSVS scalar to " + scalar);
            }
        }

        SVSComplexSubstitutionModel model;
        if (!checkConditioning) {
            model = new SVSComplexSubstitutionModel(SVS_COMPLEX_SUBSTITUTION_MODEL, dataType, freqModel, ratesParameter, indicatorParameter) {
                protected EigenSystem getDefaultEigenSystem(int stateCount) {
                    return new ComplexColtEigenSystem(stateCount, false, ColtEigenSystem.defaultMaxConditionNumber, ColtEigenSystem.defaultMaxIterations);
                }
            };
        } else {
            model = new SVSComplexSubstitutionModel(SVS_COMPLEX_SUBSTITUTION_MODEL, dataType, freqModel, ratesParameter, indicatorParameter);
        }
        boolean randomize = xo.getAttribute(RANDOMIZE, false);
        if (randomize) {
            // Randomization may need multiple tries
            int tries = 0;
            boolean valid = false;

            while (!valid && tries < maxRandomizationTries) {
                BayesianStochasticSearchVariableSelection.Utils.randomize(indicatorParameter,
                        dataType.getStateCount(), false);
                valid = !Double.isInfinite(model.getLogLikelihood());
                tries++;
            }
            Logger.getLogger("dr.app.beagle.evomodel").info("\tRandomization attempts: " + tries);
        }
        if (!xo.getAttribute(NORMALIZED, true)) {
            model.setNormalization(false);
            Logger.getLogger("dr.app.beagle.evomodel").info("\tNormalization: false");
        }
        Logger.getLogger("dr.app.beagle.evomodel").info("\t\tPlease cite: Edwards, Suchard et al. (2011)\n");
        return model;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A general irreversible model of sequence substitution for any data type.";
    }

    public Class getReturnType() {
        return SVSComplexSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new XORRule(
                    new StringAttributeRule(DataType.DATA_TYPE, "The type of sequence data",
                            DataType.getRegisteredDataTypeNames(), false),
                    new ElementRule(DataType.class),
                    true // Optional
            ),
            AttributeRule.newBooleanRule(RANDOMIZE, true),
            new XORRule(
                    new ElementRule(FREQUENCIES, FrequencyModel.class),
                    new ElementRule(ROOT_FREQUENCIES, FrequencyModel.class)),
            new ElementRule(RATES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(INDICATOR,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }, true),
            AttributeRule.newDoubleRule(BSSVS_TOLERANCE, true),
            AttributeRule.newDoubleRule(BSSVS_SCALAR, true),
            AttributeRule.newBooleanRule(CHECK_CONDITIONING, true),
            AttributeRule.newBooleanRule(NORMALIZED, true),
    };
}
