/*
 * CovarionSubstitutionModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.*;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Helen Shearman
 * @author Alexei Drummond
 *
 * @version $Id$
 */
public class CovarionSubstitutionModel extends GeneralSubstitutionModel {

    public static final String COVARION_MODEL = "covarionModel";
    public static final String ALPHA = "alpha";
    public static final String BETA = "beta";
    public static final String GAMMA = "gamma";

    /**
     * constructor
     *
     * @param dataType the data type
     * @param alphaParameter - the rate of evolution in slow mode
     * @param betaParameter - the rate of flipping between slow and fast modes
     */
    public CovarionSubstitutionModel(TwoStateCovarion dataType, FrequencyModel freqModel,
                                     Parameter alphaParameter,
                                     Parameter betaParameter, Parameter gammaParameter) {
        super(COVARION_MODEL, dataType, freqModel, 5);

        alpha = alphaParameter;
        beta = betaParameter;
        gamma = gammaParameter;

        addParameter(alpha );
        addParameter(beta);
        addParameter(gamma);
    }

    protected void setupRelativeRates() {

        relativeRates[0] = alpha.getParameterValue(0);
        relativeRates[1] = beta.getParameterValue(0);
        relativeRates[2] = 0.0;
        relativeRates[3] = 0.0;
        relativeRates[4] = gamma.getParameterValue(0);
        relativeRates[5] = 1.0;

        //for (int i = 0; i < 5; i++) {
        //    System.out.print(relativeRates[i] + " ");
        //}
        //System.out.println();
    }

    /**
     * Parses an element from an DOM document into a CovarionSubstitutionModel
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return COVARION_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter alphaParameter;
            Parameter betaParameter;
            Parameter gammaParameter;

            XMLObject cxo = (XMLObject)xo.getChild(FREQUENCIES);
            FrequencyModel freqModel = (FrequencyModel)cxo.getChild(FrequencyModel.class);

            TwoStateCovarion dataType = TwoStateCovarion.INSTANCE;  // fancy new datatype courtesy of Helen

            cxo = (XMLObject)xo.getChild(ALPHA);
            alphaParameter = (Parameter)cxo.getChild(Parameter.class);

            // alpha must be positive and less than 1.0 because the fast rate is normalized to 1.0
            alphaParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

            cxo = (XMLObject)xo.getChild(BETA);
            betaParameter = (Parameter)cxo.getChild(Parameter.class);

            cxo = (XMLObject)xo.getChild(GAMMA);
            gammaParameter = (Parameter)cxo.getChild(Parameter.class);

            if (dataType != freqModel.getDataType()) {
                throw new XMLParseException("Data type of " + getParserName() + " element does not match that of its frequencyModel.");
            }

            return new CovarionSubstitutionModel(dataType, freqModel, alphaParameter, betaParameter, gammaParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A covarion substitution model of langauge evolution with binary data and a hidden rate state with two rates.";
        }

        public Class getReturnType() { return CovarionSubstitutionModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(FREQUENCIES, FrequencyModel.class),
            new ElementRule(ALPHA,
                new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(BETA,
                new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class, true)}
            ),
            new ElementRule(GAMMA,
                new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class, true)}
            )
        };

    };


    private Parameter alpha;
    private Parameter beta;
    private Parameter gamma;
}
