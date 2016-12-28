/*
 * CovarionHKYParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodelxml.substmodel;

import dr.evolution.datatype.OldHiddenNucleotides;
import dr.oldevomodel.substmodel.AbstractCovarionDNAModel;
import dr.oldevomodel.substmodel.CovarionHKY;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a DemographicModel. Recognises
 * ConstantPopulation and ExponentialGrowth.
 */
public class CovarionHKYParser extends AbstractXMLObjectParser {
    public static final String COVARION_HKY = "CovarionHKYModel";
    public static final String KAPPA = HKYParser.KAPPA;

    public String getParserName() {
        return COVARION_HKY;
    }

    public String getParserDescription() {
        return "A covarion HKY model.";
    }

    public Class getReturnType() {
        return CovarionHKY.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(KAPPA, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the transition transversion bias")}),
            new ElementRule(AbstractCovarionDNAModel.SWITCHING_RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the rate of change between the different classes")}),
            new ElementRule(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class,
                            "A parameter representing the rates of the hidden classes relative to the first hidden class.")})
    };

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter kappaParam;
        Parameter switchingRates;
        Parameter hiddenClassRates;
        FrequencyModel freqModel;

        kappaParam = (Parameter) xo.getElementFirstChild(KAPPA);
        switchingRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.SWITCHING_RATES);
        hiddenClassRates = (Parameter) xo.getElementFirstChild(AbstractCovarionDNAModel.HIDDEN_CLASS_RATES);
        freqModel = (FrequencyModel) xo.getElementFirstChild(AbstractCovarionDNAModel.FREQUENCIES);

        if (!(freqModel.getDataType() instanceof OldHiddenNucleotides)) {
            throw new IllegalArgumentException("Datatype must be hidden nucleotides!!");
        }

        OldHiddenNucleotides dataType = (OldHiddenNucleotides) freqModel.getDataType();

        int hiddenStateCount = dataType.getHiddenClassCount();

        int switchingRatesCount = hiddenStateCount * (hiddenStateCount - 1) / 2;

        if (switchingRates.getDimension() != switchingRatesCount) {
            throw new IllegalArgumentException("switching rates parameter must have " +
                    switchingRatesCount + " dimensions, for " + hiddenStateCount +
                    " hidden categories");
        }

        CovarionHKY model = new CovarionHKY(dataType, kappaParam, hiddenClassRates, switchingRates, freqModel);
        System.out.println(model);
        return model;
    }

}
