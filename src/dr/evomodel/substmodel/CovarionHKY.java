/*
 * CovarionHKY.java
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

import dr.evolution.datatype.HiddenNucleotides;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A model with hidden states that represent different rates.
 *
 * @author Alexei Drummond
 * @version $Id: CovarionHKY.java,v 1.4 2005/05/24 20:25:58 rambaut Exp $
 */
public class CovarionHKY extends AbstractCovarionDNAModel {
    public static final String COVARION_HKY = "CovarionHKYModel";
    public static final String KAPPA = "kappa";

    /**
     * kappa
     */
    private Parameter kappaParameter;

    /**
     * @param dataType         the datatype to be used
     * @param kappaParameter   the rate of transitions versus transversions
     * @param hiddenClassRates the relative rates of the hidden categories
     *                         (first hidden category has rate 1.0 so this parameter
     *                         has dimension one less than number of hidden categories.
     *                         each hidden category.
     * @param switchingRates   rate of switching between hidden categories
     * @param freqModel        the frequencies
     */
    public CovarionHKY(HiddenNucleotides dataType, Parameter kappaParameter, Parameter hiddenClassRates, Parameter switchingRates, FrequencyModel freqModel) {

        super(COVARION_HKY, dataType, hiddenClassRates, switchingRates, freqModel);

        this.kappaParameter = kappaParameter;
        addParameter(kappaParameter);
        setupRelativeRates();
    }

    double[] getRelativeDNARates() {
        double kappa = kappaParameter.getParameterValue(0);
        return new double[]{1.0, kappa, 1.0, 1.0, kappa, 1.0};
    }

    /**
     * set kappa
     *
     * @param kappa the new value of kappa to set
     */
    public void setKappa(double kappa) {
        kappaParameter.setParameterValue(0, kappa);
    }

    /**
     * @return kappa
     */
    public double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    public String toString() {
        StringBuilder builder = new StringBuilder();
        builder.append("Covarion HKY model with ");
        builder.append(getHiddenClassCount()).append(" rate classes.\n");
        builder.append("Relative rates: \n");
        builder.append(SubstitutionModelUtils.toString(relativeRates, dataType, true, 2));
        return builder.toString();

    }

    /**
     * Parses an element from an DOM document into a DemographicModel. Recognises
     * ConstantPopulation and ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COVARION_HKY;
        }

        public String getParserDescription() {
            return "A covarion HKY model.";
        }

        public Class getReturnType() {
            return SubstitutionModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(KAPPA, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class,
                                "A parameter representing the transition transversion bias")}),
                new ElementRule(SWITCHING_RATES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class,
                                "A parameter representing the rate of change between the different classes")}),
                new ElementRule(HIDDEN_CLASS_RATES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class,
                                "A parameter representing the rates of the hidden classes relative to the first hidden class.")})
        };

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter kappaParam;
            Parameter switchingRates;
            Parameter hiddenClassRates;
            FrequencyModel freqModel;

            kappaParam = (Parameter) xo.getSocketChild(KAPPA);
            switchingRates = (Parameter) xo.getSocketChild(SWITCHING_RATES);
            hiddenClassRates = (Parameter) xo.getSocketChild(HIDDEN_CLASS_RATES);
            freqModel = (FrequencyModel) xo.getSocketChild(FREQUENCIES);

            if (!(freqModel.getDataType() instanceof HiddenNucleotides)) {
                throw new IllegalArgumentException("Datatype must be hidden nucleotides!!");
            }

            HiddenNucleotides dataType = (HiddenNucleotides) freqModel.getDataType();

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
    };
}