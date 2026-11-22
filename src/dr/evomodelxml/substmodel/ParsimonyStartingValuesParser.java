/*
 * ParsimonyStartingValuesParser.java
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

import dr.evolution.alignment.PatternList;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.ParsimonyStartingValues;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Parameter;
import dr.util.Transform;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Sets parsimony-based starting values for the rates, BSSVS indicators and clock rate of an asymmetric
 * (ComplexSubstitutionModel) or reversible (GeneralSubstitutionModel) CTMC. For a transformedRateSubstitutionModel
 * (LogRateSubstitutionModel) the transformed-rate parameter is set to the inverse transform of each rate.
 * <pre>
 * &lt;parsimonyStartingValues setClock="true" setIndicators="true" setRelativeRates="true" unobservedRate="1E-4"&gt;
 *     &lt;treeDataLikelihood idref="treeDataLikelihood"/&gt;
 *     &lt;generalSubstitutionModel idref="sm"/&gt;
 *     &lt;parameter idref="clock.rate"/&gt;    (only for setClock="true")
 *     &lt;attributePatterns idref="patterns"/&gt; (optional; default: the likelihood's BEAGLE patterns)
 * &lt;/parsimonyStartingValues&gt;
 * </pre>
 * unobservedRate (optional) is the starting rate of every pair with no parsimony changes, relative to the observed
 * rates, which average 1; by default it is 1 / (mean changes per observed pair), as if each had one change.
 *
 * @author Marc A. Suchard
 */
public class ParsimonyStartingValuesParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "parsimonyStartingValues";
    private static final String SET_CLOCK = "setClock";
    private static final String SET_INDICATORS = "setIndicators";
    private static final String SET_RELATIVE_RATES = "setRelativeRates";
    private static final String UNOBSERVED_RATE = "unobservedRate";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        boolean setClock = xo.getAttribute(SET_CLOCK, true);
        boolean setIndicators = xo.getAttribute(SET_INDICATORS, true);
        boolean setRates = xo.getAttribute(SET_RELATIVE_RATES, true);
        double unobservedRate = xo.getAttribute(UNOBSERVED_RATE, Double.NaN);
        if (!Double.isNaN(unobservedRate) && !(unobservedRate > 0.0 && unobservedRate < Double.POSITIVE_INFINITY)) {
            throw new XMLParseException(UNOBSERVED_RATE + " must be positive and finite");
        }

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        PatternList patterns = (PatternList) xo.getChild(PatternList.class);
        if (patterns == null) {
            DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
            if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
                throw new XMLParseException("Cannot reach the patterns of " + treeDataLikelihood.getId() +
                        "; add the pattern list to " + getParserName());
            }
            patterns = ((BeagleDataLikelihoodDelegate) delegate).getPatternList();
        }

        GeneralSubstitutionModel substitutionModel =
                (GeneralSubstitutionModel) xo.getChild(GeneralSubstitutionModel.class);
        boolean reversible = !(substitutionModel instanceof ComplexSubstitutionModel);

        Parameter indicators = null;
        if (setIndicators) {
            if (!(substitutionModel instanceof BayesianStochasticSearchVariableSelection)) {
                throw new XMLParseException(SET_INDICATORS + "=\"true\" needs a substitution model with a " +
                        "rateIndicator");
            }
            indicators = ((BayesianStochasticSearchVariableSelection) substitutionModel).getIndicators();
        }

        Parameter rates = null;
        Transform transform = null;
        if (setRates) {
            if (substitutionModel instanceof LogRateSubstitutionModel) {
                LogAdditiveCtmcRateProvider provider =
                        ((LogRateSubstitutionModel) substitutionModel).getRateProvider();
                if (!(provider instanceof LogAdditiveCtmcRateProvider.DataAugmented)) {
                    throw new XMLParseException("Substitution model " + substitutionModel.getId() +
                            " has no rate parameter to set (e.g. a GLM)");
                }
                rates = ((LogAdditiveCtmcRateProvider.DataAugmented) provider).getLogRateParameter();
                transform = provider.getTransform() != null ? provider.getTransform() : Transform.LOG;
            } else {
                rates = substitutionModel.getRatesParameter();
            }
            if (rates == null) {
                throw new XMLParseException("Substitution model " + substitutionModel.getId() +
                        " has no rates parameter to set");
            }
            int stateCount = substitutionModel.getDataType().getStateCount();
            int expected = reversible ? stateCount * (stateCount - 1) / 2 : stateCount * (stateCount - 1);
            if (rates.getDimension() != expected) {
                throw new XMLParseException("Rates parameter " + rates.getId() + " has dimension " +
                        rates.getDimension() + "; " + getParserName() + " expects " + expected +
                        (reversible ? " (reversible models with a rate fixed by relativeTo are not supported)" : ""));
            }
        }

        Parameter clockRate = (Parameter) xo.getChild(Parameter.class);
        if (setClock && clockRate == null) {
            throw new XMLParseException(SET_CLOCK + "=\"true\" needs a clock-rate parameter");
        }

        ParsimonyStartingValues values;
        try {
            values = new ParsimonyStartingValues(treeDataLikelihood.getTree(), patterns,
                    substitutionModel.getDataType());
            if (indicators != null) {
                values.setIndicators(indicators, reversible);
            }
            if (rates != null) {
                values.setRates(rates, setIndicators, reversible, transform, unobservedRate);
                if (!substitutionModel.getNormalization()) {
                    // An unnormalized Q carries the rate scale: make clock * (rate of Q) the parsimony clock rate
                    double clock = setClock ? values.getClockRate() :
                            (clockRate != null ? clockRate.getParameterValue(0) : 1.0);
                    values.scaleRates(rates, values.getClockRate() / (clock * rateOf(substitutionModel)), transform);
                }
            }
            if (setClock) {
                values.setClockRate(clockRate);
            }
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }

        Logger.getLogger("dr.evomodel").info("Parsimony starting values (" +
                describe(setRates, setIndicators, setClock) + "): " + values.getSummary());

        return values;
    }

    // Expected substitutions per unit time of an unnormalized Q, with the weights the model would normalize by
    private static double rateOf(GeneralSubstitutionModel model) {
        int stateCount = model.getDataType().getStateCount();
        double[] q = new double[stateCount * stateCount];
        model.getInfinitesimalMatrix(q);
        double[] weights = model instanceof ComplexSubstitutionModel ?
                ((ComplexSubstitutionModel) model).getNormalizationFrequencies() :
                model.getFrequencyModel().getFrequencies();
        double rate = 0.0;
        for (int i = 0; i < stateCount; ++i) {
            rate -= weights[i] * q[i * stateCount + i];
        }
        return rate;
    }

    private static String describe(boolean rates, boolean indicators, boolean clock) {
        StringBuilder sb = new StringBuilder();
        if (rates) sb.append("rates");
        if (indicators) sb.append(sb.length() > 0 ? ", " : "").append("indicators");
        if (clock) sb.append(sb.length() > 0 ? ", " : "").append("clock rate");
        return sb.length() > 0 ? "set " + sb : "nothing set";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newBooleanRule(SET_CLOCK, true),
            AttributeRule.newBooleanRule(SET_INDICATORS, true),
            AttributeRule.newBooleanRule(SET_RELATIVE_RATES, true),
            AttributeRule.newDoubleRule(UNOBSERVED_RATE, true),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(GeneralSubstitutionModel.class),
            new ElementRule(Parameter.class, true),
            new ElementRule(PatternList.class, true),
    };

    public String getParserDescription() {
        return "Sets parsimony-based starting values for the rates, BSSVS indicators and clock rate of an " +
                "asymmetric or reversible CTMC";
    }

    public Class getReturnType() {
        return ParsimonyStartingValues.class;
    }
}
