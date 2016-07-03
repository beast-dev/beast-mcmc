/*
 * HypermutantErrorModel.java
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

package dr.evomodel.tipstatesmodel;

import dr.evolution.alignment.HypermutantAlignment;
import dr.evolution.datatype.Nucleotides;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.*;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class HypermutantErrorModel extends TipStatesModel implements Citable {

    public static final String HYPERMUTANT_ERROR_MODEL = "hypermutantErrorModel";
    public static final String HYPERMUTATION_RATE = "hypermutationRate";
    public static final String HYPERMUTATION_INDICATORS = "hypermutationIndicators";
    public static final String UNLINKED_RATES = "unlinkedRates";

    public HypermutantErrorModel(HypermutantAlignment hypermutantAlignment, Parameter hypermutationRateParameter, Parameter hypermuationIndicatorParameter, boolean unlinkedRates) {
        super(HYPERMUTANT_ERROR_MODEL, null, null);

        this.hypermutantAlignment = hypermutantAlignment;

        this.unlinkedRates = unlinkedRates;

        this.hypermutationRateParameter = hypermutationRateParameter;
        addVariable(this.hypermutationRateParameter);

        this.hypermutationIndicatorParameter = hypermuationIndicatorParameter;

        addVariable(this.hypermutationIndicatorParameter);

        addStatistic(new TaxonHypermutatedStatistic());
        addStatistic(new TaxonHypermutationRateStatistic());
        addStatistic(new HypermutatedProportionStatistic());
    }

    protected void taxaChanged() {
        if (hypermutationIndicatorParameter.getDimension() <= 1) {
            this.hypermutationIndicatorParameter.setDimension(tree.getExternalNodeCount());
        }
        if (unlinkedRates && hypermutationRateParameter.getDimension() <= 1) {
            this.hypermutationRateParameter.setDimension(tree.getExternalNodeCount());
        }
    }

    @Override
    public Type getModelType() {
        return Type.PARTIALS;
    }

    @Override
    public void getTipStates(int nodeIndex, int[] tipStates) {
        throw new IllegalArgumentException("This model emits only tip partials");
    }

    @Override
    public void getTipPartials(int nodeIndex, double[] partials) {
        int[] states = this.states[nodeIndex];
        boolean isHypermutated = hypermutationIndicatorParameter.getParameterValue(nodeIndex) > 0.0;

        double rate = (unlinkedRates ? hypermutationRateParameter.getParameterValue(nodeIndex) : hypermutationRateParameter.getParameterValue(0));

        int k = 0;
        for (int j = 0; j < patternCount; j++) {

            switch (states[j]) {
                case Nucleotides.A_STATE: // is an A
                    partials[k] = 1.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.C_STATE: // is an C
                    partials[k] = 0.0;
                    partials[k + 1] = 1.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.G_STATE: // is an G
                    partials[k] = 0.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 1.0;
                    partials[k + 3] = 0.0;
                    break;
                case Nucleotides.UT_STATE: // is an T
                    partials[k] = 0.0;
                    partials[k + 1] = 0.0;
                    partials[k + 2] = 0.0;
                    partials[k + 3] = 1.0;
                    break;
                case Nucleotides.R_STATE: // is an A in a APOBEC context
                    if (isHypermutated) {
                        partials[k] = 1.0 - rate;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = rate;
                        partials[k + 3] = 0.0;
                    } else {
                        partials[k] = 1.0;
                        partials[k + 1] = 0.0;
                        partials[k + 2] = 0.0;
                        partials[k + 3] = 0.0;
                    }

                    break;
                default: // is an ambiguity
                    partials[k] = 1.0;
                    partials[k + 1] = 1.0;
                    partials[k + 2] = 1.0;
                    partials[k + 3] = 1.0;
            }

            k += stateCount;
        }

    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == hypermutationIndicatorParameter) {
            fireModelChanged(tree.getTaxon(index));
        } else if (variable == hypermutationRateParameter) {
            if (!unlinkedRates) {
                fireModelChanged();
            } else if (hypermutationIndicatorParameter.getValue(index)  > 0.5) {
                // only fire an update if the indicator is on....
                fireModelChanged(tree.getTaxon(index));
            }
        } else {
            throw new RuntimeException("Unknown parameter has changed in HypermutantErrorModel.handleVariableChangedEvent");
        }

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return HYPERMUTANT_ERROR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean unlinkedRates = false;
            if (xo.hasAttribute(UNLINKED_RATES)) {
                unlinkedRates = xo.getBooleanAttribute(UNLINKED_RATES);
            }

            HypermutantAlignment hypermutantAlignment = (HypermutantAlignment)xo.getChild(HypermutantAlignment.class);

            Parameter hypermutationRateParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_RATE)) {
                hypermutationRateParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_RATE);
            }

            Parameter hypermuationIndicatorParameter = null;
            if (xo.hasChildNamed(HYPERMUTATION_INDICATORS)) {
                hypermuationIndicatorParameter = (Parameter)xo.getElementFirstChild(HYPERMUTATION_INDICATORS);
            }

            HypermutantErrorModel errorModel =  new HypermutantErrorModel(hypermutantAlignment, hypermutationRateParameter, hypermuationIndicatorParameter, unlinkedRates);

            Logger.getLogger("dr.evomodel").info("Using APOBEC error model");

            return errorModel;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns a model that allows for APOBEC-type RNA editing.";
        }

        public Class getReturnType() { return HypermutantErrorModel.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                AttributeRule.newBooleanRule(UNLINKED_RATES, true),
                new ElementRule(HypermutantAlignment.class),
                new ElementRule(HYPERMUTATION_RATE, Parameter.class, "The hypermutation rate per target site per sequence"),
                new ElementRule(HYPERMUTATION_INDICATORS, Parameter.class, "A binary indicator of whether the sequence is hypermutated"),
        };
    };

    public class TaxonHypermutatedStatistic extends Statistic.Abstract {

        public TaxonHypermutatedStatistic() {
            super("isHypermutated");
        }

        public int getDimension() {
            return hypermutationIndicatorParameter.getDimension();
        }

        public String getDimensionName(int dim) {
            return taxonMap.get(dim);
        }

        public double getStatisticValue(int dim) {
            return hypermutationIndicatorParameter.getParameterValue(dim);
        }

    }

    public class TaxonHypermutationRateStatistic extends Statistic.Abstract {

        public TaxonHypermutationRateStatistic() {
            super("hypermutationRate");
        }

        public int getDimension() {
            return hypermutationRateParameter.getDimension();
        }

        public String getDimensionName(int dim) {
            return taxonMap.get(dim) + ".rate";
        }

        public double getStatisticValue(int dim) {
            return hypermutationRateParameter.getParameterValue(dim) * hypermutationIndicatorParameter.getParameterValue(dim);
        }

    }

    public class HypermutatedProportionStatistic extends Statistic.Abstract {

        public HypermutatedProportionStatistic() {
            super("proportionHypermutated");
        }

        public int getDimension() {
            return 1;
        }

        public String getDimensionName(int dim) {
            return "P(hypermutated)";
        }

        public double getStatisticValue(int dim) {
            if (mutatedContextCounts == null) {
                mutatedContextCounts = new HashMap<Integer, Integer>();
                unmutatedContextCounts = new HashMap<Integer, Integer>();

                for (int index : taxonMap.keySet()) {
                    String name = taxonMap.get(index);
                    int i = hypermutantAlignment.getTaxonIndex(name);
                    mutatedContextCounts.put(index, hypermutantAlignment.getMutatedContextCounts()[i]);
                    unmutatedContextCounts.put(index, hypermutantAlignment.getUnmutatedContextCounts()[i]);
                }
            }

            double mutatedCount = 0;
            double totalCount = 0;
            for (int i = 0; i < hypermutationIndicatorParameter.getDimension(); i++) {
                if (hypermutationIndicatorParameter.getParameterValue(i) > 0.5) {
                    mutatedCount += mutatedContextCounts.get(i);
                    totalCount += mutatedCount + unmutatedContextCounts.get(i);

                }
            }

            double r = hypermutationRateParameter.getParameterValue(0);
            return (r * mutatedCount) / totalCount;
        }

    }

    Map<Integer, Integer> mutatedContextCounts = null;
    Map<Integer, Integer> unmutatedContextCounts = null;

    private final HypermutantAlignment hypermutantAlignment;
    private final Parameter hypermutationRateParameter;
    private final Parameter hypermutationIndicatorParameter;
    private final boolean unlinkedRates;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Sequence error model";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(new Citation(
                        new Author[]{
                                new Author("A", "Rambaut"),
                                new Author("SYW", "Ho"),
                                new Author("AJ", "Drummond"),
                                new Author("B", "Shapiro"),
                        },
                        "Accommodating the effect of ancient DNA damage on inferences of demographic histories",
                        2008,
                        "Mol Biol Evol",
                        26,
                        245, 248,
                        "10.1093/molbev/msn256"
                ),
                new Citation(
                        new Author[]{
                                new Author("J", "Felsenstein"),
                        },
                        "Inferring Phylogenies",
                        2004,
                        "Sinauer Associates",
                        ""
                ));
    }
}