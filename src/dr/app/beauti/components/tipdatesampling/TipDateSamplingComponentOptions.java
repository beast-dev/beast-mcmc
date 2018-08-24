/*
 * TipDateSamplingComponentOptions.java
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

package dr.app.beauti.components.tipdatesampling;

import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.types.TipDateSamplingType;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class TipDateSamplingComponentOptions implements ComponentOptions {

    private final BeautiOptions options;

    public TipDateSamplingComponentOptions(final BeautiOptions options) {
        this.options = options;
    }

    public void createParameters(final ModelOptions modelOptions) {
        modelOptions.createNonNegativeParameterInfinitePrior("treeModel.tipDates", "date of specified tips",
                PriorScaleType.TIME_SCALE, 1.0);

        modelOptions.createScaleOperator("treeModel.tipDates", modelOptions.demoTuning, 3.0);
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY ||
                tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
            TaxonList taxa = getTaxonSet();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                Taxon taxon = taxa.getTaxon(i);

                Parameter parameter = getTipDateParameter(taxon);

                if (tipDateSamplingType != TipDateSamplingType.SAMPLE_PRECISION) {
                    // no explicit priors when sampling from precisions - the random walk operators
                    // sample uniformly from within the bounds defined by the precision.
                    params.add(parameter);
                }
            }

        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
            params.add(modelOptions.getParameter("treeModel.tipDates"));
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY ||
                tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
            TaxonList taxa = getTaxonSet();

            String description = "Random walk for the age of this tip";

            //OperatorType type = OperatorType.RANDOM_WALK_ABSORBING;
            OperatorType type = OperatorType.UNIFORM;

            if (tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
                description = "Uniform sample from precision of age of this tip";
                //type = OperatorType.UNIFORM;
            }

            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                Taxon taxon = taxa.getTaxon(i);

                Operator operator = tipDateOperators.get(taxon);

                if (operator == null) {
                    Parameter parameter = getTipDateParameter(taxon);

//                    operator = new Operator("age(" + taxon.getId() + ")", "", parameter, OperatorType.SCALE, 0.75, 1.0);
                    operator = new Operator.Builder("age(" + taxon.getId() + ")", description, parameter, type, 1.0, 1.0).build();
                    if (tipDateSamplingType == TipDateSamplingType.SAMPLE_INDIVIDUALLY) {
                        operator.setWeight(2.0);
                    }
                    tipDateOperators.put(taxon, operator);
                }
                ops.add(operator);
            }
        } else if (tipDateSamplingType == TipDateSamplingType.SAMPLE_JOINT) {
            ops.add(modelOptions.getOperator("treeModel.tipDates"));
        }
    }

    private Parameter getTipDateParameter(Taxon taxon) {
        Parameter parameter = tipDateParameters.get(taxon);

        if (parameter == null) {
            double height = 0.0;
            if (taxon.getAttribute("height") != null) {
                height = (Double)taxon.getAttribute("height");
            }

            parameter = new Parameter.Builder("age(" + taxon.getId() + ")", "sampled age of taxon, " + taxon.getId())
                    .scaleType(PriorScaleType.TIME_SCALE).prior(PriorType.UNIFORM_PRIOR).initial(height).isNonNegative(true).build();
            parameter.setPriorEdited(true);
            tipDateParameters.put(taxon, parameter);
        }

        return parameter;
    }

    public TaxonList getTaxonSet() {
        TaxonList taxa = options.taxonList;

        if (tipDateSamplingTaxonSet != null) {
            taxa = tipDateSamplingTaxonSet;
        }

        if (tipDateSamplingType == TipDateSamplingType.SAMPLE_PRECISION) {
            Taxa precisionTaxonList = new Taxa();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                Taxon taxon = taxa.getTaxon(i);
                Date date = taxon.getDate();
                if (date.getUncertainty() > 0.0) {
                    precisionTaxonList.addTaxon(taxon);
                }

            }
            taxa = precisionTaxonList;
        }

        return taxa;
    }

    public TipDateSamplingType tipDateSamplingType = TipDateSamplingType.NO_SAMPLING;
    public TaxonList tipDateSamplingTaxonSet = null;

    private Map<Taxon, Parameter> tipDateParameters = new HashMap<Taxon, Parameter>();
    private Map<Taxon, Operator> tipDateOperators = new HashMap<Taxon, Operator>();
}