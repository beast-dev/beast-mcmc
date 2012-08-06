/*
 * DiscreteTraitsComponentOptions.java
 *
 * Copyright (c) 2002-2011 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.components.discrete;

import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.evolution.datatype.GeneralDataType;
import dr.inference.operators.RateBitExchangeOperator;

import java.util.List;
import java.util.Set;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class DiscreteTraitsComponentOptions implements ComponentOptions {

    private final BeautiOptions options;

    public DiscreteTraitsComponentOptions(final BeautiOptions options) {
        this.options = options;
    }

    public void createParameters(final ModelOptions modelOptions) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(GeneralDataType.INSTANCE)) {
            String prefix = partitionData.getName() + ".";

            if (!modelOptions.parameterExists(prefix + "frequencies")) {

                modelOptions.createZeroOneParameterUniformPrior(prefix + "frequencies", "discrete state frequencies", 0.25);
                modelOptions.createCachedGammaPrior(prefix + "rates", "discrete trait instantaneous transition rates",
                        PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0, 1.0, false);

                // BSSVS
                modelOptions.createParameter(prefix + "indicators", "a vector of bits indicating non-zero rates for BSSVS", 1.0);

                // = strick clock TODO trait.mu belongs Clock Model?
//        modelOptions.createParameterExponentialPrior(prefix + "mu", "discrete trait CTMC rate parameter",
//                PriorScaleType.SUBSTITUTION_PARAMETER_SCALE, 0.1, 1.0, 0.0);

                // Poisson Prior on non zero ratesBSSVS
                modelOptions.createDiscreteStatistic(prefix + "nonZeroRates", "the number of non-zero rates for BSSVS");

                modelOptions.createScaleOperator(prefix + "frequencies", 0.75, 1.0);
                modelOptions.createOperator(prefix + "rates", OperatorType.SCALE_INDEPENDENTLY, 0.75, 1.0);
                modelOptions.createOperator(prefix + "indicators", OperatorType.BITFLIP, -1.0, 1.0);
//                modelOptions.createScaleOperator(prefix + "mu", demoTuning, 10);

                //bit Flip on clock.rate in PartitionClockModelSubstModelLink
//                modelOptions.createBitFlipInSubstitutionModelOperator(OperatorType.BITFIP_IN_SUBST.toString() + "mu", prefix + "mu",
//                        "bit Flip In Substitution Model Operator on trait.mu", getParameter("trait.mu"), this, demoTuning, 30);
                modelOptions.createOperatorUsing2Parameters(RateBitExchangeOperator.OPERATOR_NAME,
                        "(indicators, rates)",
                        "rateBitExchangeOperator (If both BSSVS and asymmetric subst selected)",
                        prefix + "indicators", prefix + "rates", OperatorType.RATE_BIT_EXCHANGE, -1.0, 6.0);

            }
        }
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(GeneralDataType.INSTANCE)) {
            String prefix = partitionData.getName() + ".";

            if (partitionData.getPartitionSubstitutionModel().isActivateBSSVS()) {
                modelOptions.getParameter(prefix + "indicators");
                Parameter nonZeroRates = modelOptions.getParameter(prefix + "nonZeroRates");

                Set<String> states = partitionData.getPartitionSubstitutionModel().getDiscreteStateSet();
                int K = states.size();
                if (partitionData.getPartitionSubstitutionModel().getDiscreteSubstType() == DiscreteSubstModelType.SYM_SUBST) {
                    nonZeroRates.mean = Math.log(2); // mean = 0.693 and offset = K-1
                    nonZeroRates.offset = K - 1;
                } else if (partitionData.getPartitionSubstitutionModel().getDiscreteSubstType() == DiscreteSubstModelType.ASYM_SUBST) {
                    nonZeroRates.mean = K - 1; // mean = K-1 and offset = 0
                    nonZeroRates.offset = 0.0;
                } else {
                    throw new IllegalArgumentException("unknown discrete substitution type");
                }

                params.add(nonZeroRates);
            }
            params.add(modelOptions.getParameter(prefix + "frequencies"));
            params.add(modelOptions.getParameter(prefix + "rates"));
//            System.err.println("Selected OK");
//               params.add(getParameter(prefix + "mu"));
        }

    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        for (AbstractPartitionData partitionData : options.getDataPartitions(GeneralDataType.INSTANCE)) {
            String prefix = partitionData.getName() + ".";

//            ops.add(modelOptions.getOperator(prefix + "frequencies")); // Usually fixed
            ops.add(modelOptions.getOperator(prefix + "rates"));

            if (partitionData.getPartitionSubstitutionModel().isActivateBSSVS()) {
                ops.add(modelOptions.getOperator(prefix + "indicators"));
            }
        }
    }


}