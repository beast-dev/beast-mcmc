/*
 * SequenceErrorModelComponentOptions.java
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

package dr.app.beauti.components.sequenceerror;

import dr.app.beauti.options.*;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.SequenceErrorType;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class SequenceErrorModelComponentOptions implements ComponentOptions {
    static public final String ERROR_MODEL = "errorModel";
    static public final String AGE_RATE = "ageRelatedErrorRate";
    static public final String BASE_RATE = "baseErrorRate";

    static public final String AGE_RATE_PARAMETER = ERROR_MODEL + "." + AGE_RATE;
    static public final String BASE_RATE_PARAMETER = ERROR_MODEL + "." + BASE_RATE;

    static public final String HYPERMUTION_RATE_PARAMETER = "hypermutation.rate";
    static public final String HYPERMUTANT_INDICATOR_PARAMETER = "hypermutant.indicator";
    static public final String HYPERMUTANT_COUNT_STATISTIC = "hypermutation.count";

    SequenceErrorModelComponentOptions() {
    }

    public void createParameters(final ModelOptions modelOptions) {
        for (AbstractPartitionData partition : sequenceErrorTypeMap.keySet()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";
            modelOptions.createNonNegativeParameterInfinitePrior(prefix + AGE_RATE_PARAMETER,"age dependent sequence error rate",
                    PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0E-8);
            modelOptions.createZeroOneParameterUniformPrior(prefix + BASE_RATE_PARAMETER,"base sequence error rate", 1.0E-8);

            modelOptions.createZeroOneParameterUniformPrior(prefix + HYPERMUTION_RATE_PARAMETER,"APOBEC editing rate per context", 1.0E-8);
            modelOptions.createParameter(prefix + HYPERMUTANT_INDICATOR_PARAMETER, "indicator parameter reflecting which sequences are hypermutated", 0.0);

            modelOptions.createDiscreteStatistic(prefix + HYPERMUTANT_COUNT_STATISTIC, "count of the number of hypermutated sequences");

            modelOptions.createScaleOperator(prefix + AGE_RATE_PARAMETER, modelOptions.demoTuning, 3.0);
            modelOptions.createOperator(prefix + BASE_RATE_PARAMETER, OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);

            modelOptions.createOperator(prefix + HYPERMUTION_RATE_PARAMETER, OperatorType.RANDOM_WALK_REFLECTING, 0.05, 3.0);
            modelOptions.createOperator(prefix + HYPERMUTANT_INDICATOR_PARAMETER, OperatorType.BITFLIP, -1.0, 10);
        }
    }

    public void selectParameters(final ModelOptions modelOptions, final List<Parameter> params) {
        for (AbstractPartitionData partition : sequenceErrorTypeMap.keySet()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";
            if (isHypermutation(partition)) {
                params.add(modelOptions.getParameter(prefix + HYPERMUTION_RATE_PARAMETER));
                params.add(modelOptions.getParameter(prefix + HYPERMUTANT_INDICATOR_PARAMETER));
                params.add(modelOptions.getParameter(prefix + HYPERMUTANT_COUNT_STATISTIC));
            }
            if (hasAgeDependentRate(partition)) {
                params.add(modelOptions.getParameter(prefix + AGE_RATE_PARAMETER));
            }
            if (hasBaseRate(partition)) {
                params.add(modelOptions.getParameter(prefix + BASE_RATE_PARAMETER));
            }
        }
    }

    public void selectStatistics(final ModelOptions modelOptions, final List<Parameter> stats) {
        // no statistics required
    }

    public void selectOperators(final ModelOptions modelOptions, final List<Operator> ops) {
        for (AbstractPartitionData partition : sequenceErrorTypeMap.keySet()) {
            String prefix = partition.getPrefix();//partition.getName() + ".";
            if (isHypermutation(partition)) {
                ops.add(modelOptions.getOperator(prefix + HYPERMUTION_RATE_PARAMETER));
                ops.add(modelOptions.getOperator(prefix + HYPERMUTANT_INDICATOR_PARAMETER));
            }
            if (hasAgeDependentRate(partition)) {
                ops.add(modelOptions.getOperator(prefix + AGE_RATE_PARAMETER));
            }
            if (hasBaseRate(partition)) {
                ops.add(modelOptions.getOperator(prefix + BASE_RATE_PARAMETER));
            }
        }
    }

    public boolean usingSequenceErrorModel() {
        for (AbstractPartitionData partition : sequenceErrorTypeMap.keySet()) {
            if (sequenceErrorTypeMap.get(partition) != SequenceErrorType.NO_ERROR) {
                return true;
            }
        }
        return false;
    }

    public boolean usingSequenceErrorModel(AbstractPartitionData partition) {
       return (getSequenceErrorType(partition) != SequenceErrorType.NO_ERROR);
    }

    public boolean hasAgeDependentRate(final AbstractPartitionData partition) {
        SequenceErrorType errorModelType = getSequenceErrorType(partition);
        return (errorModelType == SequenceErrorType.AGE_ALL) || (errorModelType == SequenceErrorType.AGE_TRANSITIONS);
    }

    public boolean hasBaseRate(final AbstractPartitionData partition) {
        SequenceErrorType errorModelType = getSequenceErrorType(partition);
        return (errorModelType == SequenceErrorType.BASE_ALL) || (errorModelType == SequenceErrorType.BASE_TRANSITIONS);
    }

    public boolean isHypermutation(final AbstractPartitionData partition) {
        SequenceErrorType errorModelType = getSequenceErrorType(partition);
        return (errorModelType == SequenceErrorType.HYPERMUTATION_ALL) ||
                (errorModelType == SequenceErrorType.HYPERMUTATION_BOTH) ||
                (errorModelType == SequenceErrorType.HYPERMUTATION_HA3G) ||
                (errorModelType == SequenceErrorType.HYPERMUTATION_HA3F);
    }

    public SequenceErrorType getSequenceErrorType(final AbstractPartitionData partition) {
        SequenceErrorType type = sequenceErrorTypeMap.get(partition);
        if (type == null) {
            type = SequenceErrorType.NO_ERROR;
            sequenceErrorTypeMap.put(partition, type);
        }
        return type;
    }

    public void setSequenceErrorType(final AbstractPartitionData partition, SequenceErrorType sequenceErrorType) {
        sequenceErrorTypeMap.put(partition, sequenceErrorType);
    }

    private final Map<AbstractPartitionData, SequenceErrorType> sequenceErrorTypeMap = new HashMap<AbstractPartitionData, SequenceErrorType>();
}