/*
 * DolloComponentOptions.java
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

package dr.app.beauti.components.dollo;

import dr.app.beauti.options.*;
import dr.app.beauti.types.PriorScaleType;

import java.util.List;

/**
 * @author Marc Suchard
 * @version $Id$
 */

public class DolloComponentOptions implements ComponentOptions {

    public static final String DEATH_RATE = "death.rate";
    public static final String DATA_NAME = "binaryDolloDataType";
    public static final String MODEL_NAME = "binaryDolloSubstModel";

	final private BeautiOptions options;

	public DolloComponentOptions(final BeautiOptions options) {
		this.options = options;
	}

	public void createParameters(ModelOptions modelOptions) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getPartitionSubstitutionModel().isDolloModel()) {
                String prefix = partition.getName() + ".";
                modelOptions.createParameterExponentialPrior(prefix + DEATH_RATE, "Stochastic Dollo death rate",
                        PriorScaleType.SUBSTITUTION_RATE_SCALE, 1.0E-5, 1.0E-4, 0.0);
                modelOptions.createScaleOperator(prefix + DEATH_RATE, modelOptions.demoTuning, 1.0);
            }
        }
	}

	public void selectOperators(ModelOptions modelOptions, List<Operator> ops) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model.isDolloModel()) {
                String prefix = partition.getName() + ".";
                ops.add(modelOptions.getOperator(prefix + DEATH_RATE));
            }
        }
	}

	public void selectParameters(ModelOptions modelOptions,
			List<Parameter> params) {
        for (AbstractPartitionData partition : options.dataPartitions) {
            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();
            if (model.isDolloModel()) {
                String prefix = partition.getName() + ".";
                params.add(modelOptions.getParameter(prefix + DEATH_RATE));
                break;
            }
        }
	}

	public void selectStatistics(ModelOptions modelOptions,
			List<Parameter> stats) {
		// Do nothing
	}

	public BeautiOptions getOptions() {
		return options;
	}

    public void setActive(boolean active) {
        this.active = active;
    }

    public boolean isActive() {
        return active;
    }

    private boolean active = false;
}