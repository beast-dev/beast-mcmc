/*
 * ClockModelGroup.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.FixRateType;
import dr.evolution.datatype.DataType;

import java.io.Serializable;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ClockModelGroup implements Serializable {

    private static final long serialVersionUID = -3034174176050520635L;
    private String name;
    private boolean fixMean = false;
    private double fixMeanRate = 1.0;
    private FixRateType rateTypeOption = FixRateType.RELATIVE_TO;

//    public ClockModelGroup() { }

    public ClockModelGroup(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public boolean isFixMean() {
        return fixMean;
    }

    public void setFixMean(boolean fixMean) {
        this.fixMean = fixMean;
    }

    public FixRateType getRateTypeOption() {
        return rateTypeOption;
    }

    public void setRateTypeOption(FixRateType rateTypeOption) {
        this.rateTypeOption = rateTypeOption;
        setFixMean(rateTypeOption == FixRateType.FIX_MEAN);
    }

    public double getFixMeanRate() {
        return fixMeanRate;
    }

    public void setFixMeanRate(double fixMeanRate, BeautiOptions options) {
        this.fixMeanRate = fixMeanRate;
        for (PartitionClockModel model : options.getPartitionClockModels(this)) {
            model.setRate(fixMeanRate, false);
        }
    }

    public boolean contain(DataType dataType, BeautiOptions options) {
        for (AbstractPartitionData pd : options.getDataPartitions(this)) {
           if (pd.getDataType().getType() == dataType.getType()) {
                return true;
            }
        }
        return false;
    }

}
