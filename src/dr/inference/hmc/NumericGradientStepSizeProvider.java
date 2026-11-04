/*
 * NumericGradientStepSizeProvider.java
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

package dr.inference.hmc;

import dr.math.MachineAccuracy;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public interface NumericGradientStepSizeProvider {
    default double getNumericGradientStepSize() {
        return StepSizeLevel.MEDIUM.getStepSizeRatio();
    }

    default void setNumericGradientStepSize(double ratio) {

    }

    String NUMERIC_STEP_SIZE = "numericStepSize";
    static double parseStepSizeRatio(XMLObject xo)  throws XMLParseException {

        if (xo.hasAttribute(NUMERIC_STEP_SIZE)) {
            String stepSizeLevelString = xo.getStringAttribute(NUMERIC_STEP_SIZE);
            for (StepSizeLevel level: StepSizeLevel.values()) {
                if (stepSizeLevelString.equalsIgnoreCase(level.getName())) {
                    return level.getStepSizeRatio();
                }
            }
        }
        return xo.getDoubleAttribute(NUMERIC_STEP_SIZE, StepSizeLevel.MEDIUM.getStepSizeRatio());
    }

    enum StepSizeLevel {
        TINY("tiny", MachineAccuracy.EPSILON),
        SMALL("small", MachineAccuracy.SQRT_EPSILON),
        MEDIUM("medium", Math.pow(MachineAccuracy.EPSILON, 0.333)),
        LARGE("large", MachineAccuracy.SQRT_SQRT_EPSILON);

        StepSizeLevel(String name, double stepSizeRatio) {
            this.name = name;
            this.stepSizeRatio = stepSizeRatio;
        }

        public String getName() {return name;}
        public double getStepSizeRatio() {return stepSizeRatio;}

        private String name;
        private double stepSizeRatio;
    }
}
