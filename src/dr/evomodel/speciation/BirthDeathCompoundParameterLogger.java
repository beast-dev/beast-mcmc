/*
 * BirthDeathCompoundParameterLogger.java
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

package dr.evomodel.speciation;

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Statistic;
import dr.math.UnivariateMinimum;

// TODO use "implements Loggable" instead of "extends Statistic.Abstract"
// TODO hard-code names
public class BirthDeathCompoundParameterLogger extends Statistic.Abstract {

    public enum BDPCompoundParameterType {
        EFFECTIVE_REPRODUCTIVE_NUMBER("effectiveReproductiveNumber") {
            public double getCompoundParameterForType(double birthRate, double deathRate, double samplingRate, double treatmentProbability, double samplingProbability) {
                return birthRate / (samplingRate * treatmentProbability + deathRate);
            }
        };

        BDPCompoundParameterType(String name) {
            this.name = name;
//            this.label = label;
        }

        public String getName() {
            return name;
        }

//        public String getLabel() {
//            return label;
//        }

        // TODO may need to be
        public abstract double getCompoundParameterForType(double birthRate, double deathRate, double samplingRate, double treatmentProbability, double samplingProbability);

        private String name;
//        private String label;
    }

    public BirthDeathCompoundParameterLogger(NewBirthDeathSerialSamplingModel bdss, BDPCompoundParameterType type) {
        this.bdss = bdss;
        this.type = type;
        this.dim = bdss.getDeathRateParameter().getDimension();
    }

    private double getCompoundParameter(int i) {
        double birth = bdss.getBirthRateParameter().getParameterValue(i);
        double death = bdss.getDeathRateParameter().getParameterValue(i);
        double sampling = bdss.getSamplingRateParameter().getParameterValue(i);
        double treatment = bdss.getTreatmentProbabilityParameter().getParameterValue(i);
        double prob = bdss.getSamplingProbabilityParameter().getParameterValue(i);

        return type.getCompoundParameterForType(birth, death, sampling, treatment, prob);
    }

    @Override
    public int getDimension() {
        return dim;
    }

    @Override
    public double getStatisticValue(int dim) {
        return getCompoundParameter(dim);
    }

    private final NewBirthDeathSerialSamplingModel bdss;
    private final int dim;
    private final BDPCompoundParameterType type;

}
