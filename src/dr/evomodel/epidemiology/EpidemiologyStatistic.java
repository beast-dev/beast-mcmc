/*
 * EpidemiologyStatistic.java
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

package dr.evomodel.epidemiology;

import dr.evolution.util.Units;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;

/**
 * A statistic that converts growth rate to a doubling time or an estimate of R0
 *
 * Approach based on the Grassly and Fraser modification of
 * Wallinga and Lipsitch (2007, Proc Roy Soc B, 274:599–604)
 *
 * R_{0}=\left(1+\frac{r}{b}\right)^{a}
 *
 * where a and b are the parameters of the gamma distribution
 * ( a = m^{2}/s^{2}  and  b = m/s^{2} where m and s are the
 * mean and standard deviation of the distribution respectively).
 *
 * @author Andrew Rambaut
 */
public class EpidemiologyStatistic extends Statistic.Abstract implements Units {
    public enum StatisticType {
        DOUBLING_TIME("doublingTime"),
        GROWTH_RATE("growthRate"),
        R0("R0");

        StatisticType(String name) {
            this.name = name;
        }

        @Override
        public String toString() {
            return name;
        }

        String name;
    }

    /**
     * Constructor for a doubling time statistic
     * @param name
     * @param growthRate
     * @param doublingTime
     * @param units
     */
    public EpidemiologyStatistic(String name, StatisticType statisticType, Parameter growthRate, Parameter doublingTime, Units.Type units) {
        this(name, statisticType, growthRate, doublingTime, units, 0, 0, null);
    }

    /**
     * Constructor for an R0 statistic with a parameterised serial interval
     * @param name
     * @param growthRate
     * @param doublingTime
     * @param units
     * @param serialIntervalMean
     * @param serialIntervalStdev
     */
    public EpidemiologyStatistic(String name, Parameter growthRate, Parameter doublingTime,  Units.Type units,
                                 double serialIntervalMean, double serialIntervalStdev) {
        this(name, StatisticType.R0, growthRate, doublingTime, units, serialIntervalMean, serialIntervalStdev, null);
    }

    private EpidemiologyStatistic(String name, StatisticType statisticType, Parameter growthRate,
                                  Parameter doublingTime,  Units.Type units,
                                  double serialIntervalMean, double serialIntervalStdev, double[] serialIntervalPDF) {
        super(name);

        this.statisticType = statisticType;

        this.growthRate = growthRate;
        this.doublingTime = doublingTime;
        if (units != Type.YEARS && units != Type.MONTHS && units != Type.WEEKS && units != Type.DAYS) {
            throw new UnsupportedOperationException("Growth rate time units should be years, months, weeks, or days");
        }
        this.units = units;
        this.a = (serialIntervalMean * serialIntervalMean) / (serialIntervalStdev * serialIntervalStdev);
        this.b = serialIntervalMean / (serialIntervalStdev * serialIntervalStdev);

        if (serialIntervalPDF != null) {
            // normalise
            this.w = new double[serialIntervalPDF.length];
            double sum = 0.0;
            for (double value : serialIntervalPDF) {
                sum += value;
            }
            for (int i = 0; i < w.length; i++) {
                this.w[i] = serialIntervalPDF[i] / sum;
            }
        } else {
            this.w = null;
        }
    }

    @Override
    public String getDimensionName(int dim) {
        return super.getDimensionName(dim);
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {

        double gr; // growth rate in original units
        if (growthRate != null) {
            gr = growthRate.getParameterValue(0);
        } else {
            // compute the growth rate from the doubling time
            double dt = doublingTime.getParameterValue(0);
            gr = Math.log(2) / dt;
        }

        double r;
        switch (units) {
            // convert the growth rate to a per day rate
            case DAYS:
                r = gr;
                break;
            case WEEKS:
                r = gr / 7;
                break;
            case MONTHS:
                r = gr / 30; // approximate
                break;
            case YEARS:
                r = gr / 365;
                break;
            default:
                throw new UnsupportedOperationException("Unsupported time unit");
        }

        if (statisticType == StatisticType.DOUBLING_TIME) {
            // doubling time in days
            double dt = Math.log(2) / r;
            // clip to zero?
//            return Math.max(0.0, dt);
            return dt;
        } else if (statisticType == StatisticType.GROWTH_RATE) {
            // growth rate in original units
            return gr;
        } else { // R0
            double R0 = 1.0;
            if (w == null) {
                R0 = Math.pow((1.0 + (r / b)), a);
            } else {
                throw new UnsupportedOperationException("R0 estimate based on empirical serial interval not implemented yet");
            }
            return R0;
        }
    }

    private final StatisticType statisticType;
    private final Parameter growthRate;
    private final Parameter doublingTime;
    private final Type units;
    private final double[] w;
    private final double a, b;

    @Override
    public Type getUnits() {
        return Type.DAYS;
    }

    @Override
    public void setUnits(Type units) {
        throw new UnsupportedOperationException("Cannot set units");
    }
}
