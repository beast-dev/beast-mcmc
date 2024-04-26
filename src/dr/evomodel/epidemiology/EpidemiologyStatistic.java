/*
 * ExternalLengthStatistic.java
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
     * @param growthRateUnits
     */
    public EpidemiologyStatistic(String name, Parameter growthRate, Units.Type growthRateUnits) {
        this(name, StatisticType.DOUBLING_TIME, growthRate, growthRateUnits, 0, 0, null);
    }

    /**
     * Constructor for an R0 statistic with a parameterised serial interval
     * @param name
     * @param growthRate
     * @param growthRateUnits
     * @param serialIntervalMean
     * @param serialIntervalStdev
     */
    public EpidemiologyStatistic(String name, Parameter growthRate, Units.Type growthRateUnits,
                                 double serialIntervalMean, double serialIntervalStdev) {
        this(name, StatisticType.R0, growthRate, growthRateUnits, serialIntervalMean, serialIntervalStdev, null);
    }

    /**
     * Constructor for an R0 statistic with a discretised serial interval
     * @param name
     * @param growthRate
     * @param growthRateUnits
     * @param serialIntervalPDF
     */
    public EpidemiologyStatistic(String name, Parameter growthRate, Units.Type growthRateUnits,
                                 double[] serialIntervalPDF) {
        this(name, StatisticType.R0, growthRate, growthRateUnits, 0, 0, serialIntervalPDF);
    }

    private EpidemiologyStatistic(String name, StatisticType statisticType, Parameter growthRate, Units.Type growthRateUnits,
                                  double serialIntervalMean, double serialIntervalStdev, double[] serialIntervalPDF) {
        super(name);

        this.statisticType = statisticType;

        this.growthRate = growthRate;
        if (growthRateUnits != Type.YEARS && growthRateUnits != Type.MONTHS && growthRateUnits != Type.WEEKS && growthRateUnits != Type.DAYS) {
            throw new UnsupportedOperationException("Growth rate time units should be years, months, weeks, or days");
        }
        this.growthRateUnits = growthRateUnits;
        this.a = (serialIntervalMean * serialIntervalMean);
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
        return super.getDimensionName(dim) + "." + statisticType.toString();
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return the height of the MRCA node.
     */
    public double getStatisticValue(int dim) {
        double r;
        switch (growthRateUnits) {
            // convert the growth rate to a per day rate
            case DAYS:
                r = growthRate.getParameterValue(0);
                break;
            case WEEKS:
                r = growthRate.getParameterValue(0) / 7;
                break;
            case MONTHS:
                r = growthRate.getParameterValue(0) / 30; // approximate
                break;
            case YEARS:
                r = growthRate.getParameterValue(0) / 365;
                break;
            default:
                throw new UnsupportedOperationException("Unsupported time unit");
        }

        if (statisticType == StatisticType.DOUBLING_TIME) {
            // doubling time in days
            return Math.log(2) / r;
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
    private final Type growthRateUnits;
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
