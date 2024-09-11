/*
 * GradientErrorLogger.java
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

import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andy Magee
 * @author Andrew Holbrook
 */
public class GradientErrorLogger implements Loggable, Reportable {

    private final GradientWrtParameterProvider source;
    private final List<Statistic> statistics;

    public GradientErrorLogger(GradientWrtParameterProvider source,
                               List<Statistic> statistics) {
        this.source = source;
        this.statistics = statistics;
    }

    private double getStatisticValue(Statistic statistic) {
        if (!gradientKnown) {
            gradient = source.getGradientLogDensity();
            reference = new GradientWrtParameterProvider.CheckGradientNumerically(source,
                    Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                    null, null
            ).getNumericalGradient();

            gradientKnown = true;
        }
        return statistic.getStatistic(gradient, reference);
    }

    private boolean gradientKnown = false;
    private double[] gradient;
    private double[] reference;

    @Override
    public LogColumn[] getColumns() {

        LogColumn[] columns = new LogColumn[statistics.size()];
        for (int i = 0; i < columns.length; ++i) {
            final Statistic statistic = statistics.get(i);
            final int index = i;
            final String name = source.getParameter().getId() + "." + statistic.getName();
            columns[i] = new NumberColumn(name) {
                @Override
                public double getDoubleValue() {
                    if (index == 0) {
                        gradientKnown = false;
                    }
                    return getStatisticValue(statistic);
                }
            };
        }

        return columns;
    }

    @Override
    public String getReport() {
        throw new RuntimeException("Not yet implemented");
    }

    public enum Statistic {
        MAX_ERROR_ABSOLUTE("maxErrorAbsolute") {
            @Override
            double getStatistic(double[] gradient, double[] reference) {
                return maxAbsDifference(gradient, reference, 1.0);
            }
        },
        MAX_ERROR_RELATIVE("maxErrorRelative") {
            @Override
            double getStatistic(double[] gradient, double[] reference) {
                double relative = Math.sqrt(
                        innerProduct(gradient, gradient) *
                        innerProduct(reference, reference));

                return maxAbsDifference(gradient, reference, relative);
            }
        },
        ANGLE("angle"){
            @Override
            double getStatistic(double[] gradient, double[] reference) {
                double innerProduct = innerProduct(gradient, reference);
                double gradientNorm = Math.sqrt(innerProduct(gradient, gradient));
                double referenceNorm = Math.sqrt(innerProduct(reference, reference));

                return Math.acos(innerProduct / (gradientNorm * referenceNorm));
            }
        };

        Statistic(String name) {
            this.name = name;
        }

        abstract double getStatistic(double[] gradient, double[] reference);

        public final String getName() { return name; }

        private final String name;

        public static Statistic parse(String name) {
            for (Statistic statistic : Statistic.values()) {
                if (statistic.name.equalsIgnoreCase(name)) {
                    return statistic;
                }
            }
            return null;
        }
    }

    private static double maxAbsDifference(double[] gradient, double[] reference, double relative) {
        assert  gradient.length == reference.length;

        double max = Math.abs(gradient[0] - reference[0]) / relative;
        for (int i = 1; i < gradient.length; ++i) {
            double abs = Math.abs(gradient[i] - reference[i]) / relative;
            if (abs > max) {
                max = abs;
            }
        }

        return max;
    }

    private static double innerProduct(double[] x, double[] y) {
        assert x.length == y.length;

        double product = 0.0;
        for (int i = 0; i < x.length; ++i) {
            product += x[i] * y[i];
        }

        return product;
    }
}
