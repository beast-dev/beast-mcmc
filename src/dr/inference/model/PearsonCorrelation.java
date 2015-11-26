/*
 * PearsonCorrelation.java
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

package dr.inference.model;

import dr.stats.DiscreteStatistics;

/**
 * @author Simon Greenhill
 * @author Alexei Drummond
 */
public class PearsonCorrelation extends Statistic.Abstract {

    Parameter X, Y;
    boolean log = false;


    public PearsonCorrelation(Parameter X, Parameter Y, boolean log) {

        if (X.getDimension() != Y.getDimension()) throw new IllegalArgumentException();
        this.X = X;
        this.Y = Y;
        this.log = log;
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {

        double[] xvalues = X.getParameterValues();
        double[] yvalues = Y.getParameterValues();

        if (log) {
            for (int i = 0; i < xvalues.length; i++) {
                xvalues[i] = Math.log(xvalues[i]);
                yvalues[i] = Math.log(yvalues[i]);
            }
        }

        double meanX = DiscreteStatistics.mean(xvalues);
        double meanY = DiscreteStatistics.mean(yvalues);
        double stdevX = DiscreteStatistics.stdev(xvalues);
        double stdevY = DiscreteStatistics.stdev(yvalues);

        double corr = 0;
        for (int i = 0; i < xvalues.length; i++) {
            double deviateX = xvalues[i] - meanX;
            double deviateY = yvalues[i] - meanY;
            corr += deviateX*deviateY;
        }
        corr /= X.getDimension();

        corr /= (stdevX*stdevY);
        return corr;


    }
}
