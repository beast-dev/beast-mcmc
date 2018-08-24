/*
 * MLLogger.java
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

package dr.inference.loggers;

import dr.inference.model.Likelihood;

/**
 * A logger that stores maximum likelihood states.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: MLLogger.java,v 1.21 2005/07/27 22:09:21 rambaut Exp $
 */
public class MLLogger extends MCLogger {

    private final Likelihood likelihood;
    private double bestLikelihood;
    private long bestState;
    private String[] bestValues = null;
    private long logEvery = 0;

    public MLLogger(Likelihood likelihood, LogFormatter formatter, long logEvery) {

        super(formatter, logEvery, false);

        this.likelihood = likelihood;
    }

    public void startLogging() {
        bestLikelihood = Double.NEGATIVE_INFINITY;
        bestState = 0;
        bestValues = new String[getColumnCount()];

        if (logEvery > 0) {
            String[] labels = new String[getColumnCount() + 1];

            labels[0] = "state";

            for (int i = 0; i < getColumnCount(); i++) {
                labels[i + 1] = getColumnLabel(i);
            }

            logLabels(labels);
        }

        super.startLogging();
    }

    public void log(long state) {

        double lik;

        lik = likelihood.getLogLikelihood();

        if (lik > bestLikelihood) {

            for (int i = 0; i < getColumnCount(); i++) {
                bestValues[i] = getColumnFormatted(i);
            }

            bestState = state;
            bestLikelihood = lik;

            if (logEvery == 1) {

                String[] values = new String[getColumnCount() + 1];

                values[0] = Long.toString(bestState);

                System.arraycopy(bestValues, 0, values, 1, getColumnCount());

                logValues(values);
            }
        }

        if (logEvery > 1 && (state % logEvery == 0)) {

            String[] values = new String[getColumnCount() + 1];

            values[0] = Long.toString(bestState);

            System.arraycopy(bestValues, 0, values, 1, getColumnCount());

            logValues(values);
        }
    }

    public void stopLogging() {
        final int columnCount = getColumnCount();
        String[] values = new String[columnCount + 2];

        values[0] = Long.toString(bestState);
        values[1] = Double.toString(bestLikelihood);

        System.arraycopy(bestValues, 0, values, 2, columnCount);

        if (logEvery > 0) {
            logValues(values);
        } else {
            String[] labels = new String[columnCount + 2];

            labels[0] = "state";
            labels[1] = "ML";

            for (int i = 0; i < columnCount; i++) {
                labels[i + 2] = getColumnLabel(i);
            }

            logLabels(labels);
            logValues(values);
        }

        super.stopLogging();
    }
}
