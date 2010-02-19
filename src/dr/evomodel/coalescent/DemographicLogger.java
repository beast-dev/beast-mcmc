/*
 * TreeLogger.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.coalescent;

import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.Logger;

import java.util.ArrayList;
import java.util.List;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class DemographicLogger implements Logger {
    private DemographicReconstructor reconstructor;

    public DemographicLogger(DemographicReconstructor reconstructor,
                             LogFormatter formatter, int logEvery) {

        addFormatter(formatter);
        this.logEvery = logEvery;

        this.reconstructor = reconstructor;
    }

    public void startLogging() {
        String title = reconstructor.getTitle();
        DemographicReconstructor.ChangeType type = reconstructor.getChangeType();
        for (LogFormatter formatter : formatters) {
            formatter.startLogging(title);
            formatter.logHeading(title);
            formatter.logHeading("type=" + type.name());
        }
    }

    public void log(int state) {

        if (logEvery > 0 && (state % logEvery == 0)) {

            double[] intervals = reconstructor.getIntervals();
            double[] popSizes = reconstructor.getPopSizes();

            final int columnCount = 1 + intervals.length + popSizes.length;

            String[] values = new String[columnCount];

            int k = 0;

            values[k] = Integer.toString(state);
            k++;

            for (double interval : intervals) {
                values[k] = Double.toString(interval);
                k++;
            }

            for (double popSize : popSizes) {
                values[k] = Double.toString(popSize);
                k++;
            }

            for (LogFormatter formatter : formatters) {
                formatter.logValues(values);
            }
        }
    }

    public void stopLogging() {
        for (LogFormatter formatter : formatters) {
            formatter.stopLogging();
        }
    }

    protected final void addFormatter(LogFormatter formatter) {
        formatters.add(formatter);
    }

    protected int logEvery = 0;
    protected List<LogFormatter> formatters = new ArrayList<LogFormatter>();

}