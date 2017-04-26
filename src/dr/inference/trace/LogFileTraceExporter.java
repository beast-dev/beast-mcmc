/*
 * LogFileTraceExporter.java
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

package dr.inference.trace;

import dr.util.TabularData;

import java.io.File;
import java.io.IOException;

/**
 * Export trace analysis data such as mean,median,HPD and ESS of trace variables.
 *
 * @author Joseph Heled
 *         Date: 25/10/2007
 */
public class LogFileTraceExporter extends TabularData {
    private final LogFileTraces analysis;
    private final String[] rows = {"mean", "median", "hpdLower", "hpdUpper", "ESS"};
    TraceCorrelation[] distributions;

    public LogFileTraceExporter(File file, int burnin) throws TraceException, IOException {

        analysis = new LogFileTraces(file.getCanonicalPath(), file);
        analysis.loadTraces();
        if (burnin >= 0) {
            analysis.setBurnIn(burnin);
        }

        distributions = new TraceCorrelation[nColumns()];
    }

    public int nColumns() {
        return analysis.getTraceCount();
    }

    public String columnName(int nColumn) {
        return analysis.getTraceName(nColumn);
    }

    public int nRows() {
        return rows.length;
    }

    public Object data(int nRow, int nColumn) {
        // read on demand
        if (distributions[nColumn] == null) {
            analysis.analyseTrace(nColumn);
            distributions[nColumn] = analysis.getCorrelationStatistics(nColumn);
        }

        TraceCorrelation distribution = distributions[nColumn];

        switch (nRow) {
            case 0: {
                return distribution.getMean();
            }
            case 1: {
                return distribution.getMedian();
            }
            case 2: {
                return distribution.getLowerHPD();
            }
            case 3: {
                return distribution.getUpperHPD();
            }
            case 4: {
                return distribution.getESS();
            }
        }

        return null;
    }

}
