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

import dr.app.tools.NexusExporter;
import dr.evolution.tree.*;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.model.Parameter;

import java.text.NumberFormat;
import java.util.*;

/**
 * A logger that logs tree and clade frequencies.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class DemographicLogger extends MCLogger {

    private Tree tree;
    private Parameter[] parameters;

    public DemographicLogger(Tree tree, Parameter[] parameters,
                      LogFormatter formatter, int logEvery) {

        super(formatter, logEvery, false);

        this.tree = tree;
        this.parameters = parameters;
    }

    public void startLogging() {
        String title = "";
        for (LogFormatter formatter : formatters) {
            formatter.startLogging(title);
        }

        if (title != null) {
            logHeading(title);
        }

        if (logEvery > 0) {
            final int columnCount = getColumnCount();
            String[] labels = new String[columnCount + 1];

            labels[0] = "state";

            for (int i = 0; i < columnCount; i++) {
                labels[i + 1] = getColumnLabel(i);
            }

            logLabels(labels);
        }
    }

    public void log(int state) {

        if (logEvery > 0 && (state % logEvery == 0)) {

            final int columnCount = getColumnCount();

            String[] values = new String[columnCount];

            values[0] = Integer.toString(state);

            for (int i = 0; i < columnCount; i++) {
                values[i + 1] = getColumnFormatted(i);
            }

            logValues(values);
        }
    }

    public void stopLogging() {
        for (LogFormatter formatter : formatters) {
            formatter.stopLogging();
        }
    }


}