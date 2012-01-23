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

package dr.evomodel.epidemiology;

import dr.app.tools.NexusExporter;
import dr.evolution.tree.BranchRates;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.evolution.tree.TreeTraitProvider;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;

import java.text.NumberFormat;
import java.util.*;

/**
 * A logger that timeseries of dynamical systems.
 *
 * @author Trevor Bedford
 * @author Andrew Rambaut
 * @version $Id: TreeLogger.java,v 1.25 2006/09/05 13:29:34 rambaut Exp $
 */
public class DynamicalLogger extends MCLogger {

    private SIRModel model;

    /**
     * Interface to indicate when to log a tree
     */
    public interface LogUpon {
        /**
         *
         * @param state
         * @return  True if log tree of this state.
         */
       boolean logNow(long state);
    }

    public DynamicalLogger(SIRModel model, LogFormatter formatter, int logEvery) {

        super(formatter, logEvery, false);
        this.model = model;

    }

    public void startLogging() {
        logLine("State");
    }

    public void log(long state) {

        final boolean doIt = (logEvery < 0 || ((state % logEvery) == 0));

        if ( doIt ) {
            StringBuffer buffer = new StringBuffer("");
            buffer.append(state);
            for (double t=0; t < 5; t += 0.01) {
                double v = model.getInfecteds(t);
                buffer.append("\t" + v);
            }
            logLine(buffer.toString());
        }
    }

    public void stopLogging() {
        super.stopLogging();
    }

    public SIRModel getModel() {
		return model;
	}

	public void setModel(SIRModel model) {
		this.model = model;
	}

}