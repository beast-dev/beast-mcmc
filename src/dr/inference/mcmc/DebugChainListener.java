/*
 * DebugChainListener.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.mcmc;

/**
 * ${CLASS_NAME}
 *
 * @author Andrew Rambaut
 * @author Guy Baele
 * @version $Id$
 *
 * $HeadURL$
 *
 * $LastChangedBy$
 * $LastChangedDate$
 * $LastChangedRevision$
 */

import dr.evomodel.tree.TreeDebugLogger;
import dr.evomodel.tree.TreeLogger;
import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Model;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;

public class DebugChainListener implements MarkovChainListener {

    private MCMC mcmc;

    public DebugChainListener(MCMC mcmc, Logger[] loggers, final long writeState, final boolean isRepeating, final String fileName) {
        this.mcmc = mcmc;
        this.loggers = loggers;

        //add TreeLoggers using TreeDebugLoggers
        int counter = 0;
        ArrayList<TreeDebugLogger> additionalLoggers = new ArrayList<TreeDebugLogger>();
        for (int i = 0; i < loggers.length; i++) {
            if (loggers[i] instanceof TreeLogger) {
                additionalLoggers.add(new TreeDebugLogger(fileName + ".tree." + counter, (TreeLogger) loggers[i]));
                counter++;
            }
        }

        Logger[] newLoggers = new Logger[this.loggers.length+additionalLoggers.size()];
        for (int i = 0; i < this.loggers.length; i++) {
            newLoggers[i] = this.loggers[i];
        }
        counter = 0;
        for (int i = this.loggers.length; i < newLoggers.length; i++) {
            newLoggers[i] = additionalLoggers.get(counter);
            counter++;
        }
        this.loggers = newLoggers;

        this.writeState = writeState;
        this.isRepeating = isRepeating;
        this.fileName = fileName;
    }

    // MarkovChainListener interface *******************************************

    /**
     * Called to update the current model keepEvery states.
     */
    public void currentState(long state, Model currentModel) {
        if (state == writeState || (isRepeating && state > 0 && (state % writeState == 0))) {
            String timeStamp = new SimpleDateFormat("yyyy.MM.dd.HH.mm.ss").format(Calendar.getInstance().getTime());
            mcmc.getMarkovChain().getLikelihood().makeDirty();
            double lnL = mcmc.getMarkovChain().getCurrentScore();

            String fileName = (this.fileName != null ? this.fileName : "beast_debug_" + timeStamp);
            DebugUtils.writeStateToFile(new File(fileName), loggers, state, lnL, mcmc.getOperatorSchedule());
        }
    }

    /**
     * Called when a new new best posterior state is found.
     */
    public void bestState(long state, Model bestModel) { }

    /**
     * Cleans up when the chain finishes (possibly early).
     */
    public void finished(long chainLength) {
        currentState(chainLength, null);
    }

    private final long writeState;
    private final boolean isRepeating;
    private final String fileName;
    private Logger[] loggers;
}
