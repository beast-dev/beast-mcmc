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

package dr.inference.checkpoint;

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

import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.mcmc.DebugUtils;
import dr.inference.mcmc.MCMC;
import dr.inference.model.Model;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class CheckPointChainListener implements MarkovChainListener {

    private MarkovChain markovChain;

    public CheckPointChainListener(MarkovChain markovChain, final long writeState, final boolean isRepeating, final String fileName) {
        this.markovChain = markovChain;
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
            markovChain.getLikelihood().makeDirty();
            //double lnL = mcmc.getMarkovChain().getLikelihood().getLogLikelihood();
            double lnL = markovChain.getCurrentScore();

            String fileName = (this.fileName != null ? this.fileName : "beast_debug_" + timeStamp);
            //DebugUtils.writeStateToFile(new File(fileName), state, lnL, markovChain.getOperatorSchedule());
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
}
