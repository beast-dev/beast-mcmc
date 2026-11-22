/*
 * IndependentChainSet.java
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

package dr.inference.mcmcmc;

import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.ChainContext;
import dr.inference.mcmc.MCMC;

import java.util.Collections;
import java.util.concurrent.Callable;
import java.util.concurrent.Semaphore;

/**
 * A single chain run alongside others. It runs in segments so that, when there are fewer
 * chain threads than chains, the chains take turns and progress together.
 *
 * @author Andrew Rambaut
 */
public class IndependentChainSet extends AbstractChainSet {

    public static final int SEGMENTS_PER_CHAIN = 100;

    public IndependentChainSet(int index, MCMC mcmc, ChainContext context, ChainContext driverContext) {
        super(index, new MCMC[]{mcmc}, new ChainContext[]{context}, driverContext);
        this.mcmc = mcmc;
        this.segmentLength = Math.max(1, mcmc.getChainLength() / SEGMENTS_PER_CHAIN);
    }

    @Override
    public void start(Semaphore permits) {
        super.start(permits);
        runOnChains(Collections.<Callable<Boolean>>singletonList(() -> {
            mcmc.start();
            return true;
        }));
    }

    @Override
    public boolean step() {
        return runOnChains(Collections.<Callable<Boolean>>singletonList(() -> mcmc.runSegment(segmentLength))).get(0);
    }

    @Override
    protected void finishChains() {
        runOnChains(Collections.<Callable<Boolean>>singletonList(() -> {
            mcmc.finish();
            return true;
        }));
    }

    @Override
    public MCMC getColdMCMC() {
        return mcmc;
    }

    @Override
    public MarkovChain getColdChain() {
        return mcmc.getMarkovChain();
    }

    @Override
    public long getCurrentState() {
        return mcmc.getMarkovChain().getCurrentLength();
    }

    @Override
    public String getStatus() {
        return "";
    }

    private final MCMC mcmc;
    private final long segmentLength;
}
