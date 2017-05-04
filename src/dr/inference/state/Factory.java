/*
 * Factory.java
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

package dr.inference.state;

import dr.inference.markovchain.MarkovChainListener;

public abstract class Factory {

    /**
     * Get an instance of StateLoader that will provide the initial state for the
     * chain. This is likely to be a previously saved state.
     * @return the StateLoader
     */
    public abstract StateLoader getInitialStateLoader();

    /**
     * Get a list of MarkovChainListeners that will save the state at a particular point or
     * on a regular interval.
     * @return the array
     */
    public abstract MarkovChainListener[] getStateSaverChainListeners();

    // Set this to a concrete instance to provide these classes to the MarkovChain
    public static Factory INSTANCE;
}
