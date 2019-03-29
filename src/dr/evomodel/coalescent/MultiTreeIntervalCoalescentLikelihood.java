/*
 * MulitTreeIntervalCoalescentLikelihood.java
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
package dr.evomodel.coalescent;

import dr.evolution.coalescent.Coalescent;
import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.util.Units;
import dr.evomodelxml.coalescent.CoalescentLikelihoodParser;
import dr.inference.model.Model;
import java.util.logging.Logger;


/**
 * A likelihood function for the coalescent. Takes an interval list and a demographic model.
 *
 * This class is intended to replace CoalescentLikelihood which took a tree.
 * Parts of that class were derived from C++ code provided by Oliver Pybus.
 *
 * @author JT McCrone
 */

public class MultiTreeIntervalCoalescentLikelihood extends AbstractCoalescentLikelihood implements Units {

    // PUBLIC STUFF

    public MultiTreeIntervalCoalescentLikelihood(MultiTreeIntervals multiTreeIntervals, DemographicModel demoModel){

        super(CoalescentLikelihoodParser.COALESCENT_LIKELIHOOD, multiTreeIntervals);
        this.multiTreeIntervals = multiTreeIntervals;
        this.demoModel = demoModel;

        addModel(demoModel);
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {

        DemographicFunction demoFunction = demoModel.getDemographicFunction();
        double lnL =  Coalescent.calculateLogLikelihood(multiTreeIntervals, demoFunction, demoFunction.getThreshold());
        if (Double.isNaN(lnL) || Double.isInfinite(lnL)) {
            Logger.getLogger("warning").severe("CoalescentLikelihood for " + demoModel.getId() + " is " + Double.toString(lnL));
        }

        return lnL;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    @Override
    public void makeDirty() {
        //multiTreeIntervals.setIntervalsUnknown();
        likelihoodKnown = false;
    }

    @Override
    public IntervalList getIntervals() {
        return multiTreeIntervals;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u)
    {
        demoModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits()
    {
        return demoModel.getUnits();
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /** The demographic model. */
    private DemographicModel demoModel = null;

    private MultiTreeIntervals multiTreeIntervals;
}

