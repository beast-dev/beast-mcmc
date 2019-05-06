/*
 * AbstractCoalescentLikelihood.java
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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.Intervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Units;
import dr.inference.model.*;
import dr.math.Binomial;


/**
 * Forms a base class for a number of coalescent likelihood calculators.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CoalescentLikelihood.java,v 1.43 2006/07/28 11:27:32 rambaut Exp $
 */
public abstract class AbstractCoalescentLikelihood extends AbstractModelLikelihood implements Units /*, CoalescentIntervalProvider*/ {

    // PUBLIC STUFF

    public AbstractCoalescentLikelihood( String name, IntervalList intervalList){
        super(name);

        this.intervalList = intervalList;

        addStatistic(new DeltaStatistic());
        if (intervalList instanceof Model) {
            addModel((Model)intervalList);
        }
    }


    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // something the likelihood is listening to has changed so flag the likelihood to update

        likelihoodKnown = false;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    protected abstract double calculateLogLikelihood();


    public IntervalList getIntervalList() {
            return intervalList;
    }

//    public double getCoalescentEventsStatisticValue(int i) {
//        if (i == 0) {
//            for (int j = 0; j < coalescentEventStatisticValues.length; j++) {
//                coalescentEventStatisticValues[j] = 0.0;
//            }
//            int counter = 0;
//            for (int j = 0; j < getCoalescentIntervalDimension(); j++) {
//                if (getCoalescentIntervalType(j) == IntervalType.COALESCENT) {
//                    this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getCoalescentIntervalLineageCount(j) * (getCoalescentIntervalLineageCount(j) - 1.0)) / 2.0;
//                    counter++;
//                } else {
//                    this.coalescentEventStatisticValues[counter] += getCoalescentInterval(j) * (getCoalescentIntervalLineageCount(j) * (getCoalescentIntervalLineageCount(j) - 1.0)) / 2.0;
//                }
//            }
//        }
//        return coalescentEventStatisticValues[i];
//    }

    public String toString() {
        return Double.toString(logLikelihood);

    }

    // ****************************************************************
    // Inner classes
    // ****************************************************************

    public class DeltaStatistic extends Statistic.Abstract {

        public DeltaStatistic() {
            super("delta");
        }

        public int getDimension() {
            return 1;
        }

        public double getStatisticValue(int i) {
            throw new RuntimeException("Not implemented");
//			return IntervalList.Utils.getDelta(intervals);
        }

    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    private IntervalList intervalList = null;

    protected double logLikelihood;
    protected double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    protected boolean storedLikelihoodKnown = false;

//    private double[] coalescentEventStatisticValues;
}