/*
 * TwoEpochDemographic.java
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

package dr.evolution.coalescent;

/**
 * This class models a demographic function that contains two epochs.
 * 
 * @author Alexei Drummond
 *
 * @version $Id: TwoEpochDemographic.java,v 1.7 2005/05/24 20:25:56 rambaut Exp $
 */
public class TwoEpochDemographic extends DemographicFunction.Abstract {

    private DemographicFunction epoch1;
    private DemographicFunction epoch2;
    private double transitionTime;

	/**
	 * Construct demographic model with default settings.
	 */
	public TwoEpochDemographic(DemographicFunction epoch1, DemographicFunction epoch2, Type units) {

		super(units);

        this.epoch1 = epoch1;
        this.epoch2 = epoch2;
	}

	/**
	 * @return the time of the transition between the two demographic functions.
	 */
	public final double getTransitionTime() { return transitionTime; }

    public final void setTransitionTime(double t) {
        if (t < 0.0 || t > Double.MAX_VALUE) {
            throw new IllegalArgumentException("transition time out of bounds.");
        }
        transitionTime = t;
    }

    public final DemographicFunction getFirstEpochDemography() {
        return epoch1;
    }

    public final DemographicFunction getSecondEpochDemography() {
        return epoch2;
    }

	// Implementation of abstract methods

	public final double getDemographic(double t) {

		if (t < transitionTime) {
            return epoch1.getDemographic(t);
        } else {
            return epoch2.getDemographic(t-transitionTime);
        }
	}

	public final double getIntensity(double t) {
	
		if (t < transitionTime) {
            return epoch1.getIntensity(t);
        } else {
            return epoch1.getIntensity(transitionTime) + epoch2.getIntensity(t-transitionTime);
        }
	}

    /**
     * By overriding this function we ensure that the integral can still be found when one of the sub-demographic models does not implement getIntensity.
     * @param start
     * @param finish
     * @return the integral of the demographic function between start and finish
     */
    public final double getIntegral(double start, double finish) {

		if (start < transitionTime) {
            if (finish < transitionTime) {
                return epoch1.getIntegral(start, finish);
            } else {
                return epoch1.getIntegral(start, transitionTime) + epoch2.getIntegral(0, finish-transitionTime);
            }
        } else {
            return epoch2.getIntegral(start-transitionTime, finish-transitionTime);
        }
	}

    /**
     * @param x a value of intensity
     * @return the time in the past that corresponds to the given intensity.
     */
	public final double getInverseIntensity(double x) {

        double time = epoch1.getInverseIntensity(x);
        if (time < transitionTime) {
            return time;
        }
        x -= epoch1.getIntensity(transitionTime);
        return transitionTime + epoch2.getInverseIntensity(x);
  	}
	
	public int getNumArguments() {
		return epoch1.getNumArguments() + epoch2.getNumArguments() + 1;
	}
	
	public final String getArgumentName(int n) {
		
		if (n < epoch1.getNumArguments()) {
            return epoch1.getArgumentName(n);
        }
        n -= epoch1.getNumArguments();
        if (n < epoch2.getNumArguments()) {
            return epoch2.getArgumentName(n);
        }
        n -= epoch2.getNumArguments();
        if (n == 0) return "transitionTime";
        throw new IllegalArgumentException();
	}
	
	public final double getArgument(int n) {

        if (n < epoch1.getNumArguments()) {
            return epoch1.getArgument(n);
        }
        n -= epoch1.getNumArguments();
        if (n < epoch2.getNumArguments()) {
            return epoch2.getArgument(n);
        }
        n -= epoch2.getNumArguments();
        if (n == 0) return transitionTime;
        throw new IllegalArgumentException();
	}
	
	public final void setArgument(int n, double value) {

        if (n < epoch1.getNumArguments()) {
            epoch1.setArgument(n, value);
        }
        n -= epoch1.getNumArguments();
        if (n < epoch2.getNumArguments()) {
            epoch2.setArgument(n, value);
        }
        n -= epoch2.getNumArguments();
        if (n == 0) transitionTime = value;
        throw new IllegalArgumentException();
	}

    /**
     * @return the lower bound of the nth argument of this function.
     */
    public final double getLowerBound(int n) {

        if (n < epoch1.getNumArguments()) {
            return epoch1.getLowerBound(n);
        }
        n -= epoch1.getNumArguments();
        if (n < epoch2.getNumArguments()) {
            return epoch2.getLowerBound(n);
        }
        n -= epoch2.getNumArguments();
        if (n == 0) return 0.0;
        throw new IllegalArgumentException();
    }

    /**
     * @return the upper bound of the nth argument of this function.
     */
    public final double getUpperBound(int n) {

        if (n < epoch1.getNumArguments()) {
            return epoch1.getUpperBound(n);
        }
        n -= epoch1.getNumArguments();
        if (n < epoch2.getNumArguments()) {
            return epoch2.getUpperBound(n);
        }
        n -= epoch2.getNumArguments();
        if (n == 0) return Double.MAX_VALUE;
        throw new IllegalArgumentException();
    }

    public final DemographicFunction getCopy() {
		TwoEpochDemographic df = new TwoEpochDemographic(epoch1, epoch2, getUnits());
		df.setTransitionTime(transitionTime);
		return df;
	}
}
