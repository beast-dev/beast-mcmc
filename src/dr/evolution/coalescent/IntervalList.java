/*
 * IntervalList.java
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

import dr.evolution.util.Units;

/**
 * An interface for a set of coalescent intevals.
 *
 * @version $Id: IntervalList.java,v 1.7 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public interface IntervalList extends Units {

	/**
	 * get number of intervals
	 */
	int getIntervalCount();

	/**
	 * get the total number of sampling events. 
	 */
	int getSampleCount();

	/**
	 * Gets an interval.
	 */
	double getInterval(int i);

	/**
	 * Returns the number of uncoalesced lineages within this interval.
	 * Required for s-coalescents, where new lineages are added as
	 * earlier samples are come across.
	 */
	int getLineageCount(int i);

	/**
	 * Returns the number coalescent events in an interval
	 */
	int getCoalescentEvents(int i);

	/**
	 * Returns the type of interval observed.
	 */
	IntervalType getIntervalType(int i);

	/**
	 * get the total duration of these intervals. 
	 */
	double getTotalDuration();

	/**
	 * Checks whether this set of coalescent intervals is fully resolved
	 * (i.e. whether is has exactly one coalescent event in each
	 * subsequent interval)
	 */
	boolean isBinaryCoalescent();

	/**
	 * Checks whether this set of coalescent intervals coalescent only
	 * (i.e. whether is has exactly one or more coalescent event in each
	 * subsequent interval)
	 */
	boolean isCoalescentOnly();


	public class Utils {	

		/**
		 * @return the number of lineages at time t.
		 * @param t the time that you are counting the number of lineages
		 */
		public static int getLineageCount(IntervalList intervals, double t) {
		
			int i = 0;
			while (i < intervals.getIntervalCount() && t > intervals.getInterval(i)) { 
				t -= intervals.getInterval(i);
				i+= 1; 
			}
			if (i == intervals.getIntervalCount()) return 1;
			return intervals.getLineageCount(i);
		}

		/**
		 * @return the delta parameter of Pybus et al (Node spread statistic)
		 * @param intervals the intervals for which the delta parameter is calculated.
		 */
		public static double getDelta(IntervalList intervals) {
			
			// Assumes ultrametric tree!
			if (!intervals.isCoalescentOnly()) {
				throw new IllegalArgumentException("Assumes ultrametric tree!");
			}
			
			int n = intervals.getIntervalCount();
						
			int numTips = n + 1;
			
			double transTreeDepth = 0.0;
			double cumInts = 0.0;
			double sum = 0.0;
			
			// transform intervals
			for (int j=0; j<n; j++) { // move from tips to root
			
				double transInt = intervals.getInterval(j) * 
									dr.math.Binomial.choose2(intervals.getLineageCount(j)); // coalescent version
				//intLenCopy[j] = getInterval(j)*getLineageCount(j); // birth-death version

				// don't include the last interval so put this before...
				sum += cumInts;

				// ...incrementing the cumInts
				cumInts += transInt;
				
				transTreeDepth += transInt;
			}
			
			double halfTreeDepth = transTreeDepth / 2.0;
						
			sum *= (1.0/(numTips-2.0));
			double top = halfTreeDepth - sum;
			double bottom = transTreeDepth * Math.sqrt((1.0/(12.0*(numTips-2.0))));

			return (top / bottom);
		}
	}
}