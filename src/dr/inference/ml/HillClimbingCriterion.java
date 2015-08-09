/*
 * HillClimbingCriterion.java
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

package dr.inference.ml;

import dr.inference.markovchain.Acceptor;

/**
 * Simple hill climbing.
 *
 * @author Alexei Drummond
 *
 * @version $Id: HillClimbingCriterion.java,v 1.1 2005/03/15 15:55:10 alexei Exp $
 */
public class HillClimbingCriterion implements Acceptor  {

	protected double bound = Double.NEGATIVE_INFINITY;

	public HillClimbingCriterion() {}

	public boolean accept(double oldScore, double newScore, double hastingsRatio, double[] logr) {
		
		if (newScore > bound) {

			bound = newScore;
			
			return true;
		}
		
		return false;
		
	}
}
