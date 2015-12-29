/*
 * GreatDelugeCriterion.java
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
 * This class encapsulates the acception criterion for a GreatDeluge search.
 *
 * @author Andrew Rambaut
 *
 * @version $Id: GreatDelugeCriterion.java,v 1.3 2006/06/13 03:50:54 alexei Exp $
 */
public class GreatDelugeCriterion implements Acceptor  {

	protected double rate;
	protected double bound = Double.NEGATIVE_INFINITY;
    protected double maxScore = Double.NEGATIVE_INFINITY;


    public GreatDelugeCriterion(double rate) {
		this.rate = rate;
	}

	public boolean accept(double oldScore, double newScore, double hastingsRatio, double[] logr) {

        // HACK HACK HACK
        if (newScore > maxScore) {
            maxScore = newScore;
            bound = maxScore - 1;
        }

        if (newScore > bound) {
			
			//if (bound == Double.NEGATIVE_INFINITY) {
			//	bound = oldScore * 2;
			//}
			
			//bound += (newScore - bound) * rate;

            return true;
		}
		
		return false;
		
	}
}
