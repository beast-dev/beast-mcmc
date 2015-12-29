/*
 * Continuous.java
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

package dr.evolution.continuous;

/**
 * a class for a single continuous value.
 *
 * @version $Id: Continuous.java,v 1.5 2006/06/18 16:20:58 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class Continuous implements Contrastable {
	
	public Continuous(double value) {
		this.value = value;
	}
	
	public double getDifference(Contrastable cont) {
	
		if (cont instanceof Continuous) {
			return value - ((Continuous)cont).value;
		} else throw new IllegalArgumentException("Expected a continuous parameter");
	}
	
	public Contrastable getWeightedMean(double weight1, Contrastable cont1, double weight2, Contrastable cont2) {
		
		double value = 0.0;
		
		value += ((Continuous)cont1).value * weight1;
		value += ((Continuous)cont2).value * weight2;
		
		value /= (weight1 + weight2);
		return new Continuous(value);
	}
	
	public double getValue() { return value; }

    public String toString() { return ""+getValue(); }

    private double value;
}
