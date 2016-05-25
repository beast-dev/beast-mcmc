/*
 * ExponentialSawtooth.java
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
 * This class models a sawtooth demographic in which exponential growth is punctuated by frequent periodic crashes.
 * 
 * @version $Id: ExponentialSawtooth.java,v 1.3 2005/04/11 11:43:03 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class ExponentialSawtooth extends ExponentialGrowth {
	
	/**
	 * Construct demographic model with default settings
	 */
	public ExponentialSawtooth(Type units) {
	
		super(units);
	}
	
	public final double getWavelength() { return wavelength; }
	
	public final void setWavelength(double t) { 
		if (t <= 0) throw new IllegalArgumentException();
		wavelength = t; 
	}
	
	public final double getOffset() { return offset; }
	
	public final void setOffset(double offset) {
		if (offset < 0 || offset >= 1.0) {
			throw new IllegalArgumentException();
		}
		this.offset = offset;
	}	
			
	// Implementation of abstract methods

	public double getDemographic(double t) {

		t += offset * wavelength;
		
		// rescale t so that 0 <= t < wavelength
		int cycle = (int)Math.floor(t / wavelength);
		t -= (cycle * wavelength);

		return super.getDemographic(t);
	}

	public double getIntensity(double t) {
	
		double absOffset = offset*wavelength;
		
		if (t < wavelength-absOffset) {
			// calculate intensity of first partial offset
			return super.getIntensity(t+absOffset) - super.getIntensity(absOffset);
		}
		
		// calculate intensity of first epoch:
		double intensity = super.getIntensity(wavelength) - super.getIntensity(absOffset);
		
		t -= (wavelength-absOffset);
		
		// calculate intensity of all full cycles
		int cycles = (int)Math.floor(t / wavelength);
		intensity += cycles * super.getIntensity(wavelength);
		t -= (cycles * wavelength);

		// calculate intensity of last partial cycle
		intensity += super.getIntensity(t);
		
		return intensity;
	}
	
	public double getInverseIntensity(double x) {
		
		throw new UnsupportedOperationException();
	}
	
	public int getNumArguments() {
		return 4;
	}
	
	public String getArgumentName(int n) {
		
		switch (n) {
			case 0: return "N0";
			case 1: return "r";
			case 2: return "wavelength";
			case 3: return "offset";
			default: throw new IllegalArgumentException();
		}
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getWavelength();
			case 3: return getOffset();
			default: throw new IllegalArgumentException();
		}
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setWavelength(value); break;
			case 3: setOffset(value); break;
			default: throw new IllegalArgumentException();
		}
	}

	public DemographicFunction getCopy() {
		ExponentialSawtooth df = new ExponentialSawtooth(getUnits());
		df.setN0(getN0());
		df.setGrowthRate(getGrowthRate());
		df.setWavelength(getWavelength());
		df.setOffset(getOffset());
		
		return df;
	}
	
	public static void main(String[] args) {
	
		double N0 = Double.parseDouble(args[0]);
		double growthRate = Double.parseDouble(args[1]);
		double wavelength = Double.parseDouble(args[2]);
		double offset = Double.parseDouble(args[3]);
		
	
		ExponentialSawtooth est = new ExponentialSawtooth(Type.SUBSTITUTIONS);
		est.setN0(N0);
		est.setGrowthRate(growthRate);
		est.setWavelength(wavelength);
		est.setOffset(offset);
	
		for (double time = 0; time < 20; time+=0.1) {
			System.out.println(time + "\t" + est.getDemographic(time) + "\t" + est.getIntensity(time));
		}
	}
	

	//
	// private stuff
	//

	private double wavelength;
	private double offset;
}
