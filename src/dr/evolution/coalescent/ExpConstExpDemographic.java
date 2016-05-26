/*
 * ExpConstExpDemographic.java
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
 * This class models a two-stage exponential growth with a constant plateau separating them.
 * 
 * @version $Id: ExpConstExpDemographic.java,v 1.4 2006/09/28 20:45:48 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class ExpConstExpDemographic extends ExponentialGrowth {
	
	//
	// Public stuff
	//

	/**
	 * Construct demographic model with default settings
	 */
	public ExpConstExpDemographic(Type units) {
	
		super(units);
	}

	/**
	 * returns the positive-valued decline rate
	 */
	public final double getGrowthRate2() { return r2; }
	
	/**
	 * sets the decline rate.
	 */
	public void setGrowthRate2(double r2) {
		this.r2 = r2;
	}

    /**
     * @return the time (back from the present) to the start of the modern growth phase
     */
    public final double getTime1() { return time1; }

    /**
     * @return the length of time in the plateau phase
     */
    public final double getPlateauTime() { return plateauTime; }

	public final void setTime1(double t) {
		if (t <= 0) throw new IllegalArgumentException();
		time1 = t;
	}

    public final void setPlateauTime(double t) {
        if (t <= 0) throw new IllegalArgumentException();
        plateauTime = t;
    }


	// Implementation of abstract methods

	public double getDemographic(double t) {

        double r = getGrowthRate();
        double r2 = getGrowthRate2();

		if (t < time1) {
			return getN0() * Math.exp(-t * r);
		}

        double plateauHeight = getN0() * Math.exp(time1 * r);

        if (t < (time1 + plateauTime)) {
		    return plateauHeight;
        } else {

			t -= (time1 + plateauTime);
		
			if (r2 == 0) {
				return plateauHeight;
			} else {
				return plateauHeight * Math.exp(-t * r2);
			}
		}
	}

	public double getIntensity(double t) {
	
		double r = getGrowthRate();
        double r2 = getGrowthRate2();
        double plateauHeight = getN0() * Math.exp(time1 * r);
        if (t < time1) {
			return super.getIntensity(t);
		} else {

            double modernGrowthPhaseIntensity = super.getIntensity(time1);

            if (t < (time1 + plateauTime)) {

                double constantIntensity = (t-time1)/plateauHeight;

                return modernGrowthPhaseIntensity + constantIntensity;

            } else {

                double constantIntensity = plateauTime/plateauHeight;

                t -= (time1 + plateauTime);

                double ancientGrowthPhaseIntensity;

                if (r2 == 0) {
                    ancientGrowthPhaseIntensity = t/plateauHeight;
                } else {

                    //(Math.exp(t*r)-1.0)/getN0()/r;
                    ancientGrowthPhaseIntensity = (Math.exp(t*r2)-1.0)/plateauHeight/r2;
                }

                double intensity = modernGrowthPhaseIntensity +
                        constantIntensity + ancientGrowthPhaseIntensity;


                if (intensity < minIntensity) {
                    //System.out.println("min intensity = " + intensity);
                    minIntensity = intensity;
                }
                if (intensity > maxIntensity) {
                    /*System.out.println("max intensity = " + intensity);

                    System.out.println("  mgf       = " + modernGrowthPhaseIntensity);
                    System.out.println("  c         = " + constantIntensity);
                    System.out.println("  agf       = " + ancientGrowthPhaseIntensity);
                    System.out.println("  r2        = " + r2);
                    System.out.println("  N1        = " + plateauHeight);
                    System.out.println("  t         = " + t);
                    System.out.println("  t*r2      = " + t*r2);
                    System.out.println("  exp(t*r2) = " + Math.exp(t*r2));
*/

                    maxIntensity = intensity;
                }

                return intensity;
            }
		}
	}
    private static double minIntensity = Double.MAX_VALUE;
    private static double maxIntensity = Double.MIN_VALUE;

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
			case 2: return "d";
			case 3: return "t";
			default: throw new IllegalArgumentException();
		}
		
	}
	
	public double getArgument(int n) {
		switch (n) {
			case 0: return getN0();
			case 1: return getGrowthRate();
			case 2: return getGrowthRate2();
			case 3: return getTime1();
            case 4: return getPlateauTime();
			default: throw new IllegalArgumentException();
		}
	}
	
	public void setArgument(int n, double value) {
		switch (n) {
			case 0: setN0(value); break;
			case 1: setGrowthRate(value); break;
			case 2: setGrowthRate2(value); break;
			case 3: setTime1(value); break;
            case 4: setPlateauTime(value); break;
			default: throw new IllegalArgumentException();
		}
	}

	public DemographicFunction getCopy() {
		ExpConstExpDemographic df = new ExpConstExpDemographic(getUnits());
		df.setN0(getN0());
		df.setGrowthRate(getGrowthRate());
		df.r2 = r2;
		df.time1 = time1;
        df.plateauTime = plateauTime;

		return df;
	}

	//
	// private stuff
	//

	private double r2;
	private double time1;
    private double plateauTime;
}
