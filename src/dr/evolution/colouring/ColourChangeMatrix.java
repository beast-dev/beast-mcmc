/*
 * ColourChangeMatrix.java
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

package dr.evolution.colouring;

/**
 * MigrationMatrix.java
 * 
 * Represents a migration matrix, both backward and forward in time.
 * 
 * @author Gerton Lunter
 *
 * @version 1
 */


public class ColourChangeMatrix {

    public ColourChangeMatrix() { }
	
    public ColourChangeMatrix( double[] changeRates, int numColours ) {

        // For now, assume two colours
        if (numColours != 2) {
            throw new IllegalArgumentException("Only 2 colours supported");
        }

        this.numColours = numColours;
        this.bwMatrix = new double[numColours][numColours];
        this.equilibrium = new double[numColours];

        // Populate matrix
        bwMatrix[0][1] = changeRates[0];
        bwMatrix[1][0] = changeRates[1];
        calculateExitRates();

        // Calculate equilibrium distribution
        calculateEquilibrium();

    }

		/*
	 * populate diagonal of bw migration matrix
	 */
	private void calculateExitRates() {
	
		for (int i=0; i<numColours; i++) {
			double exitRate = 0.0;
			for (int j=0; j<numColours; j++) {
				if (i!=j) {
					exitRate += bwMatrix[i][j];
				}
			}
			bwMatrix[i][i] = -exitRate;
		}
		
	}
	
	private void calculateEquilibrium() {
		
		if (numColours == 2) {
			
			double r01 = bwMatrix[0][1];
			double r10 = bwMatrix[1][0];
			equilibrium[0] = r10/(r01+r10);
			equilibrium[1] = r01/(r01+r10);

		} else {
		
			throw new Error("Only 2 colours supported");
			
		}
	}
	
	
	/*
	 * Returns backward migration rate from i (child) to j (parent)
	 * 
	 */
	public double getBackwardRate(int i, int j) {
		return bwMatrix[i][j];
	}
	
	/*
	 * Returns forward migration rate from (parent) i to (child) j
	 */
	public double getForwardRate(int i, int j) {
		return equilibrium[j]*bwMatrix[j][i]/equilibrium[i];
	}

    /**
     * The probability of ending in state y after time t, conditional on starting in state x,
     * according to the forward evolution matrix.
     *
     * @param x starting state
     * @param y end state
     * @param t elapsed time
     * @return probability of starting from state x and ending in state y after time t.
     */
     public double forwardTimeEvolution(int x, int y, double t) {

    	if (t<0) {
    		throw new IllegalArgumentException("Cannot go backwards in time: t="+t);
    	}

    	double m[] = {getForwardRate(0,1), getForwardRate(1,0)};
        double mt = m[0]+m[1];

        if (y==x) {
            return (m[x]*Math.exp(-mt*t) + m[1-x])/mt;
        } else {
	    return m[x]*(1.0-Math.exp(-mt*t))/mt;
        }
    }

     
     /**
      * The probability of ending in state y after time t, conditional on starting in state x,
      * according to the backward evolution matrix.
      *
      * @param x starting state (child)
      * @param y end state (parent)
      * @param t elapsed time (backwards in time)
      * @return probability of starting from state x and ending in state y after time t.
      */
      public double backwardTimeEvolution(int x, int y, double t) {

     	if (t<0) {
     		throw new IllegalArgumentException("Cannot go backwards in time: t="+t);
     	}

     	double m[] = {getBackwardRate(0,1), getBackwardRate(1,0)};
        double mt = m[0]+m[1];

        if (y==x) {
	    //             return (m[1-y]*Math.exp(-mt*t) + m[y])/mt;
             return (m[x]*Math.exp(-mt*t) + m[1-x])/mt;
         } else {
	     //             return m[y]*(1.0-Math.exp(-mt*t))/mt;
             return m[x]*(1.0-Math.exp(-mt*t))/mt;
         }
     }

     
      /**
       * 
       * @return equilibrium distribution
       */
    public double[] getEquilibrium() {
    	
    	return equilibrium.clone();
    	
    }

    
    public double getEquilibrium(int i) {
    	
    	return equilibrium[i];
    	
    }

    	
	public static void main(String[] args) {
		
		double[] pars = {1.0, 0.5};
		ColourChangeMatrix mm = new ColourChangeMatrix( pars, 2 );
		
		System.out.println( "BW 0->0, matrix:"+mm.getBackwardRate(0,0));
		System.out.println( "BW 0->1, matrix:"+mm.getBackwardRate(0,1));
		System.out.println( "BW 1->0, matrix:"+mm.getBackwardRate(1,0));
		System.out.println( "BW 1->1, matrix:"+mm.getBackwardRate(1,1));
		System.out.println( "FW 0->0, matrix:"+mm.getForwardRate(0,0));
		System.out.println( "FW 0->1, matrix:"+mm.getForwardRate(0,1));
		System.out.println( "FW 1->0, matrix:"+mm.getForwardRate(1,0));
		System.out.println( "FW 1->1, matrix:"+mm.getForwardRate(1,1));
		System.out.println( "equilibrium="+mm.equilibrium[0]+","+mm.equilibrium[1] );
		System.out.println( "BW 0->0, t=1:"+mm.backwardTimeEvolution(0,0,1.0));
		System.out.println( "BW 0->1, t=1:"+mm.backwardTimeEvolution(0,1,1.0));
		System.out.println( "BW 1->0, t=1:"+mm.backwardTimeEvolution(1,0,1.0));
		System.out.println( "BW 1->1, t=1:"+mm.backwardTimeEvolution(1,1,1.0));
		System.out.println( "FW 0->0, t=1:"+mm.forwardTimeEvolution(0,0,1.0));
		System.out.println( "FW 0->1, t=1:"+mm.forwardTimeEvolution(0,1,1.0));
		System.out.println( "FW 1->0, t=1:"+mm.forwardTimeEvolution(1,0,1.0));
		System.out.println( "FW 1->1, t=1:"+mm.forwardTimeEvolution(1,1,1.0));

		System.out.println( "BW 0->0, t=infty:"+mm.backwardTimeEvolution(0,0,1000.0));
		System.out.println( "BW 0->1, t=infty:"+mm.backwardTimeEvolution(0,1,1000.0));
		System.out.println( "BW 1->0, t=infty:"+mm.backwardTimeEvolution(1,0,1000.0));
		System.out.println( "BW 1->1, t=infty:"+mm.backwardTimeEvolution(1,1,1000.0));
		System.out.println( "FW 0->0, t=infty:"+mm.forwardTimeEvolution(0,0,1000.0));
		System.out.println( "FW 0->1, t=infty:"+mm.forwardTimeEvolution(0,1,1000.0));
		System.out.println( "FW 1->0, t=infty:"+mm.forwardTimeEvolution(1,0,1000.0));
		System.out.println( "FW 1->1, t=infty:"+mm.forwardTimeEvolution(1,1,1000.0));
		
		
	}
	
	// private stuff
	private int numColours;
	private double[][] bwMatrix;
	private double[] equilibrium;
}


