/*
 * TwoColourSampler.java
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

import dr.evolution.tree.ColourChange;
import dr.evolution.tree.Tree;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

// This is probably redundant but is left in case it is needed as a comparison.


/**
 * @author Alexei Drummond
 * @author Gerton Lunter
 *
 * @version $Id: TwoColourSampler.java,v 1.2 2006/07/03 15:15:39 gerton Exp $
 */
//public class TwoColourSampler implements ColourSampler {
abstract class TwoColourSampler implements ColourSampler {

    // m[0] is the rate of migration from population 1 to population 0
    // m[1] is the rate of migration from population 0 to population 1

    /**
     */
    public TwoColourSampler() {
    }

    
    /**
     * @return the number of colours of this sampling model (always 2).
     */
    public int getColourCount() {
        return 2;
    }

    public TreeColouring sampleTreeColouring(Tree tree, ColourChangeMatrix colourChangeMatrix, double[] N) {
        return null;
    }


    /**
     * Pr(T>v|x,y,t)
     * The probability of surviving until time v without migrating,
     * conditional on starting in state x at time 0 and being in state y at time t.
     *
     * @param x state at time 0
     * @param y state at time t
     * @param v survival time
     * @return probability of surviving until time v without migrating,
     *         conditional on starting in state x at time 0 and being in state y at time t.
     */
    private final double conditionalSurvivalProbability(int x, double v, int y, double t, double[] m) {

        // the total rate of migration
        double mt = m[0]+m[1];

        if (v < 0 || v > t) {
            throw new IllegalArgumentException("v must be non-negative and not exceed t\n v="+v+" t="+t);
        } else {

            double a,b,c,d;

            if (x==y) {

                a = m[0]*m[1];
                b = -a;
                c = m[1-x];
                d = m[x];

            } else {

                a = m[1-y];
                b = m[1-x];
                c = 1;
                d = -1;

            }

            double exponent = (b*c-a*d)/(c*d*-mt);
            double C = Math.pow(c*Math.exp(-mt*t)+d, -exponent);
            return C * Math.exp(b*v/d) * Math.pow( d + c * Math.exp(-mt*(t-v)), exponent);

        }
    }

    

    /**
     * The inverse of the conditionalSurvivalProbability function. Returns the survival time until the next
     * migration conditional on starting at state x at time 0, being in state y at time t and having a
     * survival probability of w.
     * @param x state at time 0
     * @param y state at time t
     * @param t total elapsed time between x and y
     * @param m rate parameters
     * @return the survival time to the next migration
     * conditional on starting at state x at time 0, being in state y at time t and having a
     * survival probability of w.
     */
    private double randomConditionalSurvivalTime(int x, int y, double t, double[] m) throws NoEventException {

    	double U = MathUtils.nextDouble();
        double vLeft = 0.0;
        double vRight = t;
        double fLeft = 1.0-U;
    	double fRight = conditionalSurvivalProbability(x,t,y,t, m) - U;
    	
    	if (fRight >= 0.0) {
    		if (x != y) {
    			// No event but colours at ends are not identical
    			throw new IllegalArgumentException("Problem in randomConditionalSurvivalTime for x="+x+", y="+y+" t="+t+" U="+U+"fRight="+fRight);
    		}
    		throw new NoEventException();
    	}
    	
        // We are solving the equation "conditionalSurvivalprobability(x,v,y,t,m) - U == 0"
        // for v.  The solution is bound by [vLeft,vRight], with corresponding function
        // values fLeft and fRight.  The algorithm is basically Newton-Raphson, using an approximate
    	// derivative.  Two limits bracketing the solution are kept as well, and a fallback is
    	// used if the Newton iteration jumps out of the bracket.  This gives a robust and efficient
    	// root finder.
            	
    	double vNew = vLeft;   // the bracket position that has changed
    	double vOld = vRight;  // the previous value of vNew
    	double fNew = fLeft;
    	double fOld = fRight;

    	while (Math.abs(vNew-vOld)>1e-9) {
    		
    		// Calculate proposed new v, using the changed bracket value.  
    		// This is expected to give the highest performance, but may run out of the bracket
    		double vProp = vNew - fNew/((fOld-fNew)/(vOld-vNew));
    		
    		if ((vLeft >= vProp) || (vProp >= vRight)) {
    			// Outside of bracket, so use safe option
    			vProp = vLeft - fLeft/((fRight-fLeft)/(vRight-vLeft));
    		}

    		// New function value
			double fProp = conditionalSurvivalProbability(x,vProp,y,t, m) - U;

			vOld = vNew;
			fOld = fNew;
			vNew = vProp;
			fNew = fProp;
			
			// Update bracket
			if (fProp < 0.0) {
				vRight = vProp;
				fRight = fProp;
			} else {
				vLeft = vProp;
				fLeft = fProp;
			}
    	}
    	
    	return vNew;

    }

    /**
     * Samples a migration event on a two-coloured branch, conditional on colours at both ends/
     * Migration process is forwards in (natural) time, so we are going down the tree.
     * @param currentColour indeed
     * @param currentHeight indeed
     * @param childColour colour of the branch at the child node (below current)
     * @param childHeight height of the child node (below current, i.e. lower height)
     * @param m migration parameters
     * @return ColourChange event, with the colour referring to the branchlet *below* this point
     */
    public ColourChange randomConditionalMigrationEvent(int currentColour, double currentHeight,
                                                        int childColour, double childHeight, double[] m)
            throws NoEventException {

        // Draw a valid time (or throw NoEventException)
    	if (currentHeight < childHeight) {
    		throw new IllegalArgumentException("currentHeight "+currentHeight+" is below childHeight="+childHeight);
    	}
        double time = randomConditionalSurvivalTime(currentColour, childColour, currentHeight-childHeight, m);

        // Return the corresponding event
        return new ColourChange( currentHeight-time, 1-currentColour );

    }
    
    /**
     * Samples migration events on a two-coloured branch, conditional on colours at both ends
     * Migration process is forwards in (natural) time, so we are going down the tree.
     * Returns a list of events, ordered forward in time (i.e. colour refers to branch *below* the event)
     *      *
     * @param parentColour indeed
     * @param parentHeight indeed
     * @param childColour colour of the branch at the child node (below parent)
     * @param childHeight height of the child node (below parent, i.e. lower height)
     * @param m migration parameters
     * @return List of ColourChange events, with the colour referring to the branch *below* each event
     */
    public List<ColourChange> sampleConditionalMigrationEvents2(int parentColour, double parentHeight,
    		int childColour, double childHeight, double[] m) {
    
    	List<ColourChange> colourChanges = new ArrayList<ColourChange>();

    	if (parentHeight < childHeight) {
    		throw new IllegalArgumentException("sampleConditionalMigrationEvents: parentHeight="+parentHeight+" childHeight="+childHeight+", not good.");
    	}
    	
        // Sample migration events, going from the parent (current) to the child, forward in natural time
        // until a NoEventException breaks the loop.  Migration events are returned as ColourChange-s
        try {
        	
        	int currentColour = parentColour;
        	double currentHeight = parentHeight;
        	
            while (true) { 
                ColourChange nextEvent = 
                	randomConditionalMigrationEvent(currentColour, currentHeight, childColour, childHeight, m);

                // We abuse the ColorChange interface, which is supposed to record events going up the tree,
                // (in the coalescent direction), for the duration of this loop.  So, getColourAbove should 
                // really read 'getColourBelow'
                currentHeight = nextEvent.getTime();
                currentColour = nextEvent.getColourAbove();   // colour *below* current height
            	
                // record event
                colourChanges.add(nextEvent);
                
            }
        } catch (NoEventException nee) {
            // no more events
        }

        // Reverse the list
        reverseColourChangeList( colourChanges, parentColour);
        
    	return colourChanges;
    }

    
    /**
     * Samples migration events on a two-coloured branch, conditional on colours at both ends
     * Migration process is forwards in (natural) time, so we are going down the tree.
     * Returns a list of events, ordered forward in time (i.e. colour refers to branch *below* the event)
     * This version implements a rejection algorithm that can be easily extended to more than 2 colours
     *      *
     * @param parentColour indeed
     * @param parentHeight indeed
     * @param childColour colour of the branch at the child node (below parent)
     * @param childHeight height of the child node (below parent, i.e. lower height)
     * @param m migration parameters
     * @return List of ColourChange events, with the colour referring to the branch *below* each event
     */
     public List<ColourChange> sampleConditionalMigrationEvents(int parentColour, double parentHeight,
    		int childColour, double childHeight, double[] m) {
    
    	List<ColourChange> colourChanges = new ArrayList<ColourChange>();
    	int currentColour;
    	double currentHeight;
    	
    	// Reject until we get the child colour
    	do {
    		
    		colourChanges.clear();
    		currentColour = parentColour;
    		currentHeight = parentHeight;
    		
    		// Sample events until we reach the child
    		do { 

    			// Sample a waiting time
    			double totalRate = m[ 1-currentColour ];
        		double U = MathUtils.nextDouble();
        		
        		// Neat trick (Rasmus Nielsen): 
        		// If colours of parent and child differ, sample conditioning on at least 1 event
        		if ((parentColour != childColour) && (colourChanges.size() == 0)) {
        			
        			double minU = Math.exp( -totalRate * (parentHeight-childHeight) );
        			U = minU + U*(1.0-minU);
        			
        		}
        		
        		// Calculate the waiting time, and update currentHeight
        		double time = -Math.log( U )/totalRate;
        		currentHeight -= time;
        		
        		if (currentHeight > childHeight) {
               		// Not yet reached the child.  "Sample" an event
            		currentColour = 1 - currentColour;
            		// Add it to the list
            		colourChanges.add( new ColourChange( currentHeight, currentColour ) );
        		}
        		
    		} while (currentHeight > childHeight);
    		
    	} while (currentColour != childColour);

        // Reverse the list
        reverseColourChangeList( colourChanges, parentColour);
    	    	
    	return colourChanges;

    }

     
     class NoEventException extends Exception {

		/**
		 * 
		 */
		private static final long serialVersionUID = -6860022343065166240L;

     }
    /**
     * The probability density of being in state x in the time interval [0,t) and migrating to state y
     * at time t.
     *
     * In case x==y, the function returns the probability (not: density) of surviving in state x during
     * the interval [0,t].
     *
     * @param x state at time 0
     * @param y state at time t
     * @param t total elapsed time
     * @return probability density of being in state x in the time interval [0,t) and migrating to state y
     *         at time t.
     */
    public double migrationEventProposalDensity(int x, int y, double t, double[] m) {

        // Bugfix: Did not use correct interpretation of m[] array

        if (x==y) {
            return Math.exp(-t*m[1-x]);
        } else {
            return Math.exp(-t*m[1-x])*m[1-x];
        }
    }

    /**
     * @return the equilibrium distribution of this sampling model
     */
    public final double[] equilibrium(double[] m) {

        double mt = m[0]+m[1];
        return new double[] {m[0]/mt, m[1]/mt};
    }



    /**
     * Reverse the migration events in a list, to convert from the natural sampling ordering (parent-to-child)
     * to the natural coalescent ordering.  On input, getColourAbove returns the colour of the branch *below*
     * an event.
     */
     protected final void reverseColourChangeList( List<ColourChange> colourChanges, int parentColour ) {

    	Collections.reverse( colourChanges );
    	int colour;
    	for (int i=0; i < colourChanges.size(); i++) {
    		if (i<colourChanges.size()-1) {
    			colour = (colourChanges.get( i+1 )).getColourAbove();
    		} else {
    			colour = parentColour;
    		}
    		(colourChanges.get( i )).setColourAbove( colour );
    	}
    }

 /*
    public static void main(String[] args) {

         TwoColourSampler model = new TwoColourSampler();

         double[] m = new double[] {0.2, 0.5};
         int iterations = 500000;

         long time = System.currentTimeMillis();
         for (int j = 0; j < 4; j++) {
             double[] times = new double[iterations];
             for (int i = 0; i < iterations; i++) {
                 List events = model.sampleConditionalMigrationEvents2(j/2,1.0,j%2,0.0,m);
                 if (events.size()==0) {
                 	times[i] = 1.0;
                 } else {
                	times[i] = 1.0-((ColourChange)events.get(0)).getTime();
                 }
             }
             System.out.println("Analytic:"+(j/2) + "\t" + (j%2) + "\t" + DiscreteStatistics.mean(times) + "\t" + (DiscreteStatistics.stdev(times)/Math.sqrt(100000.0)));
         }
         System.out.println("time taken:" + (System.currentTimeMillis()-time) + "ms");
         
         time = System.currentTimeMillis();
         for (int j = 0; j < 4; j++) {
             double[] times = new double[iterations];
             for (int i = 0; i < iterations; i++) {
                 List events = model.sampleConditionalMigrationEvents(j/2,1.0,j%2,0.0,m);
                 if (events.size()==0) {
                 	times[i] = 1.0;
                 } else {
                 	times[i] = 1.0-((ColourChange)events.get(0)).getTime();
                 }
             }
             System.out.println("Rejection sampling:"+(j/2) + "\t" + (j%2) + "\t" + DiscreteStatistics.mean(times) + "\t" + (DiscreteStatistics.stdev(times)/Math.sqrt(100000.0)));
         }
         System.out.println("time taken:" + (System.currentTimeMillis()-time) + "ms");
     }
     
 */
 
}
