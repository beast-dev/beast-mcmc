/*
 * DefaultBranchColouring.java
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

import java.util.ArrayList;
import java.util.List;


/**
 * BranchColouring.java
 * 
 * Represents a colouring on a single branch, both backward and forward in time.
 * 
 * @author Gerton Lunter
 * @author Andrew Rambaut
 *
 * @version $Id: BranchColouring.java,v 1.6 2006/08/12 12:55:44 gerton Exp $
 */

public class DefaultBranchColouring implements BranchColouring {

	public final static boolean checkSanity = true;
	
	public DefaultBranchColouring() {
		
		this.parentColour = -1;
		
	}
	
	public DefaultBranchColouring( int parentColour, int childColour ) {

		this.parentColour = parentColour;
		this.childColour = childColour;
		
	}
	
	public DefaultBranchColouring( DefaultBranchColouring history ) {
	
		parentColour = history.parentColour;
		childColour = history.childColour;
		colourChanges.addAll( history.colourChanges );
		changeHeights.addAll( history.changeHeights );
		finalSanityCheck();
		
	}
	
	public DefaultBranchColouring getCopy() {
		
		return new DefaultBranchColouring( this );
		
	}

	public void clear() {
		
		colourChanges.clear();
		changeHeights.clear();
	
	}

	public void reset( int parentColour, int childColour ) {
		
		this.parentColour = parentColour;
		this.childColour = childColour;
		colourChanges.clear();
		changeHeights.clear();
	
	}
	
	
	public int getParentColour() {
		return parentColour;
	}
	

	public int getChildColour() {
		return childColour;
	}

	
	/**
	 * Add a colour change event, going in the parent->child direction.
	 * 
	 * @param colour
	 * @param time
	 */
	public void addEvent( int colour, double time ) {
		
		colourChanges.add(colour);
		changeHeights.add(time);
		
		sanityCheck();
		
	}

	/**
	 * Add a number of change events below the current child
	 */
	public void addHistory( DefaultBranchColouring history ) {
		
		// require that my child colour is the new history's parent colour
		if (childColour != history.parentColour) {
			throw new Error("My child colour and the added parent colour don't match");
		}
		colourChanges.addAll( history.colourChanges );
		changeHeights.addAll( history.changeHeights );
		childColour = history.childColour;
		
		fullSanityCheck();

	}
	
	/**
	 * 
	 * @return number of events
	 */
	public int getNumEvents() {
		return colourChanges.size();
	}
	
	/**
	 * 
	 * @param i event, 0..getNumEvents inclusive.  
	 * @return Colour below event. If i==0, returns parent colour
	 */
	public int getForwardColourBelow( int i ) {
		finalSanityCheck();
		if (i==0) return parentColour;
		return colourChanges.get(i - 1);
	}

	/**
	 * 
	 * @param i event, 1..getNumEvents inclusive
	 * @return time of event
	 */
	public double getForwardTime( int i ) {
		finalSanityCheck();
		return changeHeights.get(i - 1);
	}

	/**
	 * 
	 * @param i event, 0..getNumEvents inclusive
	 * @return Colour above event.  If i==0, returns child colour
	 */
	public int getBackwardColourAbove( int i ) {
		finalSanityCheck();
		int total = colourChanges.size();
		if (i == total) return parentColour;
		return colourChanges.get(total - 1 - i);
	}

	/**
	 * 
	 * @param i event, 1..getNumEvents inclusive
	 * @return time of event
	 */
	public double getBackwardTime( int i ) {
		finalSanityCheck();
		int total = colourChanges.size();
		return changeHeights.get(total - i);
	}

	/**
	 * Returns event index corresponding to first event after time, going forward in time
	 * @param height
	 * @return 1 if first event is after time; getNumEvents+1 if time is after last event; 1 if no events exist
	 */
	public int getNextForwardEvent( double height ) {

		int i;
		for (i=0; i < changeHeights.size() && changeHeights.get(i) > height; i += 1 ) {}
		return i+1;

	}
	
    public double getTimeInColour(int colour, double parentHeight, double childHeight) {

        double totalTime = 0.0;
        int currentColour = getForwardColourBelow(0);   // parent colour
        double previousTime = parentHeight;
        for (int i = 1; i <= getNumEvents(); i++) {
            double currentTime = getForwardTime( i );
            if (currentColour == colour) {
                totalTime += previousTime - currentTime;
            }
            currentColour = getForwardColourBelow( i );
            previousTime = currentTime;
        }
        if (currentColour == colour) {
            totalTime += previousTime - childHeight;
        }
        return totalTime;

    }
	/**
	 * for backward compatibility (evolution.tree.ColouredTreePainter)
	 * 
	 * @return list of colour changes
	 */
	public List<ColourChange> getColourChanges() {

		List<ColourChange> cc = new ArrayList<ColourChange>(0);
		for (int i=1; i<=getNumEvents(); i++) {
	
			cc.add( new ColourChange( getBackwardTime(i), getBackwardColourAbove(i) ) );
	
		}
		
		return cc;

	}
	
	private void sanityCheck() {

        if (checkSanity) {
	        int total = colourChanges.size();
			if (total > 0) {
				if ( colourChanges.get(0) == parentColour) {
					throw new Error("First event does not change colour");
				}
			}
			if (total > 1) {
				if ((colourChanges.get(total-1)).intValue() == (colourChanges.get(total-2)).intValue()) {
						throw new Error("Last event does not change colour");
				}
				
				if ( changeHeights.get(total - 1) > changeHeights.get(total - 2) ) {
					throw new Error("Child event occurs before parent event");
				}
			}
        }
	}

	private void fullSanityCheck() {

        if (checkSanity) {
			sanityCheck();
	        int total = colourChanges.size();
			for (int i=1; i < total; i++) {
				if ((colourChanges.get(i)).intValue() == (colourChanges.get(i-1)).intValue()) {
						throw new Error("Event "+i+" does not change colour");
				}				
				if ( changeHeights.get(i) > changeHeights.get(i - 1) ) {
					throw new Error("Event "+i+" jumps back in time");
				}
			}
        }
	}

	
	private void finalSanityCheck() {
	
		if (checkSanity) {

			if (parentColour == -1) {
				throw new Error("Parent colour has not been set");
			}
			int total = colourChanges.size();
			if (total > 0) {
				if ( colourChanges.get(total - 1) != childColour) {
					throw new Error("Last event does not change colour into child's");
				}
			}
		}
	}

	private final List<Integer> colourChanges = new ArrayList<Integer>(0);   // new colours, forward in time
	private final List<Double> changeHeights = new ArrayList<Double>(0);   // change heights, forward in time
	private int parentColour;
	private int childColour;
	
}