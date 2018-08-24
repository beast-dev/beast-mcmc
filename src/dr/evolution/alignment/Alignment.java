/*
 * Alignment.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.sequence.SequenceList;

/**
 * interface for any alignment data.
 *
 * @version $Id: Alignment.java,v 1.12 2005/05/24 20:25:55 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public interface Alignment extends SequenceList, SiteList
{

	/**
	 * Sets the dataType of this alignment. This can be different from
	 * the actual dataTypes of the sequences - they will be translated
	 * as required.
	 */
	void setDataType(DataType dataType);

	/**
	/**
	 * Returns string representation of single sequence in
	 * alignment with gap characters included.
	 */
	String getAlignedSequenceString(int sequenceIndex);

	/**
	 * Returns string representation of single sequence in
	 * alignment with gap characters excluded.
	 */
	String getUnalignedSequenceString(int sequenceIndex);

	public abstract class Abstract implements Alignment {
	    // **************************************************************
	    // PatternList IMPLEMENTATION
	    // **************************************************************

		/**
		 * @return number of patterns
		 */
		public int getPatternCount() {
			return getSiteCount();
		}

		/**
		 * @return number of invariant sites
		 */
		public int getInvariantCount() {
			throw new RuntimeException("Not implemented yet");
		}

		/**
		 * @return number of states for this siteList
		 */
		public int getStateCount() {
			return getDataType().getStateCount();
		}

		/**
		 * Gets the length of the pattern strings which will usually be the
		 * same as the number of taxa
		 * @return the length of patterns
		 */
		public int getPatternLength() {
			return getSequenceCount();
		}

		/**
		 * Gets the pattern as an array of state numbers (one per sequence) 
		 * @return the pattern at patternIndex
		 */
		public int[] getPattern(int patternIndex) {
			return getSitePattern(patternIndex);
		}

		@Override
		public double[][] getUncertainPattern(int patternIndex) {
			return getUncertainSitePattern(patternIndex);
		}

		/**
		 * @return state at (taxonIndex, patternIndex) 
		 */
		public int getPatternState(int taxonIndex, int patternIndex) {
			return getState(taxonIndex, patternIndex);
		}

		@Override
		public double[] getUncertainPatternState(int taxonIndex, int patternIndex) {
			return getUncertainPatternState(taxonIndex, patternIndex);
		}

		/** 
		 * Gets the weight of a site pattern (always 1.0)
		 */
		public double getPatternWeight(int patternIndex) {
			return 1.0;
		}

		/**
		 * @return the array of pattern weights
		 */
		public double[] getPatternWeights() {
			int count = getSiteCount();
			double[] weights = new double[count];
			for (int i = 0; i < count; i++)
				weights[i] = 1.0;
			return weights;
		}

		/**
		 * @return the frequency of each state
		 */
		public double[] getStateFrequencies() {
			return PatternList.Utils.empiricalStateFrequencies(this);
		}

		@Override
		public boolean areUnique() {
			return false;
		}

		// **************************************************************
	    // Identifiable IMPLEMENTATION
	    // **************************************************************

		protected String id = null;

		/**
		 * @return the id.
		 */
		public String getId() {
			return id;
		}

		/**
		 * Sets the id.
		 */
		public void setId(String id) {
			this.id = id;
		}

    }	
		
}
