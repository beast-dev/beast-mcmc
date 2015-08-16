/*
 * NumberColumn.java
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

package dr.inference.loggers;

import java.text.DecimalFormat;

/**
 * An interface for a numerical column in a log.
 *
 * @version $Id: NumberColumn.java,v 1.4 2005/05/24 20:25:59 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
 
public abstract class NumberColumn extends LogColumn.Abstract {

	private int sf = -1;
	private int dp = -1;
	
	private double upperCutoff;
	private double[] cutoffTable;
	private final DecimalFormat decimalFormat = new DecimalFormat();
	private DecimalFormat scientificFormat = null;


	public NumberColumn(String label) {
        super(label);
        decimalFormat.setGroupingUsed(false); // not use comma
    }
	
	public NumberColumn(String label, int sf) {
		this(label);
		setSignificantFigures(sf);
	}
	
	/**
	 * Set the number of significant figures to display when formatted.
	 * Setting this overrides the decimal places option.
	 */
	public void setSignificantFigures(int sf) {
		this.sf = sf;
		this.dp = -1;
		
		upperCutoff = Math.pow(10,sf-1);
		cutoffTable = new double[sf];
		long num = 10;
		for (int i =0; i < cutoffTable.length; i++) {
			cutoffTable[i] = (double)num;
			num *= 10;
		}
		decimalFormat.setGroupingUsed(false);
		decimalFormat.setMinimumIntegerDigits(1);
		decimalFormat.setMaximumFractionDigits(sf-1);
		decimalFormat.setMinimumFractionDigits(sf-1);
		scientificFormat = new DecimalFormat(getPattern(sf));
	}
	
	/**
	 * Get the number of significant figures to display when formatted.
	 * Returns -1 if maximum s.f. are to be used.
	 */
	public int getSignificantFigures() { return sf; }
	
	/**
	 * Set the number of decimal places to display when formatted.
	 * Setting this overrides the significant figures option.
	 */
	public void setDecimalPlaces(int dp) {
		this.dp = dp;
		this.sf = -1;
	}
	
	/**
	 * Get the number of decimal places to display when formatted.
	 * Returns -1 if maximum d.p. are to be used.
	 */
	public int getDecimalPlaces() { return dp; }
	

    public String formatValue(double value) {
       if (dp < 0 && sf < 0) {
			// return it at full precision
			return Double.toString(value);
		}

		int numFractionDigits = 0;

		if (dp < 0) {
			
			double absValue = Math.abs(value);

			if ((absValue > upperCutoff) || (absValue < 0.1)) {

				return scientificFormat.format(value);

			} else {

				numFractionDigits = getNumFractionDigits(value);
			}

		} else {

			numFractionDigits = dp;
		}

		decimalFormat.setMaximumFractionDigits(numFractionDigits);
		decimalFormat.setMinimumFractionDigits(numFractionDigits);
		return decimalFormat.format(value);
    }

    /**
	 * Returns a string containing the current value for this column with
	 * appropriate formatting.
	 *
	 * @return the formatted string.
	 */
	protected String getFormattedValue() {
		return formatValue(getDoubleValue());
	}
	
	private int getNumFractionDigits(double value) {
		value = Math.abs(value);
		for (int i = 0; i < cutoffTable.length; i++) {
			if (value < cutoffTable[i]) return sf-i-1;
		}
		return sf - 1;
	}
	
	private String getPattern(int sf) {
		String pattern = "0.";
		 for (int i =0; i < sf-1; i++) {
		 	pattern += "#";
		 }
		 pattern += "E0";
		return pattern;
	}

	/**
	 * Returns the current value as a double.
	 */
	public abstract double getDoubleValue();
	
}
