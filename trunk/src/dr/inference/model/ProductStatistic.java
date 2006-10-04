/*
 * ProductStatistic.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.inference.model;

import dr.xml.*;

import java.util.Vector;

/**
 * @version $Id: ProductStatistic.java,v 1.2 2005/05/24 20:26:00 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class ProductStatistic extends Statistic.Abstract {
	
	public static String PRODUCT_STATISTIC = "product";

    private int dimension = 0;

	public ProductStatistic(String name) {
		super(name);
	}
	
	public void addStatistic(Statistic statistic) {
        if (dimension == 0) {
            dimension = statistic.getDimension();
        } else if (dimension != statistic.getDimension()) {
            throw new IllegalArgumentException();
        }
		statistics.add(statistic);
	}
	
	public int getDimension() { return dimension; }

	/** @return mean of contained statistics */
	public double getStatisticValue(int dim) {

        double product = 1.0;

		Statistic statistic;
		
		for (int i = 0; i < statistics.size(); i++) {
			statistic = (Statistic)statistics.get(i);
			product *= statistic.getStatisticValue(dim);
		}
		
		return product;
	}
		
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
		
		public String getParserName() { return PRODUCT_STATISTIC; }
		
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			ProductStatistic productStatistic = new ProductStatistic(PRODUCT_STATISTIC);
			
			for (int i =0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				if (child instanceof Statistic) {
					productStatistic.addStatistic((Statistic)child);
				} else {
					throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
				}
			}
				
			return productStatistic;
		}
		
		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		
		public String getParserDescription() {
			return "This element returns a statistic that is the mean of the child statistics.";
		}
		
		public Class getReturnType() { return ProductStatistic.class; }
		
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }
		
		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
		};		
	};
	

	// ****************************************************************
	// Private and protected stuff
	// ****************************************************************
	
	private Vector statistics = new Vector();
}
