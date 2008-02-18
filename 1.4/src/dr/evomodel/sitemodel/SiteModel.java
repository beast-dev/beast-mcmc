/*
 * SiteModel.java
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

package dr.evomodel.sitemodel; 

import dr.evomodel.substmodel.FrequencyModel;
import dr.inference.model.Model;

/**
 * SiteModel - Specifies how rates and substitution models vary across sites.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: SiteModel.java,v 1.77 2005/05/24 20:25:58 rambaut Exp $
 */

public interface SiteModel extends Model {

	/**
	 * Specifies whether SiteModel should integrate over the different categories at 
	 * each site. If true, the SiteModel will calculate the likelihood of each site
	 * for each category. If false it will assume that there is each site can have a
	 * different category.
	 */
	boolean integrateAcrossCategories();
	
	/**
	 * @return the number of categories of substitution processes
	 */
	int getCategoryCount();
	
	/**
	 * Get the category of a particular site. If integrateAcrossCategories is true.
	 * then throws an IllegalArgumentException.
	 * @param site the index of the site
	 * @return the index of the category
	 */
	int getCategoryOfSite(int site);
	
	/**
	 * Get the transition probability matrix for a particular category for a given time.
	 * @param category the category number
	 * @param time the branch length in units of time
	 * @param matrix an array of suitable size
	 */
	void getTransitionProbabilitiesForCategory(int category, double time, double[] matrix);
			
	/**
	 * Get the rate for a particular category.
	 * @param category the category number
	 * @return the rate.
	 */
	double getRateForCategory(int category);
			
	/**
	 * Get the number of substitutions per site for a particular category for a given time.
	 * @param category the category number
	 * @param time the branch length in units of time
	 * @return the expected number of substitutions per site.
	 */
	double getSubstitutionsForCategory(int category, double time);
			
	/**
	 * Get the transition probability matrix for a given number of substitutions per site.
	 * @param substitutions the branch length in expected number of substitutions per site
	 * @param matrix an array of suitable size
	 */
	void getTransitionProbabilities(double substitutions, double[] matrix);
						
	/**
	 * Get the expected proportion of sites in this category.
	 * @param category the category number
	 * @return the proportion.
	 */
	double getProportionForCategory(int category);
	
	/**
	 * Get an array of the expected proportion of sites in this category.
	 * @return an array of the proportion.
	 */
	double[] getCategoryProportions();
			
	/**
	 * Get the frequencyModel for this SiteModel.
	 * @return the frequencyModel.
	 */
	FrequencyModel getFrequencyModel();
			
};