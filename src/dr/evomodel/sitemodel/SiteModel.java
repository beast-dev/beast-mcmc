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
import dr.evomodel.substmodel.SubstitutionModel;
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

    public static final String SITE_MODEL = "siteModel";

    /**
     * Get this site model's substitution model
     * @return the substitution model
     */
    SubstitutionModel getSubstitutionModel();

    /**
	 * Specifies whether SiteModel should integrate over the different categories at
	 * each site. If true, the SiteModel will calculate the likelihood of each site
	 * for each category. If false it will assume that there is each site can have a
	 * different category.
     * @return the boolean
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
	 * Get the rate for a particular category. This will include the 'mu' parameter, an overall
     * scaling of the siteModel.
	 * @param category the category number
	 * @return the rate.
	 */
	double getRateForCategory(int category);

    /**
     * Get an array of the rates for all categories.
     * @return an array of rates.
     */
    double[] getCategoryRates();

	/**
	 * Get the expected proportion of sites in this category.
	 * @param category the category number
	 * @return the proportion.
	 */
	double getProportionForCategory(int category);

	/**
	 * Get an array of the expected proportion of sites for all categories.
	 * @return an array of proportions.
	 */
	double[] getCategoryProportions();

	/**
	 * Get the frequencyModel for this SiteModel.
	 * @return the frequencyModel.
	 */
	FrequencyModel getFrequencyModel();
}