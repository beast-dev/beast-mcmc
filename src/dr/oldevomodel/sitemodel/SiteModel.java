/*
 * SiteModel.java
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

package dr.oldevomodel.sitemodel;

import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.SubstitutionModel;

/**
 * SiteModel - Specifies how rates and substitution models vary across sites.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SiteModel.java,v 1.77 2005/05/24 20:25:58 rambaut Exp $
 */

public interface SiteModel extends SiteRateModel {

    public static final String SITE_MODEL = "siteModel";

    /**
     * Get this site model's substitution model
     *
     * @return the substitution model
     */
    SubstitutionModel getSubstitutionModel();

    /**
     * Specifies whether SiteModel should integrate over the different categories at
     * each site. If true, the SiteModel will calculate the likelihood of each site
     * for each category. If false it will assume that there is each site can have a
     * different category.
     *
     * @return the boolean
     */
    boolean integrateAcrossCategories();  // TODO Consider moving into SiteRateModel

    /**
     * Get the category of a particular site. If integrateAcrossCategories is true.
     * then throws an IllegalArgumentException.
     *
     * @param site the index of the site
     * @return the index of the category
     */
    int getCategoryOfSite(int site);    // TODO Consider moving into SiteRateModel

    /**
     * Get the frequencyModel for this SiteModel.
     *
     * @return the frequencyModel.
     */
    FrequencyModel getFrequencyModel();
}