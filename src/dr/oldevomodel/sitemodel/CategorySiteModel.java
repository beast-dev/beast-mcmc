/*
 * CategorySiteModel.java
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
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * CategorySiteModel - A SiteModel that has a different rate for each category;
 *
 * @author Alexei Drummond
 * @version $Id: CategorySiteModel.java,v 1.3 2004/10/01 22:40:04 alexei Exp $
 */

public class CategorySiteModel extends AbstractModel implements SiteModel {

    /**
     * @param rateParameter (relative to the rate of category 1)
     */
    public CategorySiteModel(SubstitutionModel substitutionModel,
                             Parameter muParameter,
                             Parameter rateParameter,
                             String categoryString,
                             String stateString,
                             int relativeTo) {

        super(SITE_MODEL);

        this.substitutionModel = substitutionModel;
        addModel(substitutionModel);

        this.muParameter = muParameter;
        addVariable(muParameter);
        muParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.rateParameter = rateParameter;
        addVariable(rateParameter);
        rateParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, Double.MIN_VALUE, rateParameter.getDimension()));

        states = stateString;

        if (states.length() != (rateParameter.getDimension() + 1)) {
            throw new IllegalArgumentException("States must have one more dimension than rate parameter!");
        }
        categoryCount = states.length();

        categories = new int[categoryString.length()];
        categoryWeights = new int[categoryCount];
        categoryRates = new double[categoryCount];
        for (int i = 0; i < categories.length; i++) {
            char state = categoryString.charAt(i);
            categories[i] = states.indexOf(state);
            categoryWeights[categories[i]] += 1;
        }

        siteCount = categories.length;
        this.relativeTo = relativeTo;
        ratesKnown = false;
    }

    /**
     * set mu
     */
    public void setMu(double mu) {
        muParameter.setParameterValue(0, mu);
    }

    /**
     * @return mu
     */
    public final double getMu() {
        return muParameter.getParameterValue(0);
    }


    public Parameter getMutationRateParameter() {
        return muParameter;
    }

    // *****************************************************************
    // Interface SiteModel
    // *****************************************************************

    public boolean integrateAcrossCategories() {
        return false;
    }

    public int getCategoryCount() {
        return categoryCount;
    }

    public int getCategoryOfSite(int site) {

        return categories[site];
    }

    public double getRateForCategory(int category) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        double mu = 1.0;
        if (muParameter != null) {
            mu = muParameter.getParameterValue(0);
        }

        return categoryRates[category] * mu;
    }

    public double[] getCategoryRates() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        double[] rates = new double[categoryRates.length];
        double mu = 1.0;
        if (muParameter != null) {
            mu = muParameter.getParameterValue(0);
        }
        for (int i = 0; i < rates.length; i++) {
            rates[i] = categoryRates[i] * mu;
        }

        return rates;
    }

    public void getTransitionProbabilities(double substitutions, double[] matrix) {
        substitutionModel.getTransitionProbabilities(substitutions, matrix);
    }

    /**
     * Get the expected proportion of sites in this category.
     *
     * @param category the category number
     * @return the proportion.
     */
    public double getProportionForCategory(int category) {
        throw new UnsupportedOperationException();
    }

    /**
     * Get an array of the expected proportion of sites in this category.
     *
     * @return an array of the proportion.
     */
    public double[] getCategoryProportions() {
        throw new UnsupportedOperationException();
    }

    /**
     *
     */
    private void calculateCategoryRates() {

        categoryRates[relativeTo] = 1.0;
        double total = 1.0;
        int count = 0;
        for (int i = 0; i < categoryRates.length; i++) {
            if (i != relativeTo) {
                categoryRates[i] = rateParameter.getParameterValue(count);
                total = categoryRates[i] * categoryWeights[i];
                count += 1;
            }
        }
        total /= (double) siteCount;

        // normalize so that total output rate is 1.0
        for (int i = 0; i < categoryRates.length; i++) {
            categoryRates[i] /= total;
        }

        ratesKnown = true;
    }

    /**
     * Get the frequencyModel for this SiteModel.
     *
     * @return the frequencyModel.
     */
    public FrequencyModel getFrequencyModel() {
        return substitutionModel.getFrequencyModel();
    }

    /**
     * Get the substitutionModel for this SiteModel.
     *
     * @return the substitutionModel.
     */
    public SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Substitution model has changed so fire model changed event
        listenerHelper.fireModelChanged(this, object, index);
    }

    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rateParameter) {
            ratesKnown = false;
        }
        listenerHelper.fireModelChanged(this, variable, index);
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
        ratesKnown = false;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    /**
     * the substitution model for these sites
     */
    private SubstitutionModel substitutionModel = null;

    /**
     * mutation rate parameter
     */
    private final Parameter muParameter;

    /**
     * rates parameter
     */
    private final Parameter rateParameter;

    private boolean ratesKnown;

    private final int categoryCount;

    private final double[] categoryRates;


    private final int[] categoryWeights;
    private final int[] categories;
    private final String states;
    private final int siteCount;
    private int relativeTo = 0;

}
