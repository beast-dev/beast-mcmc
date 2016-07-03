/*
 * SampleStateAndCategoryModel.java
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
import dr.oldevomodel.substmodel.YangCodonModel;
import dr.oldevomodelxml.sitemodel.SampleStateAndCategoryModelParser;
import dr.inference.model.*;

import java.util.Vector;


/**
 * SampleStateAndCategoryModel - A SiteModel that has a discrete distribution of substitutionmodels over sites
 * designed for sampling of rate categories.
 *
 * @author Roald Forsberg
 */

public class SampleStateAndCategoryModel extends AbstractModel implements SiteModel, CategorySampleModel {

    public static final double OMEGA_MAX_VALUE = 100.0;
    public static final double OMEGA_MIN_VALUE = 0.0;


    /**
     * Constructor
     */
    public SampleStateAndCategoryModel(Parameter muParameter,
                                       Parameter categoriesParameter,
                                       Vector substitutionModels) {

        super(SampleStateAndCategoryModelParser.SAMPLE_STATE_AND_CATEGORY_MODEL);


        this.substitutionModels = substitutionModels;

        for (int i = 0; i < substitutionModels.size(); i++) {
            addModel((SubstitutionModel) substitutionModels.elementAt(i));

        }

        this.categoryCount = substitutionModels.size();
        sitesInCategory = new int[categoryCount];
        //	stateCount = ((SubstitutionModel)substitutionModels.elementAt(0)).getDataType().getStateCount();

        this.muParameter = muParameter;
        addVariable(muParameter);
        muParameter.addBounds(new Parameter.DefaultBounds(1000.0, 0.0, 1));

        this.categoriesParameter = categoriesParameter;
        addVariable(categoriesParameter);

        if (categoryCount > 1) {
            for (int i = 0; i < categoryCount; i++) {
                Parameter p = (Parameter)((YangCodonModel) substitutionModels.elementAt(i)).getVariable(0);
                Parameter lower = null;
                Parameter upper = null;

                if (i == 0) {
                    upper = (Parameter)((YangCodonModel) substitutionModels.elementAt(i + 1)).getVariable(0);
                    p.addBounds(new omegaBounds(lower, upper));
                } else {
                    if (i == (categoryCount - 1)) {
                        lower = (Parameter)((YangCodonModel) substitutionModels.elementAt(i - 1)).getVariable(0);
                        p.addBounds(new omegaBounds(lower, upper));
                    } else {
                        upper = (Parameter)((YangCodonModel) substitutionModels.elementAt(i + 1)).getVariable(0);
                        lower = (Parameter)((YangCodonModel) substitutionModels.elementAt(i - 1)).getVariable(0);
                        p.addBounds(new omegaBounds(lower, upper));
                    }
                }
            }
        }
    }

    // *****************************************************************
    // Interface SiteModel
    // *****************************************************************

    public SubstitutionModel getSubstitutionModel() {
        return null;
    }

    public boolean integrateAcrossCategories() {
        return false;
    }

    public int getCategoryCount() {
        return categoryCount;
    }

    public int getCategoryOfSite(int site) {
        return (int) categoriesParameter.getParameterValue(site);
    }

    public double getRateForCategory(int category) {
        throw new RuntimeException("getRateForCategory not available in this siteModel");
    }

    public double[] getCategoryRates() {
        throw new RuntimeException("getCategoryRates not available in this siteModel");
    }

    public double getSubstitutionsForCategory(int category, double time) {
        throw new RuntimeException("getSubstitutionsForCategory not available in this siteModel");
    }

    public void getTransitionProbabilities(double substitutions, double[] matrix) {
        throw new RuntimeException("getTransitionProbabilities not available in this siteModel");
    }

    /**
     * Get the frequencyModel for this SiteModel.
     *
     * @return the frequencyModel.
     */
    public FrequencyModel getFrequencyModel() {
        return ((SubstitutionModel) substitutionModels.elementAt(0)).getFrequencyModel();
    }

    // *****************************************************************
    // Interface CategorySampleModel
    // *****************************************************************

    /**
     * provide information to the categoriesParameter
     * about the number of sites
     */
    public void setCategoriesParameter(int siteCount) {
        categoriesParameter.setDimension(siteCount);
        categoriesParameter.addBounds(new Parameter.DefaultBounds(categoryCount, 0.0, siteCount));
        for (int i = 0; i < siteCount; i++) {

            int r = (int) (Math.random() * categoryCount);
            categoriesParameter.setParameterValue(i, r);
        }

        for (int j = 0; j < categoryCount; j++) {
            sitesInCategory[j] = 0;
        }

        for (int i = 0; i < siteCount; i++) {
            int value = (int) categoriesParameter.getParameterValue(i);
            sitesInCategory[value] = sitesInCategory[value] + 1;
        }
    }

    public void addSitesInCategoryCount(int category) {
        sitesInCategory[category] = sitesInCategory[category] + 1;
    }

    public void subtractSitesInCategoryCount(int category) {
        sitesInCategory[category] = sitesInCategory[category] - 1;
    }

    public int getSitesInCategoryCount(int category) {
        return sitesInCategory[category];
    }

    public void toggleRandomSite() {
    }


    /**
     * Get the expected proportion of sites in this category.
     *
     * @param category the category number
     * @return the proportion.
     */
    public double getProportionForCategory(int category) {
        throw new IllegalArgumentException("Not integrating across categories");
    }

    /**
     * Get an array of the expected proportion of sites in this category.
     *
     * @return an array of the proportion.
     */
    public double[] getCategoryProportions() {
        throw new IllegalArgumentException("Not integrating across categories");
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Substitution model has changed so fire model changed event
        listenerHelper.fireModelChanged(this, object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

        if (variable == categoriesParameter) // instructs TreeLikelihood to set update flag for this pattern
            listenerHelper.fireModelChanged(this, this, index);
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    } // no additional state needs accepting

    public String toString() {
        StringBuffer s = new StringBuffer();

        for (int i = 0; i < categoryCount; i++) {
            s.append(sitesInCategory[i] + "\t");
        }
        /*for(int i = 0; i < categoriesParameter.getDimension(); i++){
              t = (int)(categoriesParameter.getParameterValue(i));// get result as integer
              s.append(String.valueOf(t) + "\t");
          }*/

        return s.toString();
    }

    private class omegaBounds implements Bounds<Double> {

        private final Parameter lowerOmega, upperOmega;

        public omegaBounds(Parameter lowerOmega, Parameter upperOmega) {

            this.lowerOmega = lowerOmega;
            this.upperOmega = upperOmega;
        }

        public omegaBounds(Parameter nearestOmega, boolean isUpper) {

            if (isUpper) {
                lowerOmega = nearestOmega;
                upperOmega = null;
            } else {
                lowerOmega = null;
                upperOmega = nearestOmega;
            }
        }

        /**
         * @return the upper limit of this hypervolume in the given dimension.
         */
        public Double getUpperLimit(int dimension) {

            if (dimension != 0)
                throw new RuntimeException("omega parameters have wrong dimension " + dimension);

            if (upperOmega == null)
                return OMEGA_MAX_VALUE;
            else
                return upperOmega.getParameterValue(dimension);

        }

        /**
         * @return the lower limit of this hypervolume in the given dimension.
         */
        public Double getLowerLimit(int dimension) {

            if (dimension != 0)
                throw new RuntimeException("omega parameters have wrong dimension " + dimension);

            if (lowerOmega == null)
                return OMEGA_MIN_VALUE;

            else
                return lowerOmega.getParameterValue(dimension);
        }

        /**
         * @return the dimensionality of this hypervolume.
         */
        public int getBoundsDimension() {
            return 1;
        }

    }

    /**
     * mutation rate parameter
     */
    private final Parameter muParameter;

    private final int[] sitesInCategory;

    private final Parameter categoriesParameter;

    private final Vector substitutionModels;

    private final int categoryCount;

//	private int stateCount;
}
