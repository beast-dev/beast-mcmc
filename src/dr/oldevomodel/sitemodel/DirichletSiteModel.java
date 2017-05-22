/*
 * DirichletSiteModel.java
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
import dr.xml.*;

import java.util.logging.Logger;

/**
 * DirichletSiteModel - A SiteModel that has a free rate parameter for each category.
 *                      These rate parameters are constrained to sum to 1.
 *
 * @author Benjamin Redelings
 */

public class DirichletSiteModel extends AbstractModel implements SiteModel {

    public static final String SUBSTITUTION_MODEL = "substitutionModel";
    public static final String DIRICHLET_SITE_MODEL = "DirichletSiteModel";
    public static final String CATEGORIES = "categories";
    public static final String RATES = "rates";

    /**
     * Constructor for gamma+invar distributed sites. Either shapeParameter or
     * invarParameter (or both) can be null to turn off that feature.
     */
    public DirichletSiteModel(SubstitutionModel substitutionModel, int n_categories) {
        this(substitutionModel, new Parameter.Default(new double[n_categories]));
    }

    public DirichletSiteModel(SubstitutionModel substitutionModel, Parameter ratesParameter) {
        super(DIRICHLET_SITE_MODEL);

        // Set substitution model
        this.substitutionModel = substitutionModel;
        addModel(substitutionModel);

        // Allocate array for categories and proportions
        int n_categories = ratesParameter.getDimension();
        categoryRates = new double[n_categories];
        categoryProportions = new double[n_categories];

        // We haven't calculated any rates yet.
        ratesKnown = false;

        // Initialize the rate parameter value
        for(int i=0;i<n_categories;i++) {
            ratesParameter.setParameterValue(i, 1.0/n_categories);
            categoryProportions[i] = 1.0/n_categories;
        }

        // Add the rate parameters to the model
        this.ratesParameter = ratesParameter;
        addVariable(this.ratesParameter);
        this.ratesParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, ratesParameter.getDimension()));
    }

    public Parameter getRatesParameter() {
        return ratesParameter;
    }

//    public void setRatesParameter(Parameter parameter) {
//        removeVariable(ratesParameter);
//        ratesParameter = parameter;
//        addVariable(ratesParameter);
//    }

    // *****************************************************************
    // Interface SiteModel
    // *****************************************************************

    public boolean integrateAcrossCategories() {
        return true;
    }

    public int getCategoryCount() {
        return categoryRates.length;
    }

    public int getCategoryOfSite(int site) {
        throw new IllegalArgumentException("Integrating across categories");
    }

    public double getRateForCategory(int category) {
        calculateCategoryRates();

        return categoryRates[category];
    }

    public double[] getCategoryRates() {
        calculateCategoryRates();

        double[] rates = new double[categoryRates.length];

        for (int i = 0; i < rates.length; i++) {
            rates[i] = categoryRates[i];
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
        calculateCategoryRates();

        return categoryProportions[category];
    }

    /**
     * Get an array of the expected proportion of sites in this category.
     *
     * @return an array of the proportion.
     */
    public double[] getCategoryProportions() {
        calculateCategoryRates();

        return categoryProportions;
    }

    private double get_substitution_scale()
    {
        double scale=0;
        for(int i=0;i<categoryRates.length;i++) {
            //	    System.out.println(" get_substitution_scale: i = " + i);
            scale += categoryProportions[i]*ratesParameter.getParameterValue(i);
        }
        return scale;
    }

    /**
     * Calculate the category rates from the unscaled rate variables.
     */
    private void calculateCategoryRates()
    {
        synchronized (this) {
            if (ratesKnown)
                return;
        }

        double scale = get_substitution_scale();

        double temp = 1.0/scale;

        for(int i=0;i<categoryRates.length;i++)
            categoryRates[i] = ratesParameter.getParameterValue(i)*temp;

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

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == ratesParameter) {
            ratesKnown = false;
        } else {
            // This should not happen
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String[] getParserNames() {
            return new String[] {
                    getParserName(), "beast_"+getParserName()
            };
        }

        public String getParserName() {
            return DIRICHLET_SITE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException
        {
            SubstitutionModel substitutionModel = (SubstitutionModel) xo.getElementFirstChild(SUBSTITUTION_MODEL);

            String msg = "";

            Parameter ratesParameter = null;
            if (xo.hasChildNamed("rates")) {
                XMLObject cxo = (XMLObject) xo.getChild("rates");
                ratesParameter = (Parameter) cxo.getChild(Parameter.class);
                msg += "\n  with " + ratesParameter.getDimension() + " categories.";
            }


            if (msg.length() > 0) {
                Logger.getLogger("dr.evomodel").info("Creating site model: " + msg);
            } else {
                Logger.getLogger("dr.evomodel").info("Creating site model.");
            }

            return new DirichletSiteModel(substitutionModel, ratesParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A SiteModel that has Dirichlet distributed rates across sites";
        }

        public Class getReturnType() {
            return DirichletSiteModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(SUBSTITUTION_MODEL, new XMLSyntaxRule[]{
                        new ElementRule(SubstitutionModel.class)
                }),
                new ElementRule(RATES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                })
        };
    };

    /**
     * the substitution model for these sites
     */
    private SubstitutionModel substitutionModel = null;

    /**
     * the substitution rates: these sum to 1
     */
    private Parameter ratesParameter;

    private boolean ratesKnown;

    private double[] categoryRates;

    private double[] categoryProportions;
}
