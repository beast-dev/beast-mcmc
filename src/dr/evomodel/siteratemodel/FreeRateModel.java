/*
 * FreeRateModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.siteratemodel;

import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Matthew Hall
 *
 * The parameters of this version are the differences between rates, not the rates themselves. This allows the actual
 * rates to be in ascending order.
 */

public class FreeRateModel extends AbstractModel implements SiteRateModel, Citable {

    public FreeRateModel(
            String name,
            Parameter rateDifferences,
            Parameter weights) {

        super(name);

        this.rateDifferences = rateDifferences;
        this.weights = weights;
        this.dim = Math.min(rateDifferences.getDimension(), weights.getDimension());

        addVariable(rateDifferences);
        addVariable(weights);

        ratesKnown = false;
    }

    // *****************************************************************
    // Interface SiteRateModel
    // *****************************************************************

    @Override
    public int getCategoryCount() { return dim; }

    @Override
    public double[] getCategoryRates() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        return categoryRates;
    }

    @Override
    public double[] getCategoryProportions() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        return categoryProportions;
    }

    @Override
    public double getRateForCategory(int category) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        return categoryRates[category];
    }

    @Override
    public double getProportionForCategory(int category) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }
        return categoryProportions[category];
    }


    private void calculateCategoryRates() {

        double scale = 0.0;
        double sum = 0.0;
        double[] unnormalisedRates = new double[getCategoryCount()];

        unnormalisedRates[0] = rateDifferences.getParameterValue(0);

        for (int i = 1; i < dim; i++) {
            unnormalisedRates[i] = unnormalisedRates[i-1] + rateDifferences.getParameterValue(i);
        }

        if (categoryRates == null) {
            categoryRates = new double[dim];
        }

        if (categoryProportions == null) {
            categoryProportions = new double[dim];
        }

        for (int i = 0; i < dim; i++) {
            sum += weights.getParameterValue(i);
        }

        for (int i = 0; i < dim; ++i) {
            categoryProportions[i] = weights.getParameterValue(i) / sum;
        }

        for (int i = 0; i < dim; i++) {
            scale += categoryProportions[i] * unnormalisedRates[i];
        }

        for (int i = 0; i < dim; ++i) {
            categoryRates[i] = unnormalisedRates[i] / scale;
        }

//        double checker = 0;
//        for (int i = 0; i < dim; ++i) {
//            checker += categoryRates[i] * categoryProportions[i];
//        }

        ratesKnown = true;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Substitution model has changed so fire model changed event
        listenerHelper.fireModelChanged(this, object, index);
    }

    @Override
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == rateDifferences || variable == weights) {
            ratesKnown = false;
        } else {
        	throw new RuntimeException("Unknown variable in PdfSiteRateModel.handleVariableChangedEvent");
        }
        listenerHelper.fireModelChanged(this, variable, index);
    }

    @Override
    protected void storeState() { }

    @Override
    protected void restoreState() { ratesKnown = false; }

    @Override
    protected void acceptState() { }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Discrete probability distribution free rate heterogeneity model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public final static Citation CITATION = new Citation( // TODO Update
            new Author[]{
                    new Author("Z", "Yang")
            },
            "A space-time process model for the evolution of DNA Sequences",
            1995,
            "Genetics",
            139,
            993, 1005,
            Citation.Status.PUBLISHED
    );

    private final Parameter rateDifferences;
    private final Parameter weights;
    private final int dim;

    private double[] categoryRates;
    private double[] categoryProportions;
    private boolean ratesKnown;

    // This is here solely to allow the PdfSiteModelParser to pass on the substitution model to the
    // HomogenousBranchSubstitutionModel so that the XML will be compatible with older BEAST versions. To be removed
    // at some point.
    public SubstitutionModel getSubstitutionModel() {

        return substitutionModel;
    }

    public void setSubstitutionModel(SubstitutionModel substitutionModel) {
        this.substitutionModel = substitutionModel;
    }


    private SubstitutionModel substitutionModel;
}