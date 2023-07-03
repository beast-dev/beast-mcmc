/*
 * FreeRateDelegate.java
 *
 * Copyright (c) 2002-2023 BEAST Developer Team
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

package dr.evomodel.siteratemodel;

import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.GeneralisedGaussLaguerreQuadrature;
import dr.math.distributions.GammaDistribution;
import dr.math.functionEval.GammaFunction;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

/**
 * FreeRateDelegate - A SiteModel delegate that implements the 'FreeRate' model.
 *
 * @author Andrew Rambaut
 * @version $Id: GammaSiteModel.java,v 1.31 2005/09/26 14:27:38 rambaut Exp $
 */

public class FreeRateDelegate extends AbstractModel implements SiteRateDelegate, Citable {



    /**
     * Constructor for gamma+invar distributed sites. Either shapeParameter or
     * invarParameter (or both) can be null to turn off that feature.
     */
    public FreeRateDelegate(
            String name,
            int categoryCount,
            Parameter rateParameter,
            Parameter weightParameter) {

        super(name);

        this.categoryCount = categoryCount;

        this.rateParameter = rateParameter;
        if (this.rateParameter.getDimension() == 1) {
            this.rateParameter.setDimension(categoryCount - 1);
        } else if (this.rateParameter.getDimension() != categoryCount - 1) {
            throw new IllegalArgumentException("Rate parameter should have have an initial dimension of one or category count - 1");
        }
        addVariable(this.rateParameter);
        this.rateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.weightParameter = weightParameter;
        if (this.weightParameter.getDimension() == 1) {
            this.weightParameter.setDimension(categoryCount - 1);
        } else if (this.weightParameter.getDimension() != categoryCount - 1) {
            throw new IllegalArgumentException("Weight parameter should have have an initial dimension of one or category count - 1");
        }

        addVariable(this.weightParameter);
        this.weightParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
    }

    // *****************************************************************
    // Interface SiteRateModel
    // *****************************************************************

    public int getCategoryCount() {
        return categoryCount;
    }

    public void getCategories(double[] categoryRates, double[] categoryProportions) {
        assert categoryRates != null && categoryRates.length == categoryCount;
        assert categoryProportions != null && categoryProportions.length == categoryCount;

        categoryRates[0] = 1.0;
        double sumRate = 1.0;
        double sumWeight = 0.0;
        for (int i = 0; i < categoryCount - 1; i++) {
            categoryRates[i + 1] = categoryRates[i] * rateParameter.getParameterValue(i);
            sumRate += categoryRates[i + 1];

            categoryProportions[i] = weightParameter.getParameterValue(i);
            sumWeight += categoryProportions[i];
        }
        // calculate the last value
        categoryProportions[categoryCount - 1] = 1.0 - sumWeight;

        // scale so their mean is 1
        for (int i = 0; i < categoryCount; i++) {
            categoryRates[i] = categoryCount * categoryRates[i] / sumRate;
        }
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        listenerHelper.fireModelChanged(this, object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        listenerHelper.fireModelChanged(this, variable, index);
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
    }

    protected void acceptState() {
    } // no additional state needs accepting


    /**
     * rate parameter
     */
    private final Parameter rateParameter;

    /**
     * weights parameter
     */
    private final Parameter weightParameter;


    private final int categoryCount;


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Discrete free-rate heterogeneity model";
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<>();
        return citations;
    }

    public final static Citation CITATION_YANG94 = new Citation(
            new Author[]{
                    new Author("", "")
            },
            "",
            1994,
            "J. Mol. Evol.",
            39,
            306, 314,
            Citation.Status.PUBLISHED
    );

    /**
     * Gives the category rates a mean of 1.0 and the proportions sum to 1.0
     * @param categoryRates
     * @param categoryProportions
     */
    public static void normalize(double[] categoryRates, double[] categoryProportions) {
        double mean = 0.0;
        double sum = 0.0;
        for (int i = 0; i < categoryRates.length; i++) {
            mean += categoryRates[i];
            sum += categoryProportions[i];
        }
        mean /= categoryRates.length;

        for(int i = 0; i < categoryRates.length; i++) {
            categoryRates[i] /= mean;
            categoryProportions[i] /= sum;
        }
    }


}