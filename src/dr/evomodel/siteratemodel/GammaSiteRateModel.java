/*
 * GammaSiteRateModel.java
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

import dr.inference.model.*;
import dr.math.GeneralisedGaussLaguerreQuadrature;
import dr.math.distributions.GammaDistribution;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.math.functionEval.GammaFunction;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

/**
 * GammaSiteModel - A SiteModel that has a gamma distributed rates across sites.
 *
 * @author Andrew Rambaut
 */

public class GammaSiteRateModel extends AbstractModel implements SiteRateModel, Citable {

    public GammaSiteRateModel(String name) {
        this(name,
                null,
                1.0,
                null,
                0,
                null);
    }

    public GammaSiteRateModel(String name, double alpha, int categoryCount) {
        this(name,
                null,
                1.0,
                new Parameter.Default(alpha),
                categoryCount,
                null);
    }

    public GammaSiteRateModel(String name, double alpha, int categoryCount, double pInvar) {
        this(name,
                null,
                1.0,
                new Parameter.Default(alpha),
                categoryCount,
                new Parameter.Default(pInvar));
    }

    public GammaSiteRateModel(String name, Parameter nu) {
        this(name,
                nu,
                1.0,
                null,
                -1,
                null);
    }

    /**
     * Constructor for gamma+invar distributed sites. Either shapeParameter or
     * invarParameter (or both) can be null to turn off that feature.
     */
    public GammaSiteRateModel(
            String name,
            Parameter nuParameter,
            double muWeight,
            Parameter shapeParameter, int gammaCategoryCount,
            Parameter invarParameter) {

        super(name);

        this.nuParameter = nuParameter;
        if (nuParameter != null) {
            addVariable(nuParameter);
            nuParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
        this.muWeight = muWeight;

        addStatistic(muStatistic);

        this.shapeParameter = shapeParameter;
        if (shapeParameter != null) {
            this.categoryCount = gammaCategoryCount;
            addVariable(shapeParameter);
//            shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 1E-3, 1));
            // removing the bounds on the alpha parameter - to make the prior more explicit
            shapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        } else {
            this.categoryCount = 1;
        }

        this.invarParameter = invarParameter;
        if (invarParameter != null) {
            this.categoryCount += 1;

            addVariable(invarParameter);
            invarParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));
        }

        categoryRates = new double[this.categoryCount];
        categoryProportions = new double[this.categoryCount];

        ratesKnown = false;
    }

    /**
     * set mu
     */
    public void setMu(double mu) {
        nuParameter.setParameterValue(0, mu / muWeight);
    }

    /**
     * @return mu
     */
    public final double getMu() {
        return nuParameter.getParameterValue(0) * muWeight;
    }

    /**
     * set alpha
     */
    public void setAlpha(double alpha) {
        shapeParameter.setParameterValue(0, alpha);
        ratesKnown = false;
    }

    /**
     * @return alpha
     */
    public final double getAlpha() {
        return shapeParameter.getParameterValue(0);
    }

    public void setRelativeRateParameter(Parameter nu) {
        this.nuParameter = nu;
    }

    // *****************************************************************
    // Interface SiteRateModel
    // *****************************************************************

    public int getCategoryCount() {
        return categoryCount;
    }

    public double[] getCategoryRates() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }

        for (int i = (invarParameter != null ? 1 : 0); i < categoryRates.length; i++) {
            // If a gamma rate is zero then the quantitization has failed numerically so return null.
            // This allows the likelihood to return -Inf and reject this state.
            if (categoryRates[i] == 0.0) {
                return null;
            }
        }

        return categoryRates;
    }

    public double[] getCategoryProportions() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }

        return categoryProportions;
    }

    public double getRateForCategory(int category) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }

        return categoryRates[category];
    }

    public double getProportionForCategory(int category) {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
            }
        }

        return categoryProportions[category];
    }

    /**
     * discretization of gamma distribution with equal proportions in each
     * category
     */
    private void calculateCategoryRates() {

        int offset = 0;

        if (invarParameter != null) {
            categoryRates[0] = 0.0;
            categoryProportions[0] = invarParameter.getParameterValue(0);
            offset = 1;
        }

        if (shapeParameter != null) {
            double alpha = shapeParameter.getParameterValue(0);
            final int gammaCatCount = categoryCount - offset;

            setEqualRates(categoryRates, categoryProportions, alpha, gammaCatCount, offset);
        } else if (offset > 0) {
            // just the invariant rate and variant rate
            categoryRates[offset] = 2.0;
            categoryProportions[offset] = 1.0 - categoryProportions[0];
        } else {
            categoryRates[0] = 1.0;
            categoryProportions[0] = 1.0;
        }

        if (nuParameter != null) {
            double mu = getMu();
            for (int i = 0; i < categoryCount; i++)
                categoryRates[i] *= mu;
        }

        ratesKnown = true;
    }

    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Substitution model has changed so fire model changed event
        listenerHelper.fireModelChanged(this, object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == shapeParameter) {
            ratesKnown = false;
        } else if (variable == invarParameter) {
            ratesKnown = false;
        } else if (variable == nuParameter) {
            ratesKnown = false; // MAS: I changed this because the rate parameter can affect the categories if the parameter is in siteModel and not clockModel
        } else {
            throw new RuntimeException("Unknown variable in GammaSiteRateModelParser.handleVariableChangedEvent");
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


    private Statistic muStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "mu";
        }

        public int getDimension() {
            return 1;
        }

        public String getDimensionName(int dim) {
            return getId();
        }

        public double getStatisticValue(int dim) {
            return getMu();
        }

    };


    /**
     * mutation rate parameter
     */
    private Parameter nuParameter;

    private double muWeight;

    /**
     * shape parameter
     */
    private Parameter shapeParameter;

    /**
     * invariant sites parameter
     */
    private Parameter invarParameter;

    private boolean ratesKnown;

    private int categoryCount;

    private double[] categoryRates;

    private double[] categoryProportions;


    // This is here solely to allow the GammaSiteModelParser to pass on the substitution model to the
    // HomogenousBranchSubstitutionModel so that the XML will be compatible with older BEAST versions. To be removed
    // at some point.
    public SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    public void setSubstitutionModel(SubstitutionModel substitutionModel) {
        this.substitutionModel = substitutionModel;
    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Discrete gamma-distributed rate heterogeneity model";
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<>();
        if (shapeParameter != null) {
            citations.add(CITATION_YANG94);
        }
        return citations;
    }

    public final static Citation CITATION_YANG94 = new Citation(
            new Author[]{
                    new Author("Z", "Yang")
            },
            "Maximum likelihood phylogenetic estimation from DNA sequences with variable rates over sites: approximate methods",
            1994,
            "J. Mol. Evol.",
            39,
            306, 314,
            Citation.Status.PUBLISHED
    );

    private SubstitutionModel substitutionModel;

    /**
     * set the rates as equally spaced quantiles represented by the mean as proposed by Yang 1994
     * @param categoryRates
     * @param categoryProportions
     * @param alpha
     * @param catCount
     * @param offset
     */
    public static void setEqualRates(double[] categoryRates, double[] categoryProportions, double alpha, int catCount, int offset) {
        for (int i = 0; i < catCount; i++) {
            categoryRates[i + offset] = GammaDistribution.quantile((2.0 * i + 1.0) / (2.0 * catCount), alpha, 1.0 / alpha);
            categoryProportions[i + offset] = 1.0;
        }

        normalize(categoryRates, categoryProportions);
    }

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

    public static void main(String[] argv) {
        final int catCount = 6;

        double[] categoryRates = new double[catCount];
        double[] categoryProportions = new double[catCount];

        setEqualRates(categoryRates, categoryProportions, 1.0, catCount, 0);
        System.out.println();
        System.out.println("Equal, alpha = 1.0");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
        }

    }
}