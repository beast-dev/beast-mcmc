/*
 * GammaSiteRateDelegate.java
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
import dr.inference.model.*;
import dr.math.GeneralisedGaussLaguerreQuadrature;
import dr.math.distributions.GammaDistribution;
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
 * @version $Id: GammaSiteModel.java,v 1.31 2005/09/26 14:27:38 rambaut Exp $
 */

public class GammaSiteRateDelegate extends AbstractModel implements SiteRateDelegate, Citable {

    public static final DiscretizationType DEFAULT_DISCRETIZATION = DiscretizationType.QUADRATURE;

    public enum DiscretizationType {
        EQUAL,
        QUADRATURE
    };


    /**
     * Constructor for gamma+invar distributed sites. Either shapeParameter or
     * invarParameter (or both) can be null to turn off that feature.
     */
    public GammaSiteRateDelegate(
            String name,
            Parameter shapeParameter, int gammaCategoryCount,
            DiscretizationType discretizationType,
            Parameter invarParameter) {

        super(name);
        
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

        this.discretizationType = discretizationType;
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

        int offset = 0;

        if (invarParameter != null) {
            categoryRates[0] = 0.0;
            categoryProportions[0] = invarParameter.getParameterValue(0);
            offset = 1;
        }

        if (shapeParameter != null) {
            double alpha = shapeParameter.getParameterValue(0);
            final int gammaCatCount = categoryCount - offset;

            if (discretizationType == DiscretizationType.QUADRATURE) {
                setQuadratureRates(categoryRates, categoryProportions, alpha, gammaCatCount, offset);
            } else {
                setEqualRates(categoryRates, categoryProportions, alpha, gammaCatCount, offset);
            }
        } else if (offset > 0) {
            // just the invariant rate and variant rate
            categoryRates[offset] = 2.0;
            categoryProportions[offset] = 1.0 - categoryProportions[0];
        } else {
            categoryRates[0] = 1.0;
            categoryProportions[0] = 1.0;
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
     * shape parameter
     */
    private Parameter shapeParameter;

    /**
     * invariant sites parameter
     */
    private Parameter invarParameter;

    private DiscretizationType discretizationType;

    private int categoryCount;
    

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
            if (discretizationType == DiscretizationType.QUADRATURE) {
                citations.add(CITATION_FELSENSTEIN01);
            }
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

    public final static Citation CITATION_FELSENSTEIN01 = new Citation(
            new Author[]{
                    new Author("J", "Felsenstein")
            },
            "Taking Variation of Evolutionary Rates Between Sites into Account in Inferring Phylogenies",
            2001,
            "J. Mol. Evol.",
            53,
            447, 455,
            Citation.Status.PUBLISHED
    );

    private SubstitutionModel substitutionModel;

    private static GeneralisedGaussLaguerreQuadrature quadrature = null;

    /**
     * Set the rates and proportions using a Gauss-Laguerre Quadrature, as proposed by Felsenstein 2001, JME
     *
     * @param categoryRates
     * @param categoryProportions
     * @param alpha
     * @param catCount
     * @param offset
     */
    public static void setQuadratureRates(double[] categoryRates, double[] categoryProportions, double alpha, int catCount, int offset) {
        if (quadrature == null) {
            quadrature = new GeneralisedGaussLaguerreQuadrature(catCount);
        }
        quadrature.setAlpha(alpha - 1);

        double[] abscissae = quadrature.getAbscissae();
        double[] coefficients = quadrature.getCoefficients();

        for (int i = 0; i < catCount; i++) {
            categoryRates[i + offset] = abscissae[i] / alpha;
            categoryProportions[i + offset] = coefficients[i]/GammaFunction.gamma(alpha);
        }

    }

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
        double sumRates = 0.0;
        double sumProps = 0.0;
        System.out.println();
        System.out.println("Equal, alpha = 1.0");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        setQuadratureRates(categoryRates, categoryProportions, 1.0, catCount, 0);
        sumRates = 0.0;
        sumProps = 0.0;
        System.out.println();
        System.out.println("Quadrature, alpha = 1.0");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        // Table 3 from Felsenstein 2001, JME
        // Rates and probabilities chosen by the quadrature method for six rates and coefficient of
        // variation of rates among sites 1 (a 4 1)
        // Rate     Probability
        // 0.264    0.278
        // 0.898    0.494
        // 1.938    0.203
        // 3.459    0.025
        // 5.617    0.00076
        // 8.823    0.000003

        // Output
        // Quadrature, alpha = 1.0
        // cat	rate	proportion
        // 0	0.26383406085556455	0.27765014202987454
        // 1	0.8981499048217043	0.49391058305035496
        // 2	1.938320760238456	0.20300429674372977
        // 3	3.459408283352361	0.02466882036918974
        // 4	5.617305214541558	7.6304276746353E-4
        // 5	8.822981776190357	3.1150393875275343E-6

        setEqualRates(categoryRates, categoryProportions, 0.1, catCount, 0);
        sumRates = 0.0;
        sumProps = 0.0;
        System.out.println();
        System.out.println("Equal, alpha = 0.1");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        setQuadratureRates(categoryRates, categoryProportions, 0.1, catCount, 0);
        System.out.println();
        System.out.println("Quadrature, alpha = 0.1");
        System.out.println("cat\trate\tproportion");
        sumRates = 0.0;
        sumProps = 0.0;
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        setEqualRates(categoryRates, categoryProportions, 10.0, catCount, 0);
        System.out.println();
        System.out.println("Equal, alpha = 10.0");
        System.out.println("cat\trate\tproportion");
        sumRates = 0.0;
        sumProps = 0.0;
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        setQuadratureRates(categoryRates, categoryProportions, 10.0, catCount, 0);
        sumRates = 0.0;
        sumProps = 0.0;
        System.out.println();
        System.out.println("Quadrature, alpha = 10.0");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);

        setQuadratureRates(categoryRates, categoryProportions, 100.0, catCount, 0);
        sumRates = 0.0;
        sumProps = 0.0;
        System.out.println();
        System.out.println("Quadrature, alpha = 100.0");
        System.out.println("cat\trate\tproportion");
        for (int i = 0; i < catCount; i++) {
            System.out.println(i + "\t"+ categoryRates[i] +"\t" + categoryProportions[i]);
            sumRates += categoryRates[i];
            sumProps += categoryProportions[i];
        }
        System.out.println("SUM\t"+ sumRates +"\t" + sumProps);
    }
}