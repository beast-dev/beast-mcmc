/*
 * DiscretizedSiteRateModel.java
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
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.evomodel.substmodel.SubstitutionModel;

import java.util.*;

/**
 * DiscretizedSiteRateModel - A SiteModel that has a discrete categories of rates across sites.
 *
 * @author Andrew Rambaut
 */

public class DiscretizedSiteRateModel extends AbstractModel implements SiteRateModel, Loggable {

    /**
     * Constructor for a rate homogenous (single category) SiteRateModel.
     */
    public DiscretizedSiteRateModel(String name) {
        this(name, null, 0.0, new HomogeneousRateDelegate(null));
    }

    /**
     * Constructor for a rate homogenous (single category) SiteRateModel.
     */
    public DiscretizedSiteRateModel(
            String name,
            Parameter nuParameter,
            double muWeight) {
        this(name, nuParameter, muWeight, new HomogeneousRateDelegate(null));
    }

    /**
     * Constructor for a discretized site rate model that uses a delegate to set
     * the category rates.
     */
    public DiscretizedSiteRateModel(
            String name,
            Parameter nuParameter,
            double muWeight,
            SiteRateDelegate delegate) {

        super(name);

        this.nuParameter = nuParameter;
        if (nuParameter != null) {
            addVariable(nuParameter);
            nuParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        }
        this.muWeight = muWeight;

        addStatistic(muStatistic);
        addStatistic(ratesStatistic);
        addStatistic(weightsStatistic);
        addStatistic(ratesTimesWeightsStatistic);

        this.delegate = delegate;
        addModel(delegate);

        categoryRates = new double[delegate.getCategoryCount()];
        categoryProportions = new double[delegate.getCategoryCount()];
        orderedRates = new double[delegate.getCategoryCount()][2]; // for storing ordered rate/weight pairs

        ratesKnown = false;
        orderedRatesKnown = false;
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

    // *****************************************************************
    // Interface SiteRateModel
    // *****************************************************************

    public int getCategoryCount() {
        return delegate.getCategoryCount();
    }

    public double[] getCategoryRates() {
        synchronized (this) {
            if (!ratesKnown) {
                calculateCategoryRates();
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

    private void calculateCategoryRates() {

        delegate.getCategories(categoryRates, categoryProportions);

        if (nuParameter != null) {
            double mu = getMu();
            for (int i = 0; i < getCategoryCount(); i++)
                categoryRates[i] *= mu;
        }

        ratesKnown = true;
        orderedRatesKnown = false;

    }

    private void calculateOrderedCategories() {
        for (int i = 0; i < categoryRates.length; i++) {
            orderedRates[i][0] = categoryRates[i];
            orderedRates[i][1] = categoryProportions[i];
        }
        Arrays.sort(orderedRates, Comparator.comparingDouble(a -> a[0]));
        orderedRatesKnown = true;
    }
    // *****************************************************************
    // Interface ModelComponent
    // *****************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // delegate has changed so fire model changed event
        ratesKnown = false;
        listenerHelper.fireModelChanged(this, object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        if (variable == nuParameter) {
        ratesKnown = false; // MAS: I changed this because the rate parameter can affect the categories if the parameter is in siteModel and not clockModel
//        } else {
//            throw new RuntimeException("Unknown variable in DiscretizedSiteRateModel.handleVariableChangedEvent");
//        }
        listenerHelper.fireModelChanged(this, variable, index);
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
        ratesKnown = false;
    }

    protected void acceptState() {
    } // no additional state needs accepting


    private final Statistic muStatistic = new Statistic.Abstract() {

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

    private final Statistic ratesStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "rates";
        }

        public int getDimension() {
            return getCategoryCount();
        }

        public double getStatisticValue(int dim) {
            if(!ratesKnown){
                calculateCategoryRates();
            }
            calculateOrderedCategories();

            return orderedRates[dim][0];
        }
    };

    private final Statistic weightsStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "weights";
        }

        public int getDimension() {
            return getCategoryCount();
        }

        public double getStatisticValue(int dim) {
            if(!ratesKnown){
                calculateCategoryRates();
            }
            calculateOrderedCategories();

            return orderedRates[dim][1];
        }

    };

    private final Statistic ratesTimesWeightsStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "rates_x_weights";
        }

        public int getDimension() {
            return getCategoryCount();
        }

        public double getStatisticValue(int dim) {
            if(!ratesKnown){
                calculateCategoryRates();
            }
            calculateOrderedCategories();

            return orderedRates[dim][1]* orderedRates[dim][0];
        }

    };

    /**
     * mutation rate parameter
     */
    private final Parameter nuParameter;

    private final double muWeight;

    private boolean ratesKnown, orderedRatesKnown;

    private final double[] categoryRates;
    private final double[][] orderedRates;

    private final double[] categoryProportions;

    private final SiteRateDelegate delegate;

    // This is here solely to allow the GammaSiteModelParser to pass on the substitution model to the
    // HomogenousBranchSubstitutionModel so that the XML will be compatible with older BEAST versions. To be removed
    // at some point.
    public SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    public void setSubstitutionModel(SubstitutionModel substitutionModel) {
        this.substitutionModel = substitutionModel;
    }

    private SubstitutionModel substitutionModel;

    @Override
    public LogColumn[] getColumns() {

        if (!ratesKnown) {
            calculateCategoryRates();
        }

        if (!orderedRatesKnown) {
            calculateOrderedCategories();
        }


        ArrayList<LogColumn> columns = new ArrayList<LogColumn>();

        // Columns sorted in order of increasing rate.


        Collections.addAll(columns, ratesStatistic.getColumns());
        Collections.addAll(columns, weightsStatistic.getColumns());
        Collections.addAll(columns, ratesTimesWeightsStatistic.getColumns());

        return columns.toArray(new LogColumn[columns.size()]);
    }
}