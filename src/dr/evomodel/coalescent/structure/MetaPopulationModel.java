/*
 * MetaPopulationModel.java
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

package dr.evomodel.coalescent.structure;

import dr.evolution.coalescent.structure.MetaPopulation;
import dr.evomodel.coalescent.DemographicModel;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;

/**
 * A wrapper for ConstantPopulation.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: MetaPopulationModel.java,v 1.3 2006/09/11 09:33:01 gerton Exp $
 */
public class MetaPopulationModel extends AbstractModel implements MetaPopulation {
    //
    // Public stuff
    //

    public static final String META_POPULATION_MODEL = "metaPopulationModel";

    /**
     * Construct demographic model with default settings
     */
    public MetaPopulationModel(ArrayList<DemographicModel> demographicModels, Parameter populationProportions) {

        this(META_POPULATION_MODEL, demographicModels, populationProportions);
    }

    /**
     * Construct demographic model with default settings
     */
    public MetaPopulationModel(String name, ArrayList<DemographicModel> demographicModelList, Parameter populationProportions) {

        super(name);

        this.populationProportions = populationProportions;

        if (populationProportions != null) {

            // Single demographic, each with a different weight
            addVariable(populationProportions);
            populationProportions.addBounds(new Parameter.DefaultBounds(1.0, 0.0, populationProportions.getDimension()));
            populationCount = populationProportions.getDimension() + 1;
            // Add copies of the demographicModel to the array
            for (int i = 1; i < populationCount; i++) {
                demographicModelList.add(demographicModelList.get(0));
            }
            demographicModels = demographicModelList.toArray(new DemographicModel[0]);
            addModel(demographicModels[0]);

        } else {

            // Several demographic models
            populationCount = demographicModelList.size();
            demographicModels = demographicModelList.toArray(new DemographicModel[0]);
            for (int i = 0; i < populationCount; i++) {
                addModel(demographicModels[i]);
            }

        }

        addStatistic(populationSizesStatistic);
    }

    // general functions

    private double getProportion(int population) {

        if (populationProportions == null) {
            return 1.0;
        }
        if (population > 0) {
            return populationProportions.getParameterValue(population - 1);
        }
        double proportion = 1.0;
        for (int i = 1; i < populationCount; i++) {
            proportion -= populationProportions.getParameterValue(i - 1);
        }
        return proportion;
    }


    public int getPopulationCount() {

        return populationCount;

    }

    public double[] getPopulationSizes(double time) {

        // make population size array
        double[] N = new double[populationCount];
        for (int i = 0; i < populationCount; i++) {
            N[i] = demographicModels[i].getDemographicFunction().getDemographic(time) * getProportion(i);
        }

        return N;
    }

    /* returns value of demographic function at time t  (population size; one entry of double[] getPopulationSizes)
    * (This function mirrors an equivalent function in DemographicFunction)
    */
    public double getDemographic(double time, int population) {

        return demographicModels[population].getDemographicFunction().getDemographic(time) * getProportion(population);

    }

    /* calculates the integral 1/N(x) dx from start to finish, for one of the populations
    * (This function mirrors an equivalent function in DemographicFunction)
    */
    public double getIntegral(double start, double finish, int population) {

        double integral = demographicModels[population].getDemographicFunction().getIntegral(start, finish);
        return integral / getProportion(population);

    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
    }

    protected void storeState() {
    } // no additional state needs storing

    protected void restoreState() {
        fireModelChanged();
    } // no additional state needs restoring

    protected void acceptState() {
    } // no additional state needs accepting

    public void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // nothing to do
    }

    private Statistic populationSizesStatistic = new Statistic.Abstract() {

        public String getStatisticName() {
            return "populationSizes";
        }

        public int getDimension() {
            return populationCount;
        }

        public double getStatisticValue(int dim) {

            double metaN0 = demographicModels[dim].getDemographicFunction().getDemographic(0);
            return metaN0 * getProportion(dim);
        }

    };

    /**
     * Parses an element from an DOM document into a ConstantPopulation.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MetaPopulationModel.META_POPULATION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            ArrayList<DemographicModel> demographics = new ArrayList<DemographicModel>(0);
            Parameter populationProportions = null;

            for (int i = 0; i < xo.getChildCount(); ++i) {
                final Object o = xo.getChild(i);
                if (o instanceof DemographicModel) {
                    demographics.add((DemographicModel) o);
                } else if (o instanceof Parameter) {
                    if (populationProportions != null) {
                        throw new Error("Allowed at most one Parameter in a MetaPopulationModel");
                    }
                    populationProportions = (Parameter) o;
                } else {
                    throw new Error("A MetaPopulationModel may only have children of type Parameter or DemographicModel");
                }
            }

            // Do sanity checking.
            if (populationProportions == null) {
                if (demographics.size() < 2) {
                    throw new Error("A MetaPopulationModel must have at least 2 DemographicModels (or a Parameter)");
                }
            } else {
                if (demographics.size() != 1) {
                    throw new Error("A MetaPopulationModel with a Parameter must have exactly one DemographicModel");
                }
            }

            return new MetaPopulationModel(demographics, populationProportions);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A model that represents a subdivided population.";
        }

        public Class getReturnType() {
            return MetaPopulationModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(DemographicModel.class, 1, 999),    // at least one required
                new ElementRule(Parameter.class, true)             // optional
        };
    };

    //
    // protected stuff
    //

    private int populationCount;
    private Parameter populationProportions;
    private DemographicModel[] demographicModels;
}
