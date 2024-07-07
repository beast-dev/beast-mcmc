/*
 * PiecewisePopulationSizeModel.java
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

package dr.evomodel.coalescent;

import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * Implements a population size model that is a piecewise sequence of other population size models
 *
 * @author Andrew Rambaut
 */
public class PiecewisePopulationSizeModel extends PopulationSizeModel {

    //
    // Public stuff
    //

    /**
     * Construct piecewise model with from list of epochs. The N0 parameter represents the population size
     * at the present (time = 0) and then the epochs occur going back in time for the durations given by the
     * epochDurations parameter (with the last epoch stretching to infinity). The population size at the
     * start of the next epoch is given by the previous one. Thus these epoch population size models should
     * be created without their own N0 parameters.
     *
     * @param N0Parameter the population size at the present
     * @param epochs a list of population size models
     * @param epochDurations the duration of each epoch (for all but the last epoch)
     * @param units the units
     */
    public PiecewisePopulationSizeModel(Parameter N0Parameter,
                                        List<PopulationSizeModel> epochs,
                                        Parameter epochDurations,
                                        Type units) {

        this(/*PiecewisePopulationSizeModelParser.PIECEWISE_POPULATION_SIZE_MODEL*/"", N0Parameter, epochs, epochDurations, units);
    }

    /**
     * Construct piecewise model with from list of epochs
     * @param name the name of the model
     * @param N0Parameter the population size at the present
     * @param epochs a list of population size models
     * @param epochDurations the duration of each epoch (for all but the last epoch)
     * @param units the units
     */
    public PiecewisePopulationSizeModel(String name,
                                        Parameter N0Parameter,
                                        List<PopulationSizeModel> epochs,
                                        Parameter epochDurations,
                                        Type units) {

        super(name, N0Parameter, units);

        if (epochDurations.getDimension() != epochs.size() - 1) {
            throw new IllegalArgumentException("The number of epoch durations should be one fewer than the number of epochs");
        }

        this.epochs = new ArrayList<>();
        for (PopulationSizeModel epoch : epochs) {
            // listen for updates...
            addModel(epoch);

            // add the
            this.epochs.add(epoch);
        }

        this.epochDurations = epochDurations;
        addVariable(this.epochDurations);

        addStatistic(new TransitionTimesStatistic());
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        fireModelChanged();
        updateEpochs = true;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
        updateEpochs = true;
    }

    private void setupEpochs() {
        if (updateEpochs) {
            // update the N0 for each epoch from the popsize at the end of the preceding
            // epoch.
            double t = epochDurations.getParameterValue(0);
            double N = epochs.get(0).getPopulationSizeFunction().getLogDemographic(t);
            for (int i = 1; i < epochs.size(); i++) {
                epochs.get(i).setLogN0(N);
                N = epochs.get(i).getPopulationSizeFunction().getLogDemographic(t);
            }
            updateEpochs = false;
        }
    }

    /**
     * Return the population size function as an anonymous class...
     * @return the function
     */
    public PopulationSizeFunction getPopulationSizeFunction() {

        if (populationSizeFunction == null) {
            populationSizeFunction = new PopulationSizeFunction() {
                @Override
                public double getLogDemographic(double t) {
                    setupEpochs();

                    int epochIndex = 0;
                    double epochTime = t;
                    while (epochIndex < epochDurations.getDimension() && t > epochDurations.getParameterValue(epochIndex)) {
                        epochTime -= epochDurations.getParameterValue(epochIndex);
                        epochIndex += 1;
                    }

                    return epochs.get(epochIndex).getPopulationSizeFunction().getLogDemographic(epochTime);
                }

                @Override
                public double getIntegral(double startTime, double finishTime) {
                    setupEpochs();

                    double integral = 0.0;
                    double t0 = startTime;
                    double t1 = finishTime;

                    int epochIndex = 0;

                    // crawl until we find the epoch that contains the start time
                    while (epochIndex < epochDurations.getDimension() && t0 > epochDurations.getParameterValue(epochIndex)) {
                        t0 -= epochDurations.getParameterValue(epochIndex);
                        epochIndex += 1;
                    }

                    // now crawl until we find the epoch that contains the finish time, adding
                    // all the piece-wise integrals as we go
                    while (epochIndex < epochDurations.getDimension() && t1 > epochDurations.getParameterValue(epochIndex)) {
                        integral += epochs.get(epochIndex).getPopulationSizeFunction()
                                .getIntegral(t0, epochDurations.getParameterValue(epochIndex));
                        t0 = 0;
                        t1 -= epochDurations.getParameterValue(epochIndex);
                        epochIndex += 1;
                    }

                    // add the final piece of integral
                    integral += epochs.get(epochIndex).getPopulationSizeFunction().getIntegral(t0, t1);

                    return integral;
                }

                @Override
                public Type getUnits() {
                    return PiecewisePopulationSizeModel.this.getUnits();
                }

                @Override
                public void setUnits(Type units) {
                    throw new UnsupportedOperationException();
                }
            };
        }
        return populationSizeFunction;
    }

    public class TransitionTimesStatistic extends Statistic.Abstract {

        public TransitionTimesStatistic() {
            super("transitionTimes");
        }

        public int getDimension() {
            return epochDurations.getDimension();
        }

        public double getStatisticValue(int i) {
            double time = 0;
            for (int dim = 0; dim <= i; dim ++) {
                time += epochDurations.getStatisticValue(dim);
            }
            return time;
        }

    }

    //
    // protected stuff
    //

    private final List<PopulationSizeModel> epochs;
    private final Parameter epochDurations;

    private PopulationSizeFunction populationSizeFunction = null;

    private boolean updateEpochs = true;
}
