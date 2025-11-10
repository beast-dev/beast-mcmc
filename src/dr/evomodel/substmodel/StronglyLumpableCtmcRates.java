/*
 * LogCtmcRatesFromMatrixMatrixProductParameter.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.loggers.LogColumn;
import dr.inference.model.*;
import dr.math.matrixAlgebra.WrappedVector;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * @author Xinghua Tao
 * @author Marc Suchard
 */

public class StronglyLumpableCtmcRates extends AbstractModel implements LogAdditiveCtmcRateProvider.DataAugmented, Reportable {

    // Stores SuperInfo objects, which hold transition rate information between lumped states.
    private final SuperInfo[] map;
    private final List<Lump> lumps;
    //  Represents transition rates between different lumps.
    private final Parameter acrossRates;
//    private final DataType dataType;
    // Total number of states
    private final int stateCount;
    // Total number of possible state transitions (stateCount * (stateCount - 1))
    private final int dim;
    // Stores the computed CTMC rates
    private double[] rates;
    // A backup for rates (used for caching)
    private double[] storedRates;
    //A flag to check if rates are already computed
    private boolean ratesKnown = false;
    // A flag to enable caching (currently disabled)
    private static final boolean CACHE = false;

    public StronglyLumpableCtmcRates(String name,
                                     List<Lump> lumps,
                                     Parameter acrossRates,
                                     DataType dataType) {

        super(name);

        this.lumps = lumps;
        this.acrossRates = acrossRates;
//        this.dataType = dataType;
        this.stateCount = dataType.getStateCount();
        this.dim = stateCount * (stateCount - 1);
        this.rates = new double[dim];
        this.map = buildSuperMap();

        //Handling Lump Parameters
        for (Lump lump : lumps) {
            if (lump.rates != null) {
                addVariable(lump.rates);
            }

            if (lump.proportions != null) {
                for (Proportion p : lump.proportions) {
                    addVariable(p.parameter);
                }
//                addVariable(lump.proportions);
            }
        }

        addVariable(acrossRates);

        computeRates(rates); // TODO remove
    }

    //Computing Rates
//    Uses SuperInfo to determine if the transition is within the same lump.
//    If within, retrieves the rate from the lump.
//    If not within, assigns a rate of 0 (no transition).
    private double getRate(int index) {
//        SuperInfo info = map[index];
//        final double rate;
//        if (info.within) {
//            rate = info.withinRateParameter.getParameterValue(info.withinRateIndex);
//        } else {
//            rate = info.acrossProportionParameter.getParameterValue(info.acrossProportionIndex) *
//                    acrossRates.getParameterValue(info.acrossRateIndex);
//        }
//
//        return rate;
        return map[index].getRate();
    }

    // Compute all rates
    private void computeRates(double[] rates) {
        for (int i = 0; i < dim; ++i) {
            rates[i] = getRate(i);
        }
    }

    private SuperInfo[] buildSuperMap() {

        SuperInfo[] info = new SuperInfo[dim];

        int offset = 0;
    //    info[offset] = new SuperInfo(getLumpIndex(1),getLumpIndex(4));
        for (int i = 0; i < stateCount; ++i) {
            LumpIndex ii = getLumpIndex(i);
            assert ii != null;
            for (int j = i + 1; j < stateCount; ++j) {
                LumpIndex jj = getLumpIndex(j);
                assert jj != null;
                info[offset++] = new SuperInfo(ii, jj);
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            LumpIndex jj = getLumpIndex(j);
            assert jj != null;
            for (int i = j + 1; i < stateCount; ++i) {
                LumpIndex ii = getLumpIndex(i);
                assert ii != null;
                info[offset++] = new SuperInfo(ii, jj);
            }
        }

        return info;
    }

    // XT write the function to build a map of index
    private static int[] buildMap(int dim){
        int[] map = new int[dim * dim];

        //Assign Unique Indices to Upper-Triangle Elements (i < j)
        int offset = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                map[i * dim + j] = offset++;
            }
        }

        // Assign Unique Indices to Lower-Triangle Elements (i > j)
        for (int j = 0; j < dim; ++j) {
            for (int i = j + 1; i < dim; ++i) {
                map[i * dim + j] = offset++;
            }
        }
        //map = [ ?,  0,  1,  2,
        //        6,  ?,  3,  4,
        //        7,  9,  ?,  5,
        //        8, 10, 11, ? ];

        // Mark Diagonal Entries as -1
        for (int i = 0; i < dim; ++i) {
            map[i * dim + i] = -1;
        }
        return map;
    }

    private LumpIndex getLumpIndex(int state) {
        for (int i = 0; i < lumps.size(); ++i) {
            Lump lump = lumps.get(i);

            for (int j = 0; j < lump.count; ++j) {
                if (state == lump.states[j]) {
                    for (int k = 0; k < lump.count; ++k){
                        if (state == lump.originalStates[k]) {
                            return new LumpIndex(i, j, k, lump.count);
                        }
                    }

                }
            }
        }

        return null;
    }

    @Override
    public double[] getXBeta() {
        return null;
    }

    @Override
    public LogColumn[] getColumns() {
        return null;
    }

    @Override
    public double[] getRates() {
        if (!CACHE || !ratesKnown) {
            computeRates(rates);
            ratesKnown = true;
        }
        return rates;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("No sub-models");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        ratesKnown = false;
        fireModelChanged(variable, index);
    }

    @Override
    protected void storeState() {
        if (CACHE) {
            if (storedRates == null) {
                storedRates = new double[rates.length];
            }

            System.arraycopy(rates, 0, storedRates, 0, rates.length);
        }
    }

    @Override
    protected void restoreState() {
        if (CACHE) {
            double[] tmp = rates;
            rates = storedRates;
            storedRates = tmp;

            ratesKnown = true;
        }
    }

    @Override
    protected void acceptState() {

    }

    @Override
    public String getReport() {
        double[] rates = getRates();
        return new WrappedVector.Raw(rates).toString();
    }

    @Override
    public Parameter getLogRateParameter() {
        return null;
    }


    public static class Proportion {
        final int source;
        final StateSet destination;
        final Parameter parameter;

        public Proportion(int source, StateSet destination, Parameter parameter) {
            this.source = source;
            this.destination = destination;
            this.parameter = parameter;
        }
    }

    public static class Lump {
        final String id;
        final int[] states;
        final Parameter rates;
        final List<Proportion> proportions;
        final int count;
        final int[] map;
        final int[] originalStates;

//        public Lump(String id, List<Integer> states, Parameter rates, List<Proportion> proportions) {
//            this(id, states.stream().mapToInt(i->i).toArray(), rates, proportions);
//        }

//        public Lump(String id, int[] states, Parameter rates, List<Proportion> proportions) {
//            this(id, states, (Parameter) null, (List<Proportion>) null);
//        }
//        public Lump(String id, int[] states, Parameter rates, List<Proportion> proportions) {
//            this(id, states, (Parameter) null, proportions);
//        }

//        public Lump(String id, List<Integer> states,
//                    Parameter rates, List<Proportion> proportions) {
//            this(id, states.stream().mapToInt(i->i).toArray(), rates, proportions);
//        }

        public Lump(String id, int[] states,
                    Parameter rates, List<Proportion> proportions) {
            this.id = id;
            this.states = states;
            this.rates = rates;
            this.proportions = proportions;
            this.count = states.length;

            // Save the original order of states
            this.originalStates = states.clone(); // XT

            // The states array is sorted in ascending order to ensure a consistent order for indexing.
            //XT COMMENT IT OUT
            Arrays.sort(states);

//            // Allocate map as a Flattened 2D Array
//            this.map = new int[count * count];
            // XT
            this.map = buildMap(count);

        }

    }

    private void matchInfo(SuperInfo info, int i, int j, List<int[]> matches,
                      Parameter parameter, int dimension) {
        if (info.within) {
            if (parameter == info.withinRateParameter && dimension == info.withinRateIndex) {
                matches.add(new int[]{ i, j });
            }
        } else {
            if (parameter == info.acrossProportionParameter && dimension == info.acrossProportionIndex) {
                matches.add(new int[]{ i, j});
            } else if (parameter == acrossRates && dimension == info.acrossRateIndex) {
                matches.add(new int[]{ i, j });
            }
        }
    }

    public List<int[]> searchForParameterAndDimension(Parameter parameter, int dimension) {

        List<int[]> matches = new ArrayList<>();

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                matchInfo(map[k], i, j, matches, parameter, dimension);
                ++k;
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                matchInfo(map[k], i, j, matches, parameter, dimension);
                ++k;
            }
        }

        return matches;
    }

    public class SuperInfo {

        final LumpIndex i;
        final LumpIndex j;
        final boolean within;
        final int withinRateIndex;
        final int acrossRateIndex;
        final int acrossProportionIndex;
        final Parameter withinRateParameter;
        // TX added
        final Parameter acrossProportionParameter;

        SuperInfo(LumpIndex i, LumpIndex j) {
            // in superinfo,just set all null values to -1
            this.i = i;
            this.j = j;
            // Checking if Both Indices Belong to the Same Lump
            this.within = (i.lump == j.lump);
            // retrieves the Lump object associated with i.lump
            Lump lump = lumps.get(i.lump);

            //TX comment out
            // If within is true (i.e., i and j are in the same lump), withinIndex is set using lump.map[]
            // it is the index of the rate in the parameter rate
            this.withinRateIndex = within ?
                    lump.map[i.index * i.count + j.index] : -1;
//            // different lumps) → set acrossIndex = 0;  same lump) → set acrossIndex = -1

            //XT
            int[] lumpMap = buildMap(lumps.size());
            this.acrossRateIndex = !within ?
                    lumpMap[i.lump * lumps.size() + j.lump] : -1;


            int propIndex = (i.lump < j.lump) ?
                    j.lump : j.lump+1;
            this.acrossProportionIndex = !within ?
                    j.index : -1;


            this.withinRateParameter = within ?
                    lump.rates : null;
            // XT question
//            this.acrossProportionParameter = lump.proportions.get(i.index* lumps.size()+ j.lump-1).parameter.getParameterValue(j.index);
            this.acrossProportionParameter = !within ?
                    lump.proportions.get(i.originalIndex* (lumps.size()-1)+ propIndex-1).parameter : null;

        }

        public double getRate() {
            final double rate;
            if (within) {
                rate = withinRateParameter.getParameterValue(withinRateIndex);
            } else {
                double proportion = acrossProportionParameter.getParameterValue(acrossProportionIndex);
                double across = acrossRates.getParameterValue(acrossRateIndex);
                rate = across * proportion;
            }

            return rate;
        }
    }

    // lumpindex tells us a state is in which lump which index, and the number of states in that lump
    private static class LumpIndex {
        //  Identifies which lump the state belongs to.
        final int lump;
        // The position of the state within its lump. Example: If a lump contains the states {2, 3, 5}, the state 3 would have index = 1 (since indexing typically starts at 0)
        final int index;
        // The total number of states in the lump.
        final int count;
        // XT
        final int originalIndex;

        LumpIndex(int lump, int index, int originalIndex, int count) {
            this.lump = lump;
            this.index = index;
            this.count = count;
            this.originalIndex = originalIndex;
        }

        // XT
//        public int getLump() {
//            return lump;
//        }
//        public int getIndex() {
//            return index;
//        }


    }
}
