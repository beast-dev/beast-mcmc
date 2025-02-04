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

import java.util.Arrays;
import java.util.List;

/**
 * @author Xinghua Tao
 * @author Marc Suchard
 */
public class StronglyLumpableCtmcRates extends AbstractModel implements LogAdditiveCtmcRateProvider {

    private final SuperInfo[] map;
    private final List<Lump> lumps;
    private final Parameter acrossRates;
//    private final DataType dataType;
    private final int stateCount;
    private final int dim;

    private double[] rates;
    private double[] storedRates;

    private boolean ratesKnown = false;
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

        for (Lump lump : lumps) {
            if (lump.rates != null) {
                addVariable(lump.rates);
            }

            if (lump.proportions != null) {
                addVariable(lump.proportions);
            }
        }

        addVariable(acrossRates);

        computeRates(rates); // TODO remove
    }

    private double getRate(int offset) {
        SuperInfo info = map[offset];
        final double rate;
        if (info.within) {
            rate = info.withinRateParameter.getParameterValue(info.withinIndex);
        } else {
            rate = 0;
        }

        return rate;
    }

    private void computeRates(double[] rates) {
        for (int i = 0; i < dim; ++i) {
            rates[i] = getRate(i);
        }
    }

    private SuperInfo[] buildSuperMap() {

        SuperInfo[] info = new SuperInfo[dim];

        int offset = 0;
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

    private LumpIndex getLumpIndex(int state) {
        for (int i = 0; i < lumps.size(); ++i) {
            Lump lump = lumps.get(i);

            for (int j = 0; j < lump.count; ++j) {
                if (state == lump.states[j]) {
                    return new LumpIndex(i, j, lump.count);
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

    public static class StateSet {
        final String id;
        final int[] states;
        final DataType dataType;

        public StateSet(String id,
                        List<Integer> states,
                        DataType dataType) {
            this(id, states.stream().mapToInt(i->i).toArray(), dataType);
        }

        public StateSet(String id,
                        int[] states,
                        DataType dataType) {
            this.id = id;
            this.states = states;
            this.dataType = dataType;
        }

        public int size() { return states.length; }

        public int[] states() { return states; }
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
        final Parameter proportions;
        final int count;
        final int[] map;

        public Lump(String id, List<Integer> states, Parameter rates, List<Proportion> proportions) {
            this(id, states.stream().mapToInt(i->i).toArray(), rates, proportions);
        }

        public Lump(String id, int[] states, Parameter rates, List<Proportion> proportions) {
            this(id, states, (Parameter) null, (Parameter) null);
        }

        public Lump(String id, List<Integer> states,
                    Parameter rates, Parameter proportions) {
            this(id, states.stream().mapToInt(i->i).toArray(), rates, proportions);
        }

        public Lump(String id, int[] states,
                    Parameter rates, Parameter proportions) {
            this.id = id;
            this.states = states;
            this.rates = rates;
            this.proportions = proportions;
            this.count = states.length;

            Arrays.sort(states);

            this.map = new int[count * count];

            int offset = 0;
            for (int i = 0; i < count; ++i) {
                for (int j = i + 1; j < count; ++j) {
                    map[i * count + j] = offset++;
                }
            }

            for (int j = 0; j < count; ++j) {
                for (int i = j + 1; i < count; ++i) {
                    map[i * count + j] = offset++;
                }
            }

            for (int i = 0; i < count; ++i) {
                map[i * count + i] = -1;
            }
        }
    }

    private class SuperInfo {

        final LumpIndex i;
        final LumpIndex j;
        final boolean within;
        final int withinIndex;
        final int acrossIndex;
        final int acrossProportion;

        final Parameter withinRateParameter;
        final Parameter acrossProportionParameter;

        SuperInfo(LumpIndex i, LumpIndex j) {
            this.i = i;
            this.j = j;
            this.within = (i.lump == j.lump);

            Lump lump = lumps.get(i.lump);

            this.withinIndex = within ?
                    lump.map[i.index * i.count + j.index] : -1;
            this.acrossIndex = !within ?
                    0 : -1;
            this.acrossProportion = !within ?
                    0 : -1;

            this.withinRateParameter = lump.rates;
            this.acrossProportionParameter = lump.proportions;
        }
    }

    private static class LumpIndex {
        final int lump;
        final int index;
        final int count;

        LumpIndex(int lump, int index, int count) {
            this.lump = lump;
            this.index = index;
            this.count = count;
        }
    }
}
