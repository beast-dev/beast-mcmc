/*
 * PairedDataType.java
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

package dr.evolution.datatype;

import dr.evomodel.treedatalikelihood.discrete.NodeHeightToRatiosTransformDelegate;

import java.util.HashMap;
import java.util.Map;

/**
 * PoMo data type
 *
 * @author Xiang Ji
 * @author Nicola De Maio
 * @author Ben Redelings
 * @author Marc A. Suchard
 */
public class PolymorphismAwareDataType extends DataType {

    private final DataType baseDataType;

    private final int virtualPopSize;

    public static final String DESCRIPTION = "polymorphismAware";

    private final Map<String, Integer> sequenceStateMap;

    public PolymorphismAwareDataType(DataType baseDataType, int virtualPopSize) {
        this.baseDataType = baseDataType;
        this.virtualPopSize = virtualPopSize;
        this.sequenceStateMap = new HashMap<>();
        buildSequenceStateMap();
    }

    private void buildSequenceStateMap() {

        for (int i = 0; i < baseDataType.getStateCount(); i++) {
            sequenceStateMap.put(baseDataType.getCode(i), i);
        }
        int state = baseDataType.getStateCount();
        for (int i = 0; i < baseDataType.getStateCount() - 1; i++) {
            String seq1 = baseDataType.getCode(i);
            for (int j = i + 1; j < baseDataType.getStateCount(); j++) {
                String seq2 = baseDataType.getCode(j);
                for (int k = 1; k < virtualPopSize; k++) {
                    String seq = seq1 + ":" + Integer.toString(virtualPopSize - k)
                            + "|" + seq2 + ":" + Integer.toString(k);
                    sequenceStateMap.put(seq, state);
                    state++;
                }
            }
        }
    }

    @Override
    public char[] getValidChars() {
        return new char[0];
    }

    public final int getState(int[] states, int[] counts) {
        if (states.length == 1) {
            return states[0];
        } else {
            if (states.length != 2 || counts.length != 2
                    || counts[0] + counts[1] != virtualPopSize || counts[0] < 1 || counts[1] < 1) {
                throw new RuntimeException("Illegal states or counts.");
            }
            final int firstState = states[0] < states[1] ? states[0] : states[1];
            final int secondState = states[0] < states[1] ? states[1] : states[0];
//            final int firstCount = states[0] < states[1] ? counts[0] : counts[1];
            final int secondCount = states[0] < states[1] ? counts[1] : counts[0];

            return baseDataType.getStateCount() + (virtualPopSize - 1) * (secondState - firstState - 1 + (2 * baseDataType.getStateCount() - firstState - 1) * firstState / 2)
                    + secondCount - 1;
        }
    }

    public DataType getBaseDataType() {
        return baseDataType;
    }

    @Override
    public int getStateCount() {
        final int baseStateCount = baseDataType.getStateCount();
        return baseStateCount + (virtualPopSize - 1) * baseStateCount * (baseStateCount - 1) / 2;
    }

    @Override
    public String getDescription() {
        return DESCRIPTION;
    }

    @Override
    public int getType() {
        return 0;
    }
}
