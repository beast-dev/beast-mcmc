/*
 * PolymorphismAwareDataType.java
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

import dr.evolution.sequence.UncertainSequence;
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

    public PolymorphismAwareDataType(DataType baseDataType, int virtualPopSize) {
        this.baseDataType = baseDataType;
        this.virtualPopSize = virtualPopSize;
        this.stateCount = getStateCount();
    }

    public static String getDataTypeDescription(DataType baseDataType, int virtualPopSize) {
        return DESCRIPTION + "." + baseDataType.getDescription() + ".K" + virtualPopSize;
    }

    public int getVirtualPopSize() {
        return virtualPopSize;
    }

    private int getBaseState(String stateString) {
        if (stateString.length() == 1) {
            return baseDataType.getState(stateString.charAt(0));
        } else if (baseDataType instanceof Codons) {
            return ((Codons) baseDataType).getState(stateString.charAt(0), stateString.charAt(1), stateString.charAt(2));
        } else {
            throw new RuntimeException("Unsupported data type: " + baseDataType);
        }
    }

    public int getState(UncertainSequence.UncertainCharacterList characters) {
        assert(characters.size() < 3);
        if (characters.size() == 1) {
            return getBaseState(characters.get(0).getSequenceString());
        } else if (characters.size() == 2) {
            final int firstState = getBaseState(characters.get(0).getSequenceString());
            final int secondState = getBaseState(characters.get(1).getSequenceString());
            final int firstCount = (int) characters.get(0).getWeight();
            final int secondCount = (int) characters.get(1).getWeight();
            return getState(new int[]{firstState, secondState}, new int[]{firstCount, secondCount});
        } else {
            throw new RuntimeException("Unexpected more than 3 characters in sequence");
        }
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    public final int getState(int[] states, int[] counts) {
        if (states.length == 1) {
            return states[0];
        } else {
            if (states.length != 2 || counts.length != 2
                    || counts[0] + counts[1] != virtualPopSize || counts[0] < 0 || counts[1] < 0) {
                throw new RuntimeException("Illegal states or counts.");
            }
            if (counts[0] == 0) {
                return states[1];
            }
            if (counts[1] == 0) {
                return states[0];
            }
            final int firstState = states[0] < states[1] ? states[0] : states[1];
            final int secondState = states[0] < states[1] ? states[1] : states[0];
            final int firstCount = states[0] < states[1] ? counts[0] : counts[1];
//            final int secondCount = states[0] < states[1] ? counts[1] : counts[0];

            final int state = baseDataType.getStateCount() + ((2 * baseDataType.getStateCount() - firstState - 1) * firstState / 2 + secondState - firstState - 1) * (virtualPopSize - 1) + firstCount - 1;

            return state;
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
