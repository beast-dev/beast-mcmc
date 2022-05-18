/*
 * StateDependentNucleotides.java
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

package dr.evolution.datatype;

public class StateDependentNucleotides extends Nucleotides implements StateDependentDataType {

    public static final String DESCRIPTION = "stateDependentNucleotide";

    static final StateDependentNucleotides NUCLEOTIDE_STATE_DEPENDENT_1 = new StateDependentNucleotides(1);
    static final StateDependentNucleotides NUCLEOTIDE_STATE_DEPENDENT_2 = new StateDependentNucleotides(2);
    static final StateDependentNucleotides NUCLEOTIDE_STATE_DEPENDENT_3 = new StateDependentNucleotides(3);
    static final StateDependentNucleotides NUCLEOTIDE_STATE_DEPENDENT_4 = new StateDependentNucleotides(4);
    static final StateDependentNucleotides NUCLEOTIDE_STATE_DEPENDENT_8 = new StateDependentNucleotides(8);

    /**
     * Private constructor - DEFAULT_INSTANCE provides the only instance
     */
    private StateDependentNucleotides(int stateDependentClassCount) {
        super();
        this.stateDependentClassCount = stateDependentClassCount;
    }

    public int getState(int state, int dependentClass) {
        return StateDependentDataType.Utils.getState(state, dependentClass);
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state, int dependentClass) {
        return StateDependentDataType.Utils.getStateSet(state, stateCount, stateDependentClassCount, Nucleotides.INSTANCE);
    }

    public String getCode(int state, int dependentClass) {
        return StateDependentDataType.getCodeImpl(state, stateCount, super::getCode);
    }

    public String getCodeWithoutDependentState(int state) {
        return StateDependentDataType.getCodeWithoutDependentStateImpl(state, stateCount, super::getCode);
    }

    public int getStateCount() {
        return stateCount * stateDependentClassCount;
    }

    private int stateDependentClassCount;

    public int getDependentClassCount() {
        return stateDependentClassCount;
    }

    public String getDescription() {
        return DESCRIPTION + stateDependentClassCount;
    }
}