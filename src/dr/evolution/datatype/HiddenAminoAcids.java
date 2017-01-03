/*
 * HiddenAminoAcids.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

/**
 * @author Marc A. Suchard
 */

public class HiddenAminoAcids extends AminoAcids implements HiddenDataType {

    public static final String DESCRIPTION = "hiddenAminoAcid";

    public static final HiddenAminoAcids AMINO_ACIDS_HIDDEN_1 = new HiddenAminoAcids(1);
    public static final HiddenAminoAcids AMINO_ACIDS_HIDDEN_2 = new HiddenAminoAcids(2);
    public static final HiddenAminoAcids AMINO_ACIDS_HIDDEN_3 = new HiddenAminoAcids(3);
    public static final HiddenAminoAcids AMINO_ACIDS_HIDDEN_4 = new HiddenAminoAcids(4);

    /**
     * Private constructor - DEFAULT_INSTANCE provides the only instance
     */
    private HiddenAminoAcids(int hiddenClassCount) {
        super();
        this.hiddenClassCount = hiddenClassCount;
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {
        return Utils.getStateSet(state, stateCount, hiddenClassCount, AminoAcids.INSTANCE);
    }

    public int getStateCount() {
        return stateCount * hiddenClassCount;
    }

    public int getHiddenClassCount() {
        return hiddenClassCount;
    }

    private int hiddenClassCount;
}