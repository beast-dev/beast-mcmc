/*
 * ExtractPairs.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;

/**
 * @author Alexei Drummond
 *
 * @version $Id: ExtractPairs.java,v 1.3 2005/04/20 21:26:18 rambaut Exp $
 */
public class ExtractPairs {

    Alignment alignment;

    public ExtractPairs(Alignment alignment) {

        this.alignment = alignment;
    }

    public Alignment getPairAlignment(int x, int y) {

        SimpleAlignment pairAlignment = new SimpleAlignment();

        StringBuffer sequence0 = new StringBuffer();
        StringBuffer sequence1 = new StringBuffer();

        DataType dataType = alignment.getDataType();
        int stateCount = dataType.getStateCount();

        for (int i = 0; i < alignment.getSiteCount(); i++) {

            int s0 = alignment.getState(x,i);
            int s1 = alignment.getState(y,i);

            char c0 = dataType.getChar(s0);
            char c1 = dataType.getChar(s1);

            if (s0 < stateCount || s1 < stateCount) {
                sequence0.append(c0);
                sequence1.append(c1);
            }
        }

        // trim hanging ends on left
        int left = 0;
        while (
            (dataType.getState(sequence0.charAt(left)) >= stateCount) ||
            (dataType.getState(sequence1.charAt(left)) >= stateCount)) {
            left += 1;
        }

        // trim hanging ends on right
        int right = sequence0.length()-1;
        while (
            (dataType.getState(sequence0.charAt(right)) >= stateCount) ||
            (dataType.getState(sequence1.charAt(right)) >= stateCount)) {
            right -= 1;
        }

        if (right < left) return null;

        String sequenceString0 = sequence0.substring(left,right+1);
        String sequenceString1 = sequence1.substring(left,right+1);

        pairAlignment.addSequence(new Sequence(alignment.getTaxon(x),sequenceString0));
        pairAlignment.addSequence(new Sequence(alignment.getTaxon(y),sequenceString1));

        return pairAlignment;
    }

}
