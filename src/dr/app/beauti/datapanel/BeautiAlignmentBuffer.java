/*
 * BeautiAlignmentBuffer.java
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

package dr.app.beauti.datapanel;


import java.util.List;

import dr.app.beauti.alignmentviewer.AlignmentBuffer;
import dr.app.beauti.alignmentviewer.AlignmentBufferListener;
import dr.evolution.alignment.Alignment;
import dr.evolution.datatype.DataType;
import dr.evolution.sequence.Sequence;

/**
 * @author Andrew Rambaut
 * @version $Id: SimpleAlignmentBuffer.java,v 1.2 2005/12/11 22:41:25 rambaut Exp $
 */
public class BeautiAlignmentBuffer implements AlignmentBuffer {

    public BeautiAlignmentBuffer(Alignment alignment) {
        this.alignment = alignment;
        DataType type = alignment.getDataType();

        stateTable = new String[type.getAmbiguousStateCount()];
        for (int i = 0; i < stateTable.length; i++) {
            stateTable[i] = Character.toString(type.getChar(i));
        }

        gapState = (byte)type.getGapState();
    }

    public int getSequenceCount() {
        return alignment.getSequenceCount();
    }

    public int getSiteCount() {
        return alignment.getSiteCount();
    }

    public String getTaxonLabel(int i) {
        return alignment.getTaxonId(i);
    }

    public String[] getStateTable() {
        return stateTable;
    }

    public void getStates(int sequenceIndex, int fromSite, int toSite, byte[] states) {
        Sequence sequence = alignment.getSequence(sequenceIndex);
        int j = 0;
        for (int i = fromSite; i <= toSite; i++) {
            states[j] = (byte)sequence.getState(i);
            j++;
        }
    }

    public void addAlignmentBufferListener(AlignmentBufferListener listener) {
    }

    private final Alignment alignment;
    private final String[] stateTable;
    private final byte gapState;

}
