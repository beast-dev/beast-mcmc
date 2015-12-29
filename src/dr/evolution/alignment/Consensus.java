/*
 * Consensus.java
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
import dr.evolution.util.Taxon;

/**
 * @author Alexei Drummond
 *
 * @version $Id: Consensus.java,v 1.6 2005/04/20 21:26:18 rambaut Exp $
 */
public class Consensus {

    int[] counts;
    int[] consensus;
    int total;
    String name;
    DataType dataType = null;

    public Consensus(String name, Alignment alignment, boolean ignoreGaps) {

        this.name = name;
        dataType = alignment.getDataType();

        int[][] frequencies = new int[alignment.getSiteCount()][dataType.getAmbiguousStateCount()];
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            for (int j = 0; j < alignment.getSiteCount(); j++) {
                int state = alignment.getState(i, j);
                if (ignoreGaps) {
                    if (state < dataType.getStateCount()) {
                        frequencies[j][state] += 1;
                    }
                } else {
                    frequencies[j][state] += 1;
                }
            }
        }

        counts = new int[alignment.getSiteCount()];
        total = alignment.getSequenceCount();
        consensus = new int[alignment.getSiteCount()];
        for (int i = 0; i < alignment.getSiteCount(); i++) {

            int maxState = 0;
            int maxFreq = frequencies[i][0];
            for (int j = 1; j < frequencies[i].length; j++) {
                int freq = frequencies[i][j];
                if (freq > maxFreq) {
                    maxState = j;
                    maxFreq = freq;
                }
            }
            consensus[i] = maxState;
            counts[i] = maxFreq;
        }
    }

    /**
     * @param site
     * @return the probability that a sequences in the original alignment has
     * the same value as the consensus sequences at the given site.
     */
    public double getReliability(int site) {
        return (double)counts[site]/(double)total;
    }

    public int getState(int site) {
        return consensus[site];    
    }

    public final Sequence getConsensusSequence() {
        StringBuffer buffer = new StringBuffer();
        for (int i = 0; i < consensus.length; i++) {
            buffer.append(dataType.getChar(getState(i)));
        }
        Sequence sequence = new Sequence(new Taxon(name),buffer.toString());
        sequence.setDataType(dataType);
        return sequence;
    }
}
