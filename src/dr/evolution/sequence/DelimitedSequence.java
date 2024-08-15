/*
 * DelimitedSequence.java
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

package dr.evolution.sequence;

import dr.evolution.datatype.DataType;
import dr.evolution.util.Taxon;

/**
 * @author Marc A. Suchard
 */
public class DelimitedSequence extends Sequence {

    public DelimitedSequence(Taxon taxon, String sequence, DataType dataType) {
        super(taxon, sequence);
        this.dataType = dataType;
        codes = sequence.split(dataType.getDelimiter());
        original = sequence;
        sequenceString = null;
    }

    @Override
    public char getChar(int index) {
        throw new RuntimeException("Not single character available");
    }

    @Override
    public String getSequenceString() {
        return original;
    }

    @Override
    public int getLength() {
        return codes.length;
    }

    @Override
    public int getState(int index) {
        return dataType.getState(codes[index]);
    }

    @Override
    public void setState(int index, int state) {
        throw new RuntimeException("Not implemented");
    }

    private final String[] codes;

    private final String original;
}
