/*
 * CharSetAlignment.java
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

import dr.app.beauti.util.NexusApplicationImporter;
import dr.evolution.sequence.Sequence;

/**
 * Created by IntelliJ IDEA.
 * User: alexei
 * Date: 1/08/11
 * Time: 8:04 PM
 * To change this template use File | Settings | File Templates.
 */
public class CharSetAlignment extends SimpleAlignment {

    public CharSetAlignment(NexusApplicationImporter.CharSet charset, Alignment parentAlignment) {
        setId(charset.getName());

        for (int i = 0; i < parentAlignment.getSequenceCount(); i++) {
            Sequence sequence = parentAlignment.getSequence(i);
            String sequenceString = parentAlignment.getAlignedSequenceString(i);

            String filteredSequence = filter(charset, sequenceString);

            addSequence(new Sequence(sequence.getTaxon(), filteredSequence));
        }

        setDataType(parentAlignment.getDataType());
    }

    private String filter(NexusApplicationImporter.CharSet charset, String sequenceString) {
        StringBuilder filtered = new StringBuilder();
        for (NexusApplicationImporter.CharSetBlock block : charset.getBlocks()) {
            int toSite = block.getToSite();
            if (toSite <= 0) {
                toSite = sequenceString.length();
            }
            for (int i = block.getFromSite(); i <= toSite; i += block.getEvery()) {
                // the -1 comes from the fact that charsets are indexed from 1 whereas strings are indexed from 0
                filtered.append(sequenceString.charAt(i - 1));
            }
        }
        return filtered.toString();
    }

}
