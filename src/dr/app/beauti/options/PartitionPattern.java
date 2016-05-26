/*
 * PartitionPattern.java
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

package dr.app.beauti.options;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class PartitionPattern extends AbstractPartitionData {
    private static final long serialVersionUID = 6631884346312113086L; // microsatellite

    private final Patterns patterns;

    public PartitionPattern(BeautiOptions options, String name, String fileName, Patterns patterns) {
        super(options, name, fileName);
        this.patterns = patterns;

        this.traits = null;

        calculateMeanDistance(patterns);
    }

    public Patterns getPatterns() {
        return patterns;
    }

    public TaxonList getTaxonList() {
        return getPatterns();
    }

    public int getSiteCount() {
        return 1;
    }

    public DataType getDataType() {
        if (patterns != null) {
            return patterns.getDataType();
        } else {
            throw new RuntimeException("patterns should not be null");
        }
    }

    public String getDataDescription() {
        if (patterns != null) {
            return patterns.getDataType().getDescription();
        } else {
            throw new RuntimeException("patterns should not be null");
        }
    }

    public String getPrefix() {
        String prefix = "";
        if (options.getPartitionPattern().size() > 1) { // getPartitionPattern() already excludes traits and PartitionData
            prefix += getName() + ".";
        }
        return prefix;
    }

    public void setName(String name) {
        this.name = name;
        options.microsatelliteOptions.initParametersAndOperators();
    }
}
