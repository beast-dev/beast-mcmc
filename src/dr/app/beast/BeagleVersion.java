/*
 * BeastVersion.java
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

package dr.app.beast;

import beagle.BeagleInfo;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Arrays;
import java.util.List;

/**
 * This class provides a mechanism for returning the version number and
 * citation of the BEAGLE high-performance computational library.
 *
 * @author Guy Baele
 *
 * $Id$
 */

public class BeagleVersion implements Citable {

    public static final BeagleVersion INSTANCE = new BeagleVersion();

    public String getVersion() {
        return BeagleInfo.getVersion();
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.COMPUTATIONAL_LIBRARY;
    }

    @Override
    public String getDescription() {
        return "BEAGLE citation";
    }

    @Override
    public List<Citation> getCitations() {
        return Arrays.asList(CITATIONS);
    }

    public static Citation[] CITATIONS = new Citation[] {
            new Citation(
                    new Author[]{
                            new Author("DL", "Ayres"),
                            new Author("MP", "Cummings"),
                            new Author("G", "Baele"),
                            new Author("AE", "Darling"),
                            new Author("PO", "Lewis"),
                            new Author("DL", "Swofford"),
                            new Author("JP", "Huelsenbeck"),
                            new Author("P", "Lemey"),
                            new Author("A", "Rambaut"),
                            new Author("MA", "Suchard")
                    },
                    "BEAGLE 3: Improved performance, scaling, and usability for a high-performance computing library for statistical phylogenetics",
                    2019,
                    "Systematic Biology",
                    68, 1052, 1061,
                    "10.1093/sysbio/syz020")
    };

}
