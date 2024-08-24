/*
 * BeagleVersion.java
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

package dr.app.beast;

import beagle.BeagleInfo;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Arrays;
import java.util.Collections;
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
        return Collections.singletonList(CommonCitations.AYRES_2019_BEAGLE);
    }

}
