/*
 * StratifiedTraitOutputFormat.java
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

package dr.evomodel.substmodel;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for enumerating different robust counting output formats in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 */

public enum StratifiedTraitOutputFormat {
    SUM_OVER_SITES("sumOverAllSites", false),
    SUM_OVER_SITES_WITH_CONDITIONED("sumOverAllSitesWithUnconditioned", true),
    PER_SITE("perSite", false),
    PER_SITE_WITH_UNCONDITIONED("perSiteWithUnconditioned", true),
    ARBITRARY_SITES("arbitrarySites", false),
    ARBITRARY_SITES_WITH_UNCONDITIONED("arbitrarySitesWithUnconditioned", true);

    private StratifiedTraitOutputFormat(String text, boolean supportsUnconditioned) {
        this.text = text;
        this.supportsUnconditioned = supportsUnconditioned;
    }

    public String getText() {
        return text;
    }

    public boolean getSupportsUnconditioned() {
        return supportsUnconditioned;
    }

    public static StratifiedTraitOutputFormat parseFromString(String text) {
        for (StratifiedTraitOutputFormat format : StratifiedTraitOutputFormat.values()) {
            if (format.getText().compareToIgnoreCase(text) == 0)
                return format;
        }
        return null;
    }

    private final String text;
    private final boolean supportsUnconditioned;
}
