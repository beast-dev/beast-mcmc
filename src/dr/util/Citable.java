/*
 * Citable.java
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

package dr.util;

import java.util.Collection;
import java.util.List;
import java.util.Map;

/**
 * Interface for associating a list of citations with an object
 *
 * @author Marc A. Suchard
 */

public interface Citable {

    Citation.Category getCategory();

    String getDescription();

    /**
     * @return a list of citations associated with this object
     */
    List<Citation> getCitations();

    class Utils {
        public static String getCitationString(Citable citable, String prepend, String postpend) {
            if (citable.getCitations().size() == 0) {
                return null;
            }
            StringBuilder builder = new StringBuilder();
            for (Citation citation : citable.getCitations()) {
                builder.append(prepend);
                builder.append(citation.toString());
                builder.append(postpend);
            }
            return builder.toString();
        }

        public static String getCitationString(Citable citable) {
            return getCitationString(citable, DEFAULT_PREPEND, DEFAULT_POSTPEND);
        }

        public static final String DEFAULT_PREPEND = "\t\t";
        public static final String DEFAULT_POSTPEND = "\n";
    }
}
