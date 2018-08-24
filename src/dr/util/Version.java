/*
 * Version.java
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

/**
 * Version last changed 2004/05/07 by AER
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: Version.java,v 1.13 2005/07/11 14:06:25 rambaut Exp $
 */
public interface Version {

    String getVersion();

	String getVersionString();

	String getBuildString();

	String getDateString();

    String[] getCredits();

    String getHTMLCredits();

    class Utils {
        /**
         * Is version1 more recent (higher) than version2?
         * @param version1
         * @param version2
         * @return
         */
        public static boolean isMoreRecent(String version1, String version2) {
            String[] v1 = version1.split("\\.");
            String[] v2 = version2.split("\\.");

            for (int i = 0; i < Math.min(v1.length, v2.length); i++) {
                if (Integer.parseInt(v1[i]) < Integer.parseInt(v2[i])) {
                    return false;
                } else if (Integer.parseInt(v1[i]) > Integer.parseInt(v2[i])) {
                    return true;
                }
            }

            return false;
        }
    }
}
