/*
 * Version.java
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

package dr.util;

/**
 * Version last changed 2004/05/07 by AER
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public interface Version {

    String getVersion();

    String getVersionString();

    String getBuildString();

    String getDateString();

    String[] getCredits();

    String getHTMLCredits();

    class Utils {
        private static String[] splitVersion(String version) {
            if (version.contains("-")) {
                // in case it has a -beta or similar suffix
                version = version.substring(0, version.indexOf("-"));
            }
            return version.trim().split("\\.");
        }

        /**
         * Is version1 more recent (higher) than version2?
         * @param version1
         * @param version2
         * @return
         */
        public static boolean isMoreRecent(String version1, String version2) {
            String[] v1 = splitVersion(version1);
            String[] v2 = splitVersion(version2);

            try {
                for (int i = 0; i < Math.min(v1.length, v2.length); i++) {
                    if (Integer.parseInt(v1[i]) < Integer.parseInt(v2[i])) {
                        return false;
                    } else if (Integer.parseInt(v1[i]) > Integer.parseInt(v2[i])) {
                        return true;
                    }
                }
            } catch (NumberFormatException nfe) {
                // if the numbers can't be parsed just allow the parsing to proceed.
                return false;
            }

            return false;
        }
    }
}
