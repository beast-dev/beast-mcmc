/*
 * ImportReport.java
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

package dr.app.beauti.util;

import java.util.ArrayList;
import java.util.List;

/**
 * What an import of a BEAST XML file was able to make of it. The import never fails over a part it does not
 * understand; it notes it here instead, so that the user can be told which settings they still have to
 * check.
 *
 * @author Andrew Rambaut
 */
public class ImportReport {

    private final List<String> recovered = new ArrayList<String>();
    private final List<String> notRecovered = new ArrayList<String>();

    /**
     * Notes a setting taken from the file, as it should appear in a list ("the substitution model, HKY").
     */
    public void recovered(String description) {
        recovered.add(description);
    }

    /**
     * Notes something the file carried that could not be turned into a setting, with the reason.
     */
    public void notRecovered(String description) {
        notRecovered.add(description);
    }

    public List<String> getRecovered() {
        return recovered;
    }

    public List<String> getNotRecovered() {
        return notRecovered;
    }

    /**
     * @return true if nothing was recovered and nothing was missed, i.e. no settings were looked for
     */
    public boolean isEmpty() {
        return recovered.isEmpty() && notRecovered.isEmpty();
    }

    /**
     * @return true if everything the file carried was understood
     */
    public boolean isComplete() {
        return notRecovered.isEmpty();
    }

    /**
     * @return the report as HTML for a dialog, or null if there is nothing to say
     */
    public String toHTML() {
        if (isEmpty()) {
            return null;
        }
        StringBuilder html = new StringBuilder("<html>");
        if (!recovered.isEmpty()) {
            html.append("<p>These settings were read from the file:</p><ul>");
            for (String line : recovered) {
                html.append("<li>").append(escape(line)).append("</li>");
            }
            html.append("</ul>");
        }
        if (!notRecovered.isEmpty()) {
            html.append("<p>These were not, and are left as BEAUti's defaults:</p><ul>");
            for (String line : notRecovered) {
                html.append("<li>").append(escape(line)).append("</li>");
            }
            html.append("</ul><p>Check them before generating a file.</p>");
        }
        return html.append("</html>").toString();
    }

    private static String escape(String text) {
        return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;");
    }

    @Override
    public String toString() {
        return recovered.size() + " settings read, " + notRecovered.size() + " not";
    }
}
