/*
 * BeastVersion.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.util.Version;

/**
 * This class provides a mechanism for returning the version number of the
 * dr software. It relies on the administrator of the dr source using the
 * module tagging system in CVS. The method getVersionString() will return
 * the version of dr under the following condition: <BR>
 * 1. the dr source has been checked out *by tag* before being packaged for
 * distribution.
 * <p/>
 * Version last changed 2009/08/1 by AER
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class BeastVersion implements Version {

    /**
     * Version string: assumed to be in format x.x.x
     */
    private static final String VERSION = "1.5.1";
    private static final String DATE_STRING = "2002-2009";

    /**
     * this used to parse the CVS ID string but there is no equivalent
     * for SVN. We must increment this manually.
     */
    private static final String BUILD_ID = "Build r2015";

    public String getVersionString() {
        return "v" + VERSION;
    }

    public String getDateString() {
        return DATE_STRING;
    }

    public String[] getCredits() {
        return new String[] {
                "Designed and developed by",
                "Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard",
                "",
                "Department of Computer Science",
                "University of Auckland",
                "alexei@cs.auckland.ac.nz",
                "",
                "Institute of Evolutionary Biology",
                "University of Edinburgh",
                "a.rambaut@ed.ac.uk",
                "",
                "David Geffen School of Medicine",
                "University of California, Los Angeles",
                "msuchard@ucla.edu",
                "",
                "Downloads, Help & Resources:",

                "\thttp://beast.bio.ed.ac.uk/",
                "",
                "Source code distributed under the GNU Lesser General Public License:",
                "\thttp://code.google.com/p/beast-mcmc/",
                "",
                "BEAST developers:",
                "\tAlex Alekseyenko, Erik Bloomquist, Joseph Heled, Sebastian Hoehna, ",
                "\tPhilippe Lemey, Gerton Lunter, Sidney Markowitz, Vladimir Minin, ",
                "\tMichael Defoin Platel, Oliver Pybus, Walter Xie",
                "",
                "Thanks to:",
                "\tRoald Forsberg, Beth Shapiro and Korbinian Strimmer"};
    }

    public String getHTMLCredits() {
        return
                "<p>Designed and developed by<br>" +
                        "Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard</p>" +
                        "<p>Department of Computer Science, University of Auckland<br>" +
                        "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                        "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                        "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                        "<p>David Geffen School of Medicine, University of California, Los Angeles<br>" +
                        "<a href=\"mailto:msuchard@ucla.edu\">msuchard@ucla.edu</a></p>" +
                        "<p><a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                        "<p>Source code distributed under the GNU LGPL:<br>" +
                        "<a href=\"http://code.google.com/p/beast-mcmc/\">http://code.google.com/p/beast-mcmc/</a></p>" +
                        "<p>BEAST developers:<br>" +
                        "Alex Alekseyenko, Erik Bloomquist, Joseph Heled, Sebastian Hoehna, Philippe Lemey,<br>" +
                        "Gerton Lunter, Sidney Markowitz, Vladimir Minin, Michael Defoin Platel, Oliver Pybus,<br>" +
                        "Walter Xie</p>" +
                        "<p>Thanks to Roald Forsberg, Beth Shapiro and Korbinian Strimmer</p>";
    }

    public String getBuildString() {
        return BUILD_ID;
    }
}
