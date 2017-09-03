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

import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.Version;

import java.util.Collections;
import java.util.List;

/**
 * This class provides a mechanism for returning the version number of the
 * dr software.
 *
 * This is manually updated as required. The REVISION string is no longer used
 * since switching to GitHub.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * $Id$
 */
public class BeastVersion implements Version, Citable {

    public static final BeastVersion INSTANCE = new BeastVersion();

    /**
     * Version string: assumed to be in format x.x.x
     */
    private static final String VERSION = "1.10.0";

    private static final String DATE_STRING = "2002-2017";

    private static final boolean IS_PRERELEASE = true;

    // this is now being manually updated since the move to GitHub. 7 digits of GitHub hash.
    private static final String REVISION = "VEME2017";

    public String getVersion() {
        return VERSION;
    }

    public String getVersionString() {
        return "v" + VERSION + (IS_PRERELEASE ? " Prerelease #" + REVISION : "");
    }

    public String getDateString() {
        return DATE_STRING;
    }

    public String[] getCredits() {
        return new String[]{
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

                "\thttp://beast.community",
                "",
                "Source code distributed under the GNU Lesser General Public License:",
                "\thttp://github.com/beast-dev/beast-mcmc",
                "",
                "BEAST developers:",
                "\tAlex Alekseyenko, Guy Baele, Trevor Bedford, Filip Bielejec, Erik Bloomquist, Matthew Hall,",
                "\tJoseph Heled, Sebastian Hoehna, Denise Kuehnert, Philippe Lemey, Wai Lok Sibon Li,",
                "\tGerton Lunter, Sidney Markowitz, Vladimir Minin, Michael Defoin Platel,",
                "\tOliver Pybus, Chieh-Hsi Wu, Walter Xie",
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
                        "<p><a href=\"http://beast.community\">http://beast.community</a></p>" +
                        "<p>Source code distributed under the GNU LGPL:<br>" +
                        "<a href=\"http://github.com/beast-dev/beast-mcmc\">http://github.com/beast-dev/beast-mcmc</a></p>" +
                        "<p>BEAST developers:<br>" +
                        "Alex Alekseyenko, Guy Baele, Trevor Bedford, Filip Bielejec, Erik Bloomquist, Matthew Hall,<br>"+
                        "Joseph Heled, Sebastian Hoehna, Denise Kuehnert, Philippe Lemey, Wai Lok Sibon Li,<br>"+
                        "Gerton Lunter, Sidney Markowitz, Vladimir Minin, Michael Defoin Platel,<br>"+
                        "Oliver Pybus, Chieh-Hsi Wu, Walter Xie</p>" +
                        "<p>Thanks to Roald Forsberg, Beth Shapiro and Korbinian Strimmer</p>";
    }

    public String getBuildString() {
        return "https://github.com/beast-dev/beast-mcmc/commit/" + REVISION;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "BEAST primary citation";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("AJ", "Drummond"),
                    new Author("MA", "Suchard"),
                    new Author("Dong", "Xie"),
                    new Author("A", "Rambaut")
            },
            "Bayesian phylogenetics with BEAUti and the BEAST 1.7",
            2012,
            "Mol Biol Evol",
            29, 1969, 1973,
            "10.1093/molbev/mss075");

}
