/*
 * BeastVersion.java
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

import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.Version;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

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
    private static final String VERSION = "10.5.0";

    private static final String DATE_STRING = "2002-2025";

    private static final boolean IS_PRERELEASE = false;

    public String getVersion() {
        return VERSION;
    }

    public String getVersionString() {
        return "v" + VERSION + (IS_PRERELEASE ? " Prerelease #" + getRevision() : "");
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
                "Institute of Ecology and Evolution",
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
                "\tAlex Alekseyenko, Daniel Ayres, Guy Baele, Trevor Bedford, Filip Bielejec, Erik Bloomquist, ",
                "\tLuiz Max Carvalho, Gabriela Cybis, Mandev Gill, Matthew Hall, Gabe Hassler, Joseph Heled, ",
                "\tSebastian Hoehna, Xiang Ji, Michael Karcher, Denise Kuehnert, Philippe Lemey, Wai Lok Sibon Li, ",
                "\tGerton Lunter, Andy Magee, Sidney Markowitz, JT McCrone, Volodymyr Minin, Julia Palacios, ",
                "\tMichael Defoin Platel, Oliver Pybus, Yucai Shao, Max Tolkoff, Chieh-Hsi Wu, Walter Xie and Zhenyu Zhang",
                "",
                "Thanks to:",
                "\tRoald Forsberg, Beth Shapiro and Korbinian Strimmer"};
    }

    public String getHTMLCredits() {
        return
                "<p>Designed and developed by<br>" +
                        "Alexei J. Drummond, Andrew Rambaut and Marc A. Suchard</p>" +
                        "<p>Department of Computer Science, University of Auckland<br>" +
                        "alexei@cs.auckland.ac.nz</p>" +
                        "<p>Institute of Ecology and Evolution, University of Edinburgh<br>" +
                        "a.rambaut@ed.ac.uk</p>" +
                        "<p>David Geffen School of Medicine, University of California, Los Angeles<br>" +
                        "msuchard@ucla.edu</p>" +
                        "<p>BEAST website and documentation:<br>" +
                        "<em>http://beast.community</em></p>" +
                        "<p>Source code distributed under the GNU LGPL:<br>" +
                        "<em>http://github.com/beast-dev/beast-mcmc</em></p>" +
                        "<p>BEAST developers:<br>" +
                        "Alex Alekseyenko, Daniel Ayres, Guy Baele, Trevor Bedford, Filip Bielejec, Erik Bloomquist,<br>"+
                        "Luiz Max Carvalho, Gabriela Cybis, Mandev Gill, Matthew Hall, Gabe Hassler, Joseph Heled,<br>"+
                        "Sebastian Hoehna, Xiang Ji, Michael Karcher, Denise Kuehnert, Philippe Lemey, Wai Lok Sibon Li,<br>"+
                        "Gerton Lunter, Andy Magee, Sidney Markowitz, JT McCrone, Volodymyr Minin, Julia Palacios,<br>"+
                        "Michael Defoin Platel, Oliver Pybus, Yucai Shao, Max Tolkoff, Chieh-Hsi Wu, Walter Xie and Zhenyu Zhang</p>" +
                        "<p>Thanks to Roald Forsberg, Beth Shapiro and Korbinian Strimmer</p>";
    }

    public String getBuildString() {
        // I think having the tag release is more useful than the last commit...
    //    return "https://github.com/beast-dev/beast-mcmc/commit/" + getRevision();
        return "https://github.com/beast-dev/beast-mcmc/releases/tag/v" + getVersion();
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
        return Arrays.asList(CITATIONS);
    }

    public static Citation[] CITATIONS = new Citation[] {
            new Citation(
                    new Author[]{
                            new Author("G", "Baele"),
                            new Author("X", "Ji"),
                            new Author("GW", "Hassler"),
                            new Author("JT", "McCrone"),
                            new Author("Y", "Shao"),
                            new Author("Z", "Zhang"),
                            new Author("AJ", "Holbrook"),
                            new Author("P", "Lemey"),
                            new Author("AJ", "Drummond"),
                            new Author("A", "Rambaut"),
                            new Author("MA", "Suchard"),
                    },
                    "BEAST X for Bayesian phylogenetic, phylogeographic and phylodynamic inference",
                    2025,
                    "Nature Methods",
                    "",
                    "10.1038/s41592-025-02751-x"),
    };
    public static void main(String[] args) {
        System.out.println(getRevision());
    }
    public static String getRevision() {
        try {
            try (InputStream in = BeastVersion.class.getResourceAsStream("/revision.txt")) {

                if (in != null) {

                    List<String> lines = new BufferedReader(
                            new InputStreamReader(in, StandardCharsets.UTF_8))
                            .lines()
                            .collect(Collectors.toList());
                    String revision = lines.get(1); //"commit-dirty" -dirty is only output if there are uncommited changes
                    if (revision.endsWith("-dirty")) {
                        revision = revision.substring(0, revision.length() - "-dirty".length());
                    }
                    return revision;
                } else {
                    return "unknown";
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("No revision file found. Try running `ant revision` to make it");
        }
    }
}
