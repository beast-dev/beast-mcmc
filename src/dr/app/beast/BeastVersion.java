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

    private static final String DATE_STRING = "2002-2024";

    private static final boolean IS_PRERELEASE = true;

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
                "\tAlex Alekseyenko, Daniel Ayres, Guy Baele, Trevor Bedford, Filip Bielejec, Erik Bloomquist, ",
                "\tLuiz Max Carvalho, Gabriela Cybis, Mandev Gill, Matthew Hall, Gabe Hassler, Joseph Heled, ",
                "\tSebastian Hoehna, Xiang Ji, Michael Karcher, Denise Kuehnert, Philippe Lemey, Wai Lok Sibon Li, ",
                "\tGerton Lunter, Andy Magee, Sidney Markowitz, JT McCrone, Volodymyr Minin, Julia Palacios, ",
                "\tMichael Defoin Platel, Oliver Pybus, Yucai Shao, Max Tolkoff, Chieh-Hsi Wu, Walter Xie and ",
                "\tZhenyu Zhang",
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
                        "Michael Defoin Platel, Oliver Pybus, Yucai Shao, Max Tolkoff, Chieh-Hsi Wu, Walter Xie and<br>"+
                        "Zhenyu Zhang</p>" +
                        "<p>Thanks to Roald Forsberg, Beth Shapiro and Korbinian Strimmer</p>";
    }

    public String getBuildString() {
        return "https://github.com/beast-dev/beast-mcmc/commit/" + getRevision();
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
                            new Author("MA", "Suchard"),
                            new Author("P", "Lemey"),
                            new Author("G", "Baele"),
                            new Author("DL", "Ayres"),
                            new Author("AJ", "Drummond"),
                            new Author("A", "Rambaut")
                    },
                    "Bayesian phylogenetic and phylodynamic data integration using BEAST 1.10",
                    2018,
                    "Virus Evolution",
                    4, "vey016",
                    "10.1093/ve/vey016"),
//            new Citation(
//                    new Author[]{
//                            new Author("AJ", "Drummond"),
//                            new Author("MA", "Suchard"),
//                            new Author("Dong", "Xie"),
//                            new Author("A", "Rambaut")
//                    },
//                    "Bayesian phylogenetics with BEAUti and the BEAST 1.7",
//                    2012,
//                    "Mol Biol Evol",
//                    29, 1969, 1973,
//                    "10.1093/molbev/mss075")
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
                    return lines.get(1); //"commit-dirty" -dirty is only output if there are uncommited changes
                } else {
                    return "unknown";
                }
            }
        } catch (IOException e) {
            throw new RuntimeException("No revision file found. Try running `ant revision` to make it");
        }
    }
}
