/*
 * TimeTreeSearchMergeParser.java
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

package dr.evomodelxml.tree.search;

import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * {@code <timeTreeSearchMerge>}: after the independent starts of a topology search (plan 7 §7.5), each a TreeLime run
 * in its own BEAST process, picks the start with the highest F (read from its .search.trees), writes a table of every
 * start's F and RF distances to the best ({@code <stem>.starts.tsv}), and copies the best start's outputs
 * ({@code <stem>.start<k>.*}) to {@code <stem>.best.*}. A start whose tree file is missing (it failed) is skipped.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeSearchMergeParser extends AbstractXMLObjectParser {

    public static final String NAME = "timeTreeSearchMerge";
    public static final String STEM = "stem";
    public static final String STARTS = "starts";

    private static final Pattern F = Pattern.compile("&F=([-0-9.Ee+]+)");

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final String stem = xo.getStringAttribute(STEM);
        final int starts = xo.getIntegerAttribute(STARTS);
        try {
            return merge(stem, starts);
        } catch (IOException e) {
            throw new XMLParseException("Merging the starts of " + stem + " failed: " + e.getMessage());
        }
    }

    /**
     * @return the best start's index
     */
    public static int merge(String stem, int starts) throws IOException {
        final List<Integer> found = new ArrayList<Integer>();
        final List<Double> values = new ArrayList<Double>();
        final List<Tree> trees = new ArrayList<Tree>();
        for (int k = 0; k < starts; ++k) {
            final File file = new File(startStem(stem, k) + ".search.trees");
            if (!file.exists()) {
                Logger.getLogger("dr.evomodel").warning("No " + file + ": start " + k + " is left out");
                continue;
            }
            final String text = new String(Files.readAllBytes(file.toPath()), "UTF-8");
            final Matcher m = F.matcher(text);
            if (!m.find()) {
                throw new IOException(file + " has no F");
            }
            final NexusImporter importer = new NexusImporter(new FileReader(file));
            try {
                trees.add(importer.importNextTree());
            } catch (dr.evolution.io.Importer.ImportException e) {
                throw new IOException(file + ": " + e.getMessage());
            }
            found.add(k);
            values.add(Double.parseDouble(m.group(1)));
        }
        if (found.isEmpty()) {
            throw new IOException("no start has a tree file");
        }
        int best = 0;
        for (int i = 1; i < found.size(); ++i) {
            if (values.get(i) > values.get(best)) {
                best = i;
            }
        }

        final PrintWriter table = new PrintWriter(new File(stem + ".starts.tsv"));
        table.println("start\tF\tdF_to_best\trootedRF_to_best\tunrootedRF_to_best\tbest");
        for (int i = 0; i < found.size(); ++i) {
            table.println(found.get(i) + "\t" + values.get(i) + "\t" + (values.get(i) - values.get(best)) + "\t" +
                    SplitDistance.rooted(trees.get(best), trees.get(i)) + "\t" +
                    SplitDistance.unrooted(trees.get(best), trees.get(i)) + "\t" + (i == best));
        }
        table.close();

        // the best start's outputs as <stem>.best.*
        final String prefix = new File(startStem(stem, found.get(best))).getName() + ".";
        File directory = new File(stem).getAbsoluteFile().getParentFile();
        File[] files = directory.listFiles();
        if (files != null) {
            for (File file : files) {
                if (file.getName().startsWith(prefix) && !file.getName().endsWith(".xml")) {
                    File copy = new File(directory, new File(stem).getName() + ".best." +
                            file.getName().substring(prefix.length()));
                    Files.copy(file.toPath(), copy.toPath(), StandardCopyOption.REPLACE_EXISTING);
                }
            }
        }
        Logger.getLogger("dr.evomodel").info("Best of " + found.size() + " starts: start " + found.get(best) +
                ", F = " + values.get(best) + "; outputs copied to " + stem + ".best.*, all starts in " + stem +
                ".starts.tsv");
        return found.get(best);
    }

    /**
     * @return the file-name stem of start k
     */
    public static String startStem(String stem, int k) {
        return stem + ".start" + k;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(STEM, "the analysis's file-name stem; start k wrote <stem>.start<k>.*"),
            AttributeRule.newIntegerRule(STARTS, false, "the number of starts"),
    };

    @Override
    public String getParserDescription() {
        return "Picks the best of a topology search's independent starts and copies its outputs to <stem>.best.*";
    }

    @Override
    public Class getReturnType() {
        return Integer.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }
}
