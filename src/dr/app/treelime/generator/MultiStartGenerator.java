/*
 * MultiStartGenerator.java
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

package dr.app.treelime.generator;

import dr.app.treelime.components.ComponentFactory;
import dr.app.treelime.options.PartitionTreeModel;
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.app.treelime.types.StartingTreeType;
import dr.evomodelxml.tree.search.TimeTreeSearchMergeParser;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.LinkedHashMap;
import java.util.Map;

/**
 * Independent starts of the topology search (plan 7 §7.5), as separate BEAST runs: start 0 from TreeLime's starting
 * tree, every other from a random starting tree (BEAST's -seed draws it) with its own search seed; then a merge
 * ({@code <timeTreeSearchMerge>}) picks the best and copies its outputs to {@code <stem>.best.*}. Separate processes
 * because BEAST's random numbers and parser state are global.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class MultiStartGenerator extends Generator {

    public static final String SCRIPT_FILE = "run_starts.sh";
    public static final String MERGE_FILE = "search.merge.xml";

    private final ComponentFactory[] components;

    public MultiStartGenerator(TreeLimeOptions options, ComponentFactory[] components) {
        super(options, components);
        this.components = components;
    }

    /**
     * Writes the starts' XML files, the merge and the run script into the directory.
     */
    public void generate(File directory, int starts) throws Exception {
        if (starts < 2) {
            throw new GeneratorException("Choose at least two starts.");
        }
        TopologySearchOptions search = options.maximizationOptions.getTopologySearch();
        if (!search.isOn()) {
            throw new GeneratorException("Independent starts are starts of the topology search: turn it on.");
        }
        if (!directory.isDirectory() && !directory.mkdirs()) {
            throw new IOException("Cannot create the directory " + directory);
        }
        final String stem = options.fileNameStem;
        final long seed = search.seed;
        final Map<PartitionTreeModel, StartingTreeType> startingTrees =
                new LinkedHashMap<PartitionTreeModel, StartingTreeType>();
        for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
            startingTrees.put(tree, tree.getStartingTreeType());
        }
        try {
            for (int k = 0; k < starts; ++k) {
                options.fileNameStem = TimeTreeSearchMergeParser.startStem(stem, k);
                search.seed = seed + k;
                for (PartitionTreeModel tree : options.getPartitionTreeModels()) {
                    tree.setStartingTreeType(k == 0 ? startingTrees.get(tree) : StartingTreeType.RANDOM);
                }
                new BeastGenerator(options, components).generateXML(new File(directory,
                        new File(options.fileNameStem).getName() + ".xml"));
            }
        } finally {
            options.fileNameStem = stem;
            search.seed = seed;
            for (Map.Entry<PartitionTreeModel, StartingTreeType> entry : startingTrees.entrySet()) {
                entry.getKey().setStartingTreeType(entry.getValue());
            }
        }
        writeMerge(directory, starts);
        writeScript(directory, starts);
    }

    private void writeMerge(File directory, int starts) throws IOException {
        PrintWriter out = new PrintWriter(new FileWriter(new File(directory, MERGE_FILE)));
        out.println("<?xml version=\"1.0\" standalone=\"yes\"?>");
        out.println("<!-- The best of the " + starts + " starts of " + options.fileNameStem + "'s topology search -->");
        out.println("<beast>");
        out.println("\t<" + TimeTreeSearchMergeParser.NAME + " " + TimeTreeSearchMergeParser.STEM + "=\"" +
                options.fileNameStem + "\" " + TimeTreeSearchMergeParser.STARTS + "=\"" + starts + "\"/>");
        out.println("</beast>");
        out.close();
    }

    private void writeScript(File directory, int starts) throws IOException {
        File script = new File(directory, SCRIPT_FILE);
        PrintWriter out = new PrintWriter(new FileWriter(script));
        try {
            int jobs = options.maximizationOptions.getBootstrapJobs();
            int threads = Math.max(0, options.maximizationOptions.bootstrapThreads);
            out.println("#!/bin/sh");
            out.println("# Runs the " + starts + " starts of " + options.fileNameStem + "'s topology search with BEAST, " +
                    jobs + " at a time (start 0 from TreeLime's starting tree, the others from random trees),");
            out.println("# then copies the best start's outputs to " + options.fileNameStem + ".best.* and lists them " +
                    "all in " + options.fileNameStem + ".starts.tsv.");
            out.println("# BEAST: the BEAST command if it is not 'beast'; JOBS: the starts run at once;");
            out.println("# THREADS: the threads each run may use (-threads and -beagle_threads; 0 lets BEAST decide).");
            out.println("cd \"$(dirname \"$0\")\" || exit 1");
            out.println("BEAST=\"${BEAST:-beast}\"");
            out.println("JOBS=\"${JOBS:-" + jobs + "}\"");
            out.println("THREADS=\"${THREADS:-" + threads + "}\"");
            out.println("THREAD_OPTIONS=\"\"");
            out.println("if [ \"$THREADS\" -gt 0 ]; then THREAD_OPTIONS=\"-threads $THREADS -beagle_threads $THREADS\"; fi");
            out.println("export BEAST THREAD_OPTIONS");
            out.println("# each start in its own BEAST, JOBS at a time; the seed draws a random start's tree");
            out.println("export STEM=\"" + options.fileNameStem + "\"");
            out.println("seq 0 " + (starts - 1) + " | xargs -n 1 -P \"$JOBS\" sh -c " +
                    "'xml=\"$STEM.start$1.xml\"; log=\"$STEM.start$1.log\"; " +
                    "$BEAST $THREAD_OPTIONS -overwrite -seed \"$(($1 + 1))\" \"$xml\" > \"$log\" 2>&1 || " +
                    "echo \"$xml failed: see $log\"' start");
            out.println();
            out.println("# the best start");
            out.println("$BEAST -overwrite " + MERGE_FILE + " > search.merge.log 2>&1 && " +
                    "echo \"Best start's outputs copied to $STEM.best.*; all starts in $STEM.starts.tsv\" || " +
                    "echo \"the merge failed: see search.merge.log\"");
        } finally {
            out.close();
        }
        script.setExecutable(true);
    }
}
