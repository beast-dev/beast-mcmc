/*
 * TimeTreeSearchPerturbationTest.java
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

package test.dr.app.treelime;

import dr.app.beast.BeastParser;
import dr.app.beast.BeastVersion;
import dr.app.treelime.generator.MultiStartGenerator;
import dr.app.treelime.options.TopologySearchOptions;
import dr.app.treelime.options.TreeLimeOptions;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.CandidateSet;
import dr.evomodel.tree.search.OperatorPerturbation;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.evomodel.tree.search.TopologyHash;
import dr.evomodelxml.tree.search.TimeTreeSearchMergeParser;
import dr.inference.model.Likelihood;
import dr.math.MathUtils;
import dr.xml.XMLObject;
import junit.framework.TestCase;

import java.io.File;
import java.io.FileReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import static test.dr.app.treelime.TimeTreeSearchTest.Setup;
import static test.dr.app.treelime.TimeTreeSearchTest.joint;
import static test.dr.app.treelime.TimeTreeSearchTest.parse;
import static test.dr.app.treelime.TimeTreeSearchTest.search;
import static test.dr.app.treelime.TimeTreeSearchTest.tree;

/**
 * Escape and multi-start of the topology search (plan 7 M5): perturbation by BEAST's own topology operators, the
 * candidate set, nStop, nniDelta, and independent starts as separate runs merged afterwards. Run with assertions on,
 * as TimeTreeSearchTest.
 */
public class TimeTreeSearchPerturbationTest extends TestCase {

    public TimeTreeSearchPerturbationTest(String name) {
        super(name);
    }

    // each operator alone keeps a valid time tree, keeps a constrained clade, and is kept at least
    // sometimes; the objective stays finite
    public void testOperators() throws Exception {
        for (String operator : new String[]{"nni", "narrow", "wide", "leap", "wb"}) {
            Map<String, XMLObject> objects = parse(new Setup().clade(new int[]{1, 2, 3}).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            Likelihood joint = joint(objects);
            TopologyConstraints constraints = new TopologyConstraints(tree, joint, null);
            MathUtils.setSeed(3);
            OperatorPerturbation perturbation = new OperatorPerturbation(tree, search(objects).getEvaluator(),
                    constraints, operator, 0.1);
            final String before = TreeUtils.uniqueNewick(tree, tree.getRoot());
            int kept = perturbation.perturb(10);
            assertTrue(operator, kept > 0);
            assertTrue(operator, constraints.areSatisfied(tree));
            assertTrue(operator, joint.getLogLikelihood() > Double.NEGATIVE_INFINITY);
            for (int n = 0; n < tree.getNodeCount(); ++n) {
                NodeRef node = tree.getNode(n);
                if (!tree.isRoot(node)) {
                    assertTrue(operator, tree.getNodeHeight(node) < tree.getNodeHeight(tree.getParent(node)));
                }
            }
            if (!operator.equals("narrow")) { // narrow exchange may only reject on a UPGMA tree
                assertFalse(operator, before.equals(TreeUtils.uniqueNewick(tree, tree.getRoot())));
            }
            // the state is consistent: a full recomputation gives the same F
            final double f = joint.getLogLikelihood();
            search(objects).getEvaluator().verify(f);
        }
    }

    // the search with perturbation: the best F never decreases, it stops on nStop, the candidate set is sorted and
    // distinct and its best is the final tree's topology, and the same seed gives the same log
    public void testPerturbedSearch() throws Exception {
        File directory = Files.createTempDirectory("perturbation").toFile();
        File candidatesFile = new File(directory, "run.candidates.trees");
        File treeFile = new File(directory, "run.search.trees");
        String attributes = "nStop=\"4\" seed=\"11\" candidatesFileName=\"" + candidatesFile.getAbsolutePath() +
                "\" treeFileName=\"" + treeFile.getAbsolutePath() + "\"";
        List<String> logs = new ArrayList<String>();
        for (int run = 0; run < 2; ++run) {
            Map<String, XMLObject> objects = parse(new Setup().relaxed(true).search(attributes).xml());
            if (objects == null) {
                return;
            }
            TreeModel tree = tree(objects);
            TimeTreeSearch search = search(objects);
            search.run();
            assertEquals(TimeTreeSearch.StopReason.N_STOP, search.getStopReason());
            assertTrue(search.getPerturbations() >= 4);

            double best = Double.NEGATIVE_INFINITY;
            boolean perturbRow = false;
            StringBuilder sb = new StringBuilder();
            for (String row : search.getLog().getRows()) {
                String[] fields = row.split("\t", -1);
                double fBest = Double.parseDouble(fields[9]);
                assertTrue(row, fBest >= best - 1E-9 * Math.abs(fBest));
                best = fBest;
                perturbRow |= fields[2].equals("perturb");
                sb.append(row, 0, row.lastIndexOf('\t')).append('\n');
            }
            assertTrue(perturbRow);
            logs.add(sb.toString());

            CandidateSet candidates = search.getCandidates();
            assertTrue(candidates.size() > 1);
            Set<Long> hashes = new HashSet<Long>();
            double previous = Double.POSITIVE_INFINITY;
            for (CandidateSet.Candidate candidate : candidates.getCandidates()) {
                assertTrue(hashes.add(candidate.hash));
                assertTrue(candidate.f <= previous);
                previous = candidate.f;
            }
            assertEquals(candidates.getCandidates().get(0).hash, TopologyHash.of(tree));
            assertTrue(search.getFinalValue() >= candidates.getCandidates().get(0).f - 1E-6);

            // nniDelta at every internal node but the root, in the tree file; the candidates file has every tree
            double[] delta = search.getNniDelta();
            for (int i = 0; i < tree.getInternalNodeCount(); ++i) {
                NodeRef node = tree.getInternalNode(i);
                assertEquals(tree.isRoot(node), Double.isNaN(delta[node.getNumber()]));
            }
            String trees = new String(Files.readAllBytes(treeFile.toPath()), "UTF-8");
            assertTrue(trees.contains("nniDelta="));
            dr.evolution.io.NexusImporter importer = new dr.evolution.io.NexusImporter(new FileReader(candidatesFile));
            assertEquals(candidates.size(), importer.importTrees(null).size());
            assertTrue(search.getReport().contains("Perturbations: "));
        }
        assertEquals(logs.get(0), logs.get(1));
    }

    // independent starts: start 0 from the starting tree, the others from random trees with their own seeds, run here
    // one after the other as the script would (BEAST's -seed k + 1), then merged: the best F is picked and its outputs
    // copied
    public void testStarts() throws Exception {
        TreeLimeOptions options = TreeLimeSearchGeneratorTest.dengue(false);
        if (options == null) {
            return;
        }
        options.maximizationOptions.printToScreen = false;
        options.maximizationOptions.getTopologySearch().mode = TopologySearchOptions.Mode.NNI;
        File directory = Files.createTempDirectory("starts").toFile();
        String stem = new File(directory, "dengue").getAbsolutePath();
        options.fileNameStem = stem;
        options.selectParameters();
        options.selectOperators();
        new MultiStartGenerator(options, TreeLimeHarness.COMPONENTS).generate(directory, 3);
        assertEquals(stem, options.fileNameStem); // restored
        String script = new String(Files.readAllBytes(new File(directory, MultiStartGenerator.SCRIPT_FILE).toPath()),
                "UTF-8");
        assertTrue(script.contains("-seed"));
        assertTrue(script.contains(MultiStartGenerator.MERGE_FILE));

        double[] values = new double[3];
        for (int k = 0; k < 3; ++k) {
            File xml = new File(stem + ".start" + k + ".xml");
            String text = new String(Files.readAllBytes(xml.toPath()), "UTF-8");
            assertEquals(k == 0, !text.contains("<coalescentSimulator"));
            assertTrue(text.contains("seed=\"" + (1 + k) + "\""));
            MathUtils.setSeed(k + 1);
            Map<String, XMLObject> objects = new BeastParser(new String[]{}, null, false, false, false,
                    new BeastVersion()).parse(new java.io.StringReader(text), false);
            values[k] = ((TimeTreeSearch) objects.get("treeSearch").getNativeObject()).getFinalValue();
        }
        int best = TimeTreeSearchMergeParser.merge(stem, 3);
        for (int k = 0; k < 3; ++k) {
            assertTrue(values[best] >= values[k]);
        }
        List<String> rows = Files.readAllLines(new File(stem + ".starts.tsv").toPath());
        assertEquals(4, rows.size());
        assertTrue(new File(stem + ".best.search.trees").exists());
        assertTrue(new File(stem + ".best.timetree.nex").exists());
        assertEquals(new String(Files.readAllBytes(new File(stem + ".start" + best + ".search.trees").toPath())),
                new String(Files.readAllBytes(new File(stem + ".best.search.trees").toPath())));
    }
}
