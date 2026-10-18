/*
 * EpochSpectralPreOrderReportTest.java
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

package test.dr.evomodel.treedatalikelihood;

import beagle.Beagle;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchmodel.EpochBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EpochEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.DiscretePreOrderReport;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Parameter;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Checks that DiscretePreOrderReport's BOTTOM_SPECTRAL/TOP_SPECTRAL display rotates each node's raw pre-order
 * partial into the eigenbasis of THAT node's own segment (the substitution model on the branch immediately above
 * it), rather than always the first substitution model -- which matters for any multi-epoch model, and is the only
 * choice that is even well defined once a branch may carry several segments (augmented degree-2 nodes), each in a
 * different epoch.
 * <p>
 * Runs against the standard (non-spectral) BEAGLE CPU implementation: the eigenbasis rotation is plain Java math,
 * independent of which BEAGLE CPU implementation computed the raw (un-rotated) partial that it rotates.
 *
 * @author Marc A Suchard
 */
public class EpochSpectralPreOrderReportTest extends TraceCorrelationAssert {

    private static final int STATES = 4;

    // "node <id>[ augmented, branch <b>, epoch <e>][ root][ taxon=...]" header lines
    private static final Pattern NODE_HEADER = Pattern.compile("^node (\\d+)[^\\n]*$", Pattern.MULTILINE);

    // "  category <c> pattern <p> start=[..] end=[..]" data rows
    private static final Pattern DATA_ROW =
            Pattern.compile("  category (\\d+) pattern (\\d+) start=\\[(.*?)\\] end=\\[(.*?)\\]");

    public EpochSpectralPreOrderReportTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    public void testBottomSpectralRotatesPerSegmentModel() {
        checkRotation(DiscretePartialsType.BOTTOM_SPECTRAL, true);
    }

    public void testTopSpectralRotatesPerSegmentModel() {
        checkRotation(DiscretePartialsType.TOP_SPECTRAL, false);
    }

    private void checkRotation(DiscretePartialsType displayType, boolean bottom) {
        for (int seed = 0; seed < 3; ++seed) {

            Fixture f = new Fixture(seed);
            DiscretePreOrderReport report = new DiscretePreOrderReport(f.likelihood, displayType, 0.0);
            String text = report.getReport();

            final Beagle beagle = f.delegate.getBeagleInstance();
            final double[] raw = new double[f.categories * f.patternCount * STATES];

            // header positions, in order, so each header's block runs up to the next header (or end of text)
            List<Integer> headerStarts = new ArrayList<Integer>();
            List<Integer> headerIds = new ArrayList<Integer>();
            Matcher header = NODE_HEADER.matcher(text);
            while (header.find()) {
                headerStarts.add(header.end());
                headerIds.add(Integer.parseInt(header.group(1)));
            }

            int checked = 0;
            int rowsChecked = 0;
            int augmentedChecked = 0;
            for (int h = 0; h < headerIds.size(); ++h) {
                int id = headerIds.get(h);
                int blockEnd = (h + 1 < headerStarts.size()) ? headerStarts.get(h + 1) : text.length();
                String block = text.substring(headerStarts.get(h), blockEnd);

                beagle.getPartials(f.delegate.getPreOrderPartialIndex(id), Beagle.NONE, raw); // the raw (un-rotated) partial

                final boolean isRoot = (id == f.tree.getRoot().getNumber());
                double[] rotation = null;
                if (!isRoot) {
                    EigenDecomposition ed = f.delegate.getEvolutionaryProcessDelegate()
                            .getSubstitutionModelForBranch(id).getEigenDecomposition();
                    EigenDecomposition edT = ed.transpose();
                    rotation = bottom ? edT.getEigenVectors() : edT.getInverseEigenVectors();

                    if (f.registry.isAugmented(id)) {
                        ++augmentedChecked;
                    }
                }

                Matcher row = DATA_ROW.matcher(block);
                int rowsForNode = 0;
                while (row.find()) {
                    int c = Integer.parseInt(row.group(1));
                    int p = Integer.parseInt(row.group(2));
                    // getReport() always prints a "start=" (requested as TOP) and "end=" (requested as BOTTOM)
                    // column; only the one matching beagleType was actually computed -- the other is zero-filled
                    double[] reported = parseVector(bottom ? row.group(4) : row.group(3));

                    int offset = (c * f.patternCount + p) * STATES;
                    double[] expected = new double[STATES];
                    if (isRoot) {
                        // no branch above the root: the report leaves its partial (the root state
                        // frequencies) unrotated
                        System.arraycopy(raw, offset, expected, 0, STATES);
                    } else {
                        for (int i = 0; i < STATES; ++i) {
                            double sum = 0.0;
                            for (int j = 0; j < STATES; ++j) {
                                sum += rotation[i * STATES + j] * raw[offset + j];
                            }
                            expected[i] = sum;
                        }
                    }

                    assertEquals("node " + id + ", category " + c + ", pattern " + p, STATES, reported.length);
                    for (int s = 0; s < STATES; ++s) {
                        assertEquals("node " + id + ", category " + c + ", pattern " + p + ", state " + s,
                                expected[s], reported[s], 1E-9 * Math.max(1.0, Math.abs(expected[s])));
                    }
                    ++rowsForNode;
                    ++rowsChecked;
                }
                assertEquals("node " + id, f.categories * f.patternCount, rowsForNode);
                ++checked;
            }

            // every original node but the root, plus every augmented node, appears in the report
            int expectedCount = f.tree.getNodeCount();
            for (int n = 0; n < f.tree.getNodeCount(); ++n) {
                if (!f.tree.isRoot(f.tree.getNode(n))) {
                    expectedCount += f.registry.getChainLength(n);
                }
            }
            assertEquals("seed " + seed, expectedCount, checked);
            assertTrue("seed " + seed + ": no augmented node was exercised", augmentedChecked > 0);
            assertTrue("seed " + seed + ": patternCount/categoryCount too small to exercise the row loop",
                    rowsChecked > checked); // i.e., more than one (category, pattern) row per node was checked
        }
    }

    /**
     * Confirms the fix actually matters: an augmented node whose own model differs from model 0 is NOT correctly
     * rotated by model 0's eigenvectors (the bug this test guards against -- DiscretePreOrderReport previously
     * always used substitutionModels.get(0), regardless of which epoch a node's segment belonged to).
     */
    public void testDifferentEpochsGiveDifferentRotations() {
        Fixture f = new Fixture(0);
        DiscretePreOrderReport report = new DiscretePreOrderReport(f.likelihood, DiscretePartialsType.BOTTOM_SPECTRAL,
                0.0);
        String text = report.getReport();

        final Beagle beagle = f.delegate.getBeagleInstance();
        final double[] raw = new double[f.categories * f.patternCount * STATES];

        // an augmented node whose epoch is not epoch 0 (model 0), if there is one
        int[] chain = new int[f.registry.getBoundaryCount()];
        int lateId = -1;
        for (int n = 0; n < f.tree.getNodeCount() && lateId < 0; ++n) {
            if (f.tree.isRoot(f.tree.getNode(n))) {
                continue;
            }
            int count = f.registry.copyChain(n, chain);
            for (int j = count - 1; j >= 0; --j) {
                if (f.registry.getMatrixEpoch(chain[j]) != 0) {
                    lateId = chain[j];
                    break;
                }
            }
        }
        assertTrue("fixture has no node past epoch 0", lateId >= 0);

        double[] reported = null;
        Matcher header = NODE_HEADER.matcher(text);
        while (header.find()) {
            if (Integer.parseInt(header.group(1)) == lateId) {
                Matcher row = DATA_ROW.matcher(text.substring(header.end()));
                row.find(); // category 0, pattern 0
                reported = parseVector(row.group(4)); // BOTTOM_SPECTRAL -> the "end=" column
                break;
            }
        }
        assertNotNull(reported);

        beagle.getPartials(f.delegate.getPreOrderPartialIndex(lateId), Beagle.NONE, raw);

        // model 0's rotation -- what the pre-fix code would have applied regardless of lateId's actual epoch
        EigenDecomposition wrongEd = f.models.get(0).getEigenDecomposition().transpose();
        double[] wrong = multiply(wrongEd.getEigenVectors(), raw);

        boolean allClose = true;
        for (int s = 0; s < STATES; ++s) {
            if (Math.abs(wrong[s] - reported[s]) > 1E-6) {
                allClose = false;
            }
        }
        assertFalse("model 0's rotation should not (in general) match the report's, which correctly uses " +
                "node " + lateId + "'s own epoch " + f.registry.getMatrixEpoch(lateId), allClose);
    }

    private static double[] parseVector(String csv) {
        String[] parts = csv.split(",\\s*");
        double[] out = new double[parts.length];
        for (int i = 0; i < parts.length; ++i) {
            out[i] = Double.parseDouble(parts[i]);
        }
        return out;
    }

    private static double[] multiply(double[] matrix, double[] vector) {
        double[] out = new double[STATES];
        for (int i = 0; i < STATES; ++i) {
            double sum = 0.0;
            for (int j = 0; j < STATES; ++j) {
                sum += matrix[i * STATES + j] * vector[j];
            }
            out[i] = sum;
        }
        return out;
    }

    private class Fixture {

        final DefaultTreeModel tree;
        final BeagleDataLikelihoodDelegate delegate;
        final TreeDataLikelihood likelihood;
        final AugmentedNodeRegistry registry;
        final List<SubstitutionModel> models = new ArrayList<SubstitutionModel>();

        final int categories;
        final int patternCount;

        Fixture(long seed) {
            Random random = new Random(seed);
            categories = 1 + random.nextInt(4);

            tree = new DefaultTreeModel(randomTree(random));

            final int epochCount = 3 + random.nextInt(4);
            final double[] times = new double[epochCount - 1];
            double time = 0.0;
            for (int i = 0; i < times.length; ++i) {
                time += 0.05 + 0.3 * random.nextDouble();
                times[i] = time;
            }
            Parameter epochTimes = new Parameter.Default(times);

            for (int i = 0; i < epochCount; ++i) {
                double[] frequencies = new double[STATES];
                double sum = 0.0;
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] = 0.5 + random.nextDouble();
                    sum += frequencies[j];
                }
                for (int j = 0; j < STATES; ++j) {
                    frequencies[j] /= sum;
                }
                // distinct kappas, so different epochs really do have different (non-proportional) eigenbases
                models.add(new HKY(0.5 + 4.0 * random.nextDouble(), new FrequencyModel(Nucleotides.INSTANCE,
                        frequencies)));
            }

            EpochBranchModel branchModel = new EpochBranchModel(tree, models, epochTimes);
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            SitePatterns patterns = new SitePatterns(alignment, null, 0, -1, 1, true);
            patternCount = patterns.getPatternCount();

            // pre-order, standard (non-spectral) BEAGLE representation, augmented epoch nodes
            delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel, siteRateModel,
                    false, false, PartialsRescalingScheme.NONE, false,
                    new PreOrderSettings(true, true, false, false, false, false, true));

            likelihood = new TreeDataLikelihood(delegate, tree, new DefaultBranchRateModel());
            registry = ((EpochEvolutionaryProcessDelegate) delegate.getEvolutionaryProcessDelegate())
                    .getAugmentedNodeRegistry();
        }

        private SimpleTree randomTree(Random random) {
            List<SimpleNode> active = new ArrayList<SimpleNode>();
            for (int i = 0; i < taxa.length; ++i) {
                SimpleNode tip = new SimpleNode();
                tip.setTaxon(taxa[i]);
                tip.setHeight(0.0);
                active.add(tip);
            }

            while (active.size() > 1) {
                SimpleNode left = active.remove(random.nextInt(active.size()));
                SimpleNode right = active.remove(random.nextInt(active.size()));

                SimpleNode parent = new SimpleNode();
                parent.setHeight(Math.max(left.getHeight(), right.getHeight()) + 0.02 + 0.2 * random.nextDouble());
                parent.addChild(left);
                parent.addChild(right);
                active.add(parent);
            }

            return new SimpleTree(active.get(0));
        }
    }
}
