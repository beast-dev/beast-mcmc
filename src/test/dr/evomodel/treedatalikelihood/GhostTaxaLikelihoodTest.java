/*
 * GhostTaxaLikelihoodTest.java
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

import beagle.BeagleInfo;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.ExternalInternalBranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.AbstractBranchRateModel;
import dr.evolution.tree.Tree;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.BeagleFunctionality;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import test.dr.inference.trace.TraceCorrelationAssert;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;

/**
 * Taxa in the tree that a partition has no sequences for are "ghosts" in that partition. Their likelihood must equal
 * that of the same analysis with the ghosts' sequences written as all N, for every arrangement of the ghosts in the
 * tree, whichever way a node with one ghost child is computed, and through a run of proposals that move ghosts into
 * and out of ghost subtrees and are accepted or rejected.
 *
 * @author Andrew Rambaut
 */
public class GhostTaxaLikelihoodTest extends TraceCorrelationAssert {

    private static final double TOLERANCE = 1E-9;

    private enum Arrangement { SINGLES, CHERRY, CLADE, ROOT_SISTER, ONE_REAL, TWO_REAL, MIXED, HALF }

    public GhostTaxaLikelihoodTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
        createAlignment(DENGUE4_TAXON_SEQUENCE, Nucleotides.INSTANCE);
    }

    @Override
    protected void tearDown() throws Exception {
        System.clearProperty(BeagleDataLikelihoodDelegate.GHOST_ONE_CHILD_PROPERTY);
        System.clearProperty(BeagleDataLikelihoodDelegate.GHOST_MERGE_PROPERTY);
        super.tearDown();
    }

    /**
     * The all-unknown tip that stands in for a ghost subtree works with any BEAGLE.
     */
    public void testWithAnAllUnknownTip() {
        System.setProperty(BeagleDataLikelihoodDelegate.GHOST_ONE_CHILD_PROPERTY, "false");
        compareThroughProposals(new PartialsRescalingScheme[]{PartialsRescalingScheme.NONE,
                PartialsRescalingScheme.DYNAMIC, PartialsRescalingScheme.ALWAYS, PartialsRescalingScheme.AUTO});
    }

    /**
     * One-child operations need native BEAGLE 4.2.0 or later; with auto rescaling the delegate does not use them.
     */
    public void testWithOneChildOperations() {
        if (!BeagleFunctionality.IS_MATRIX_VECTOR_SUPPORTED()) {
            System.err.println("SKIPPED " + getClass().getSimpleName() + "." + getName() +
                    ": needs native BEAGLE 4.2.0 or later, found " + BeagleInfo.getVersion());
            return;
        }
        compareThroughProposals(new PartialsRescalingScheme[]{PartialsRescalingScheme.NONE,
                PartialsRescalingScheme.DYNAMIC, PartialsRescalingScheme.ALWAYS});
    }

    /**
     * Without joining the branches through nodes with one ghost child, each of those nodes is computed.
     */
    public void testWithAnAllUnknownTipWithoutJoining() {
        System.setProperty(BeagleDataLikelihoodDelegate.GHOST_MERGE_PROPERTY, "false");
        testWithAnAllUnknownTip();
    }

    public void testWithOneChildOperationsWithoutJoining() {
        System.setProperty(BeagleDataLikelihoodDelegate.GHOST_MERGE_PROPERTY, "false");
        testWithOneChildOperations();
    }

    /**
     * A branch model with a rate matrix for external branches and another for internal ones: the branches through a
     * node with one ghost child cannot be joined, so the delegate computes those nodes instead.
     */
    public void testTwoSubstitutionModels() {
        assertNull(BeagleDataLikelihoodDelegate.ghostMergeProblem(
                new HomogeneousBranchModel(new HKY(2.0, new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.25, 0.25, 0.25, 0.25}))), PreOrderSettings.getDefault()));
        externalInternal = true;
        try {
            compareThroughProposals(new PartialsRescalingScheme[]{PartialsRescalingScheme.DYNAMIC,
                    PartialsRescalingScheme.ALWAYS});
        } finally {
            externalInternal = false;
        }
    }

    /**
     * A move that touches ghosts alone changes nothing that is computed, so it must cost no partials operations and
     * leave the likelihood as it was; proposals of every kind after it, accepted or rejected, must still match all N,
     * and a full evaluation must recompute every computed node.
     */
    public void testMovesAmongGhostsComputeNothing() {
        for (boolean oneChild : new boolean[]{true, false}) {
            if (oneChild && !BeagleFunctionality.IS_MATRIX_VECTOR_SUPPORTED()) {
                continue;
            }
            System.setProperty(BeagleDataLikelihoodDelegate.GHOST_ONE_CHILD_PROPERTY, Boolean.toString(oneChild));
            for (Arrangement arrangement : Arrangement.values()) {
                for (PartialsRescalingScheme scheme : new PartialsRescalingScheme[]{PartialsRescalingScheme.NONE,
                        PartialsRescalingScheme.DYNAMIC, PartialsRescalingScheme.ALWAYS}) {
                    final String label = arrangement + ", " + scheme.getText() + (oneChild ? ", one child" : ", tip");
                    Fixture ghosts = new Fixture(0, arrangement, true, scheme, 4, false);
                    Fixture allN = new Fixture(0, arrangement, false, scheme, 4, false);
                    assertSame(label, allN, ghosts);

                    int ghostMoves = 0;
                    for (int iteration = 0; iteration < 80; ++iteration) {
                        final double accepted = allN.likelihood.getLogLikelihood();
                        ghosts.store();
                        allN.store();

                        if (iteration % 2 == 1) {
                            final long before = ghosts.likelihood.getTotalCalculationCount();
                            final boolean moved = ghosts.proposeAmongGhosts();
                            assertEquals(moved, allN.proposeAmongGhosts());
                            if (moved) {
                                ghostMoves++;
                                assertEquals(label + ", iteration " + iteration + ", unchanged", accepted,
                                        ghosts.likelihood.getLogLikelihood(), TOLERANCE * Math.abs(accepted));
                                assertEquals(label + ", iteration " + iteration + ", partials computed", before,
                                        ghosts.likelihood.getTotalCalculationCount());
                            }
                        } else {
                            ghosts.propose();
                            allN.propose();
                        }
                        assertSame(label + ", iteration " + iteration, allN, ghosts);

                        if (iteration % 3 == 0) {
                            ghosts.likelihood.acceptModelState();
                            allN.likelihood.acceptModelState();
                        } else {
                            ghosts.likelihood.restoreModelState();
                            allN.likelihood.restoreModelState();
                            assertEquals(label + ", restored", accepted, ghosts.likelihood.getLogLikelihood(),
                                    TOLERANCE * Math.abs(accepted));
                        }
                    }
                    if (ghosts.realTipCount() > 1 && ghosts.ghostIds.size() > 1) {
                        assertTrue(label + ": no move among ghosts was possible", ghostMoves > 0);
                    }

                    // a full evaluation recomputes every computed node: at least one per real tip less one
                    final long before = ghosts.likelihood.getTotalCalculationCount();
                    ghosts.likelihood.makeDirty();
                    allN.likelihood.makeDirty();
                    assertSame(label + ", full evaluation", allN, ghosts);
                    assertTrue(label + ", full evaluation recomputes", ghosts.likelihood.getTotalCalculationCount() -
                            before >= ghosts.realTipCount() - 1);
                }
            }
        }
    }

    private boolean externalInternal = false;

    private void compareThroughProposals(PartialsRescalingScheme[] schemes) {
        for (Arrangement arrangement : Arrangement.values()) {
            for (PartialsRescalingScheme scheme : schemes) {
                for (boolean useAmbiguities : new boolean[]{false, true}) {
                    for (int seed = 0; seed < 2; ++seed) {
                        final int categories = (seed == 0) ? 1 : 4;
                        final String label = arrangement + ", " + scheme.getText() + ", " +
                                (useAmbiguities ? "ambiguities" : "states") + ", " + categories + " categories";

                        Fixture ghosts = new Fixture(seed, arrangement, true, scheme, categories, useAmbiguities);
                        Fixture allN = new Fixture(seed, arrangement, false, scheme, categories, useAmbiguities);

                        assertSame(label, allN, ghosts);

                        for (int iteration = 0; iteration < 60; ++iteration) {
                            final double accepted = allN.likelihood.getLogLikelihood();

                            ghosts.store();
                            allN.store();

                            ghosts.propose();
                            allN.propose();

                            assertSame(label + ", iteration " + iteration, allN, ghosts);

                            if (iteration % 3 == 0) {
                                ghosts.likelihood.acceptModelState();
                                allN.likelihood.acceptModelState();
                            } else {
                                ghosts.likelihood.restoreModelState();
                                allN.likelihood.restoreModelState();

                                assertEquals(label + ", restored", accepted, ghosts.likelihood.getLogLikelihood(),
                                        TOLERANCE * Math.abs(accepted));
                            }
                        }
                    }
                }
            }
        }
    }

    private void assertSame(String label, Fixture expected, Fixture actual) {
        final double logL = expected.likelihood.getLogLikelihood();
        assertEquals(label, logL, actual.likelihood.getLogLikelihood(), TOLERANCE * Math.abs(logL));
    }

    // ---- the XML: two partitions with different ghosts, and the errors -------------------------------------------

    /**
     * One tree, two partitions each missing different taxa, against the same partitions with those taxa as all N.
     */
    public void testTwoPartitionsWithDifferentGhosts() throws Exception {
        Set<String> missing1 = setOf(taxa[0].getId(), taxa[1].getId(), taxa[2].getId());
        Set<String> missing2 = setOf(taxa[5].getId(), taxa[9].getId());

        double ghosts = logLikelihoodOf(xml(missing1, missing2, false, ""));
        double allN = logLikelihoodOf(xml(missing1, missing2, true, ""));

        assertEquals(allN, ghosts, TOLERANCE * Math.abs(allN));
    }

    /**
     * A taxon that a partition has and the tree does not is still an error.
     */
    public void testATaxonNotInTheTreeIsAnError() throws Exception {
        final String sequence = alignment.getSequence(0).getSequenceString().substring(0, alignment.getSiteCount() / 2);
        String xml = xml(setOf(taxa[0].getId()), setOf(), false, "")
                .replace("</taxa>", "<taxon id=\"outsider\"/>\n</taxa>")
                .replace("<alignment id=\"gene1\" dataType=\"nucleotide\">\n",
                        "<alignment id=\"gene1\" dataType=\"nucleotide\">\n<sequence><taxon idref=\"outsider\"/>" +
                                sequence + "</sequence>\n");
        try {
            logLikelihoodOf(xml);
            fail("A partition with a taxon that is not in the tree was accepted");
        } catch (XMLParseException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("outsider"));
        }
    }

    public void testGhostsWithPreOrderAreAnError() throws Exception {
        try {
            logLikelihoodOf(xml(setOf(taxa[0].getId()), setOf(), false,
                    " usePreOrder=\"true\""));
            fail("Ghost taxa were accepted with pre-order partials");
        } catch (XMLParseException e) {
            assertTrue(e.getMessage(), e.getMessage().contains("usePreOrder"));
        }
    }

    private static Set<String> setOf(String... ids) {
        Set<String> set = new HashSet<String>();
        for (String id : ids) {
            set.add(id);
        }
        return set;
    }

    // the parsers the XML needs, as BEAST registers them (BeastParser itself reads its list from the jar)
    private static final String[] PARSERS = {
            "dr.inference.model.ParameterParser",
            "dr.evoxml.TaxonParser",
            "dr.evoxml.TaxaParser",
            "dr.evoxml.SequenceParser",
            "dr.evoxml.AlignmentParser",
            "dr.evoxml.SitePatternsParser",
            "dr.evoxml.NewickParser",
            "dr.evomodelxml.tree.TreeModelParser",
            "dr.evomodelxml.substmodel.FrequencyModelParser",
            "dr.evomodelxml.substmodel.HKYParser",
            "dr.evomodelxml.siteratemodel.OldGammaSiteModelParser",
            "dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser",
            "dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser"
    };

    private double logLikelihoodOf(String xml) throws Exception {
        XMLParser parser = new XMLParser(true, true, true, null);
        for (String name : PARSERS) {
            parser.addXMLObjectParser((XMLObjectParser) Class.forName(name).getDeclaredConstructor().newInstance());
        }
        Map<String, XMLObject> objects = parser.parse(new StringReader(xml), true);
        Likelihood likelihood = (Likelihood) objects.get("treeLikelihood").getNativeObject();
        return likelihood.getLogLikelihood();
    }

    /**
     * Two partitions of the alignment's first and second halves. Each leaves out its missing taxa, or keeps them with
     * all-N sequences.
     */
    private String xml(Set<String> missing1, Set<String> missing2, boolean asAllN, String likelihoodAttributes) {
        final int half = alignment.getSiteCount() / 2;

        StringBuilder sb = new StringBuilder("<beast>\n<taxa id=\"taxa\">\n");
        for (Taxon taxon : taxa) {
            sb.append("<taxon id=\"").append(taxon.getId()).append("\"/>\n");
        }
        sb.append("</taxa>\n");

        appendAlignment(sb, "gene1", 0, half, missing1, asAllN);
        appendAlignment(sb, "gene2", half, alignment.getSiteCount(), missing2, asAllN);

        sb.append("<patterns id=\"gene1.patterns\" from=\"1\"><alignment idref=\"gene1\"/></patterns>\n");
        sb.append("<patterns id=\"gene2.patterns\" from=\"1\"><alignment idref=\"gene2\"/></patterns>\n");

        sb.append("<newick id=\"startingTree\" usingDates=\"false\">").append(newick()).append("</newick>\n");
        sb.append("<treeModel id=\"treeModel\"><newick idref=\"startingTree\"/>" +
                "<rootHeight><parameter id=\"treeModel.rootHeight\"/></rootHeight>" +
                "<nodeHeights internalNodes=\"true\"><parameter id=\"treeModel.internalNodeHeights\"/></nodeHeights>" +
                "<nodeHeights internalNodes=\"true\" rootNode=\"true\">" +
                "<parameter id=\"treeModel.allInternalNodeHeights\"/></nodeHeights></treeModel>\n");

        sb.append("<HKYModel id=\"hky\"><frequencies><frequencyModel dataType=\"nucleotide\"><frequencies>" +
                "<parameter id=\"frequencies\" value=\"0.3 0.2 0.2 0.3\"/></frequencies></frequencyModel>" +
                "</frequencies><kappa><parameter id=\"kappa\" value=\"3.0\"/></kappa></HKYModel>\n");
        sb.append("<siteModel id=\"siteModel\"><substitutionModel><HKYModel idref=\"hky\"/></substitutionModel>" +
                "<gammaShape gammaCategories=\"4\"><parameter id=\"alpha\" value=\"0.5\"/></gammaShape></siteModel>\n");

        sb.append("<treeDataLikelihood id=\"treeLikelihood\" useAmbiguities=\"false\"")
                .append(likelihoodAttributes).append(">\n");
        sb.append("<partition><patterns idref=\"gene1.patterns\"/><siteModel idref=\"siteModel\"/></partition>\n");
        sb.append("<partition><patterns idref=\"gene2.patterns\"/><siteModel idref=\"siteModel\"/></partition>\n");
        sb.append("<treeModel idref=\"treeModel\"/>\n");
        sb.append("<strictClockBranchRates id=\"branchRates\"><rate><parameter id=\"clock.rate\" value=\"1.0\"/>" +
                "</rate></strictClockBranchRates>\n");
        sb.append("</treeDataLikelihood>\n</beast>\n");
        return sb.toString();
    }

    private void appendAlignment(StringBuilder sb, String id, int from, int to, Set<String> missing, boolean asAllN) {
        sb.append("<alignment id=\"").append(id).append("\" dataType=\"nucleotide\">\n");
        for (int i = 0; i < alignment.getSequenceCount(); i++) {
            String taxon = alignment.getTaxonId(i);
            String sequence = alignment.getSequence(i).getSequenceString().substring(from, to);
            if (missing.contains(taxon)) {
                if (!asAllN) {
                    continue;
                }
                sequence = sequence.replaceAll(".", "N");
            }
            sb.append("<sequence><taxon idref=\"").append(taxon).append("\"/>").append(sequence)
                    .append("</sequence>\n");
        }
        sb.append("</alignment>\n");
    }

    private String newick() {
        SimpleTree tree = randomTree(new Random(7));
        return newick(tree, tree.getRoot()) + ";";
    }

    private static String newick(SimpleTree tree, NodeRef node) {
        double length = tree.isRoot(node) ? 0.0 : tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
        if (tree.isExternal(node)) {
            return tree.getTaxonId(node.getNumber()) + ":" + length;
        }
        return "(" + newick(tree, tree.getChild(node, 0)) + "," + newick(tree, tree.getChild(node, 1)) + "):" +
                length;
    }

    private SimpleTree randomTree(Random random) {
        List<SimpleNode> active = new ArrayList<SimpleNode>();
        for (Taxon taxon : taxa) {
            SimpleNode tip = new SimpleNode();
            tip.setTaxon(taxon);
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

    // ---- a likelihood over the alignment, with ghosts left out or written as all N --------------------------------

    /**
     * A tree over every taxon, ghosts chosen from it by the arrangement, and the likelihood of the alignment with the
     * ghosts left out of the patterns or kept as all N. Two fixtures with the same seed and arrangement have the same
     * tree and make the same proposals.
     */
    private class Fixture {

        final DefaultTreeModel tree;
        final TreeDataLikelihood likelihood;
        final Set<String> ghostIds;
        private final HKY hky;
        private final NodeRates rates;
        private final GammaSiteRateModel siteRateModel;

        private final Random proposals;

        Fixture(long seed, Arrangement arrangement, boolean leaveOut, PartialsRescalingScheme scheme,
                int categories, boolean useAmbiguities) {
            Random random = new Random(seed);
            proposals = new Random(seed + 1000);

            tree = new DefaultTreeModel(randomTree(random));

            Set<String> ghosts = chooseGhosts(arrangement, random);
            ghostIds = ghosts;

            SimpleAlignment data = new SimpleAlignment();
            data.setDataType(Nucleotides.INSTANCE);
            for (int i = 0; i < alignment.getSequenceCount(); i++) {
                Sequence original = alignment.getSequence(i);
                String id = original.getTaxon().getId();
                if (ghosts.contains(id)) {
                    if (leaveOut) {
                        continue;
                    }
                    Sequence allN = new Sequence(original.getSequenceString().replaceAll(".", "N"));
                    allN.setTaxon(original.getTaxon());
                    allN.setDataType(Nucleotides.INSTANCE);
                    data.addSequence(allN);
                } else {
                    data.addSequence(original);
                }
            }

            SitePatterns patterns = new SitePatterns(data, null, 0, -1, 1, true);

            // root frequencies that differ from the stationary ones, so that a root with one ghost child is only
            // right if its real child is carried up its own branch
            HKY hky = new HKY(2.0 + random.nextDouble(),
                    new FrequencyModel(Nucleotides.INSTANCE, new double[]{0.3, 0.2, 0.2, 0.3}));
            FrequencyModel rootFrequencies = new FrequencyModel(Nucleotides.INSTANCE,
                    new double[]{0.1, 0.4, 0.4, 0.1});
            this.hky = hky;
            final BranchModel branchModel;
            if (externalInternal) {
                HKY internal = new HKY(6.0, new FrequencyModel(Nucleotides.INSTANCE,
                        new double[]{0.2, 0.3, 0.3, 0.2}));
                branchModel = new ExternalInternalBranchModel(tree, hky, internal);
                assertNotNull(BeagleDataLikelihoodDelegate.ghostMergeProblem(branchModel,
                        PreOrderSettings.getDefault()));
            } else {
                branchModel = new HomogeneousBranchModel(hky, rootFrequencies);
            }
            siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            // a rate for each branch, so that the branches joined through a node with one ghost child differ in rate
            double[] initialRates = new double[tree.getNodeCount()];
            for (int i = 0; i < initialRates.length; i++) {
                initialRates[i] = 0.5 + random.nextDouble();
            }
            rates = new NodeRates(initialRates);

            BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel,
                    siteRateModel, useAmbiguities, false, scheme, false, PreOrderSettings.getDefault(), null, null,
                    leaveOut);

            likelihood = new TreeDataLikelihood(delegate, tree, rates);
        }

        private Set<String> chooseGhosts(Arrangement arrangement, Random random) {
            Set<String> ghosts = new HashSet<String>();
            final int n = tree.getExternalNodeCount();

            switch (arrangement) {
                case SINGLES:
                    while (ghosts.size() < 4) {
                        ghosts.add(tree.getTaxonId(random.nextInt(n)));
                    }
                    break;
                case CHERRY:
                    addTips(cladeOfSize(2, 2), ghosts);
                    break;
                case CLADE:
                    addTips(cladeOfSize(3, n - 3), ghosts);
                    break;
                case ROOT_SISTER: {
                    NodeRef root = tree.getRoot();
                    NodeRef left = tree.getChild(root, 0);
                    NodeRef right = tree.getChild(root, 1);
                    addTips(tipCount(left) <= tipCount(right) ? left : right, ghosts);
                    break;
                }
                case ONE_REAL:
                    for (int i = 0; i < n; i++) {
                        ghosts.add(tree.getTaxonId(i));
                    }
                    ghosts.remove(tree.getTaxonId(random.nextInt(n)));
                    break;
                case TWO_REAL: {
                    // two real tips, with long runs of nodes with one ghost child above each
                    for (int i = 0; i < n; i++) {
                        ghosts.add(tree.getTaxonId(i));
                    }
                    ghosts.remove(tree.getTaxonId(random.nextInt(n)));
                    String other;
                    do {
                        other = tree.getTaxonId(random.nextInt(n));
                    } while (!ghosts.contains(other));
                    ghosts.remove(other);
                    break;
                }
                case HALF:
                    while (ghosts.size() < n / 2) {
                        ghosts.add(tree.getTaxonId(random.nextInt(n)));
                    }
                    break;
                case MIXED:
                    addTips(cladeOfSize(3, n - 4), ghosts);
                    while (ghosts.size() < n - 3) {
                        String id = tree.getTaxonId(random.nextInt(n));
                        ghosts.add(id);
                        if (ghosts.size() >= n / 2 + 2) {
                            break;
                        }
                    }
                    break;
            }
            return ghosts;
        }

        private NodeRef cladeOfSize(int min, int max) {
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                int size = tipCount(node);
                if (size >= min && size <= max) {
                    return node;
                }
            }
            throw new IllegalStateException("No clade of " + min + " to " + max + " tips");
        }

        private int tipCount(NodeRef node) {
            if (tree.isExternal(node)) {
                return 1;
            }
            return tipCount(tree.getChild(node, 0)) + tipCount(tree.getChild(node, 1));
        }

        private void addTips(NodeRef node, Set<String> ids) {
            if (tree.isExternal(node)) {
                ids.add(tree.getTaxonId(node.getNumber()));
            } else {
                addTips(tree.getChild(node, 0), ids);
                addTips(tree.getChild(node, 1), ids);
            }
        }

        void store() {
            likelihood.storeModelState();
        }

        /**
         * Exchanges two subtrees, which carries ghosts into and out of ghost subtrees, or moves a node height
         */
        void propose() {
            final int move = proposals.nextInt(8);
            if (move == 3) {
                // a new rate on one branch
                NodeRef node;
                do {
                    node = tree.getNode(proposals.nextInt(tree.getNodeCount()));
                } while (tree.isRoot(node));
                rates.setRate(node.getNumber(), 0.5 + proposals.nextDouble());
                return;
            }
            if (move == 4) {
                if (proposals.nextBoolean()) {
                    // a change to the rate matrix of every branch
                    hky.setKappa(1.5 + 3.0 * proposals.nextDouble());
                } else {
                    // a change to the rates of the site categories, so to the matrices of every branch
                    siteRateModel.setAlpha(0.2 + proposals.nextDouble());
                }
                return;
            }
            if (move == 5 && moveAboveTheRoot()) {
                return;
            }
            if (move == 6 && moveTheRoot()) {
                return;
            }
            if (move == 7 && proposeAmongGhosts()) {
                return;
            }
            if (move == 0) {
                NodeRef node;
                do {
                    node = tree.getInternalNode(proposals.nextInt(tree.getInternalNodeCount()));
                } while (tree.isRoot(node));

                double low = Math.max(tree.getNodeHeight(tree.getChild(node, 0)),
                        tree.getNodeHeight(tree.getChild(node, 1)));
                double high = tree.getNodeHeight(tree.getParent(node));
                tree.setNodeHeight(node, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
                return;
            }

            final int nodeCount = tree.getNodeCount();
            for (int attempt = 0; attempt < 1000; attempt++) {
                NodeRef i = tree.getNode(proposals.nextInt(nodeCount));
                NodeRef j = tree.getNode(proposals.nextInt(nodeCount));
                if (tree.isRoot(i) || tree.isRoot(j) || i == j) {
                    continue;
                }
                NodeRef iP = tree.getParent(i);
                NodeRef jP = tree.getParent(j);
                if (iP == jP || isAncestor(i, j) || isAncestor(j, i)) {
                    continue;
                }
                if (tree.getNodeHeight(j) >= tree.getNodeHeight(iP) ||
                        tree.getNodeHeight(i) >= tree.getNodeHeight(jP)) {
                    continue;
                }

                tree.beginTreeEdit();
                tree.removeChild(iP, i);
                tree.removeChild(jP, j);
                tree.addChild(jP, i);
                tree.addChild(iP, j);
                tree.endTreeEdit();
                return;
            }
        }

        /**
         * Wilson-Balding to above the root: a subtree with its parent is moved there, so its parent is the new root
         * and the old root is not the root any more.
         */
        private boolean moveAboveTheRoot() {
            final NodeRef root = tree.getRoot();
            for (int attempt = 0; attempt < 1000; attempt++) {
                NodeRef i = tree.getNode(proposals.nextInt(tree.getNodeCount()));
                if (tree.isRoot(i) || tree.isRoot(tree.getParent(i))) {
                    continue;
                }
                NodeRef iP = tree.getParent(i);
                NodeRef sibling = tree.getChild(iP, 0) == i ? tree.getChild(iP, 1) : tree.getChild(iP, 0);
                NodeRef grandparent = tree.getParent(iP);

                tree.beginTreeEdit();
                tree.setNodeHeight(iP, tree.getNodeHeight(root) + 0.02 + 0.2 * proposals.nextDouble());
                tree.removeChild(iP, sibling);
                tree.removeChild(grandparent, iP);
                tree.addChild(iP, root);
                tree.addChild(grandparent, sibling);
                tree.setRoot(iP);
                tree.endTreeEdit();
                return true;
            }
            return false;
        }

        /**
         * Wilson-Balding of the root: one of its children takes it down to a branch in the other child's subtree, and
         * the other child is the new root.
         */
        private boolean moveTheRoot() {
            final NodeRef root = tree.getRoot();
            final NodeRef i = tree.getChild(root, proposals.nextInt(2));
            final NodeRef other = tree.getChild(root, 0) == i ? tree.getChild(root, 1) : tree.getChild(root, 0);
            if (tree.isExternal(other)) {
                return false;
            }
            for (int attempt = 0; attempt < 1000; attempt++) {
                NodeRef j = tree.getNode(proposals.nextInt(tree.getNodeCount()));
                if (j == other || tree.isRoot(j) || !isAncestor(other, j)) {
                    continue;
                }
                NodeRef k = tree.getParent(j);
                double low = Math.max(tree.getNodeHeight(j), tree.getNodeHeight(i));
                double high = tree.getNodeHeight(k);
                if (low >= high) {
                    continue;
                }

                tree.beginTreeEdit();
                tree.setNodeHeight(root, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
                tree.removeChild(k, j);
                tree.removeChild(root, other);
                tree.addChild(root, j);
                tree.addChild(k, root);
                tree.setRoot(other);
                tree.endTreeEdit();
                return true;
            }
            return false;
        }

        int realTipCount() {
            return tree.getExternalNodeCount() - ghostIds.size();
        }

        /**
         * @return whether every tip below a node is a ghost
         */
        boolean isGhost(NodeRef node) {
            if (tree.isExternal(node)) {
                return ghostIds.contains(tree.getTaxonId(node.getNumber()));
            }
            return isGhost(tree.getChild(node, 0)) && isGhost(tree.getChild(node, 1));
        }

        /**
         * A move that touches ghosts alone: two ghost subtrees exchanged, the height of a ghost node, the rate of a
         * ghost branch, or a ghost subtree whose parent is a ghost node moved with it onto another ghost branch.
         *
         * @return false if the chosen kind of move was not possible
         */
        boolean proposeAmongGhosts() {
            final int kind = proposals.nextInt(4);
            final int nodeCount = tree.getNodeCount();
            for (int attempt = 0; attempt < 1000; attempt++) {
                final NodeRef i = tree.getNode(proposals.nextInt(nodeCount));
                if (tree.isRoot(i) || !isGhost(i)) {
                    continue;
                }
                final NodeRef iP = tree.getParent(i);
                if (kind == 0) {
                    // the rate of a ghost branch
                    rates.setRate(i.getNumber(), 0.5 + proposals.nextDouble());
                    return true;
                }
                if (kind == 1) {
                    // the height of a ghost node
                    if (tree.isExternal(i)) {
                        continue;
                    }
                    double low = Math.max(tree.getNodeHeight(tree.getChild(i, 0)),
                            tree.getNodeHeight(tree.getChild(i, 1)));
                    double high = tree.getNodeHeight(iP);
                    tree.setNodeHeight(i, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
                    return true;
                }
                final NodeRef j = tree.getNode(proposals.nextInt(nodeCount));
                if (tree.isRoot(j) || !isGhost(j) || j == i) {
                    continue;
                }
                if (kind == 2) {
                    // two ghost subtrees exchanged
                    final NodeRef jP = tree.getParent(j);
                    if (iP == jP || isAncestor(i, j) || isAncestor(j, i) ||
                            tree.getNodeHeight(j) >= tree.getNodeHeight(iP) ||
                            tree.getNodeHeight(i) >= tree.getNodeHeight(jP)) {
                        continue;
                    }
                    tree.beginTreeEdit();
                    tree.removeChild(iP, i);
                    tree.removeChild(jP, j);
                    tree.addChild(jP, i);
                    tree.addChild(iP, j);
                    tree.endTreeEdit();
                    return true;
                }
                // a ghost subtree and its parent, a ghost node, moved onto another ghost branch
                if (tree.isRoot(iP) || !isGhost(iP) || isAncestor(iP, j)) {
                    continue;
                }
                final NodeRef sibling = tree.getChild(iP, 0) == i ? tree.getChild(iP, 1) : tree.getChild(iP, 0);
                final NodeRef grandparent = tree.getParent(iP);
                final NodeRef jP = (j == sibling) ? grandparent : tree.getParent(j);
                final double low = Math.max(tree.getNodeHeight(j), tree.getNodeHeight(i));
                final double high = tree.getNodeHeight(jP);
                if (low >= high) {
                    continue;
                }
                tree.beginTreeEdit();
                tree.removeChild(iP, sibling);
                tree.removeChild(grandparent, iP);
                tree.addChild(grandparent, sibling);
                tree.removeChild(jP, j);
                tree.addChild(iP, j);
                tree.addChild(jP, iP);
                tree.setNodeHeight(iP, low + (high - low) * (0.01 + 0.98 * proposals.nextDouble()));
                tree.endTreeEdit();
                return true;
            }
            return false;
        }

        private boolean isAncestor(NodeRef ancestor, NodeRef node) {
            for (NodeRef n = node; n != null; n = tree.getParent(n)) {
                if (n == ancestor) {
                    return true;
                }
            }
            return false;
        }
    }

    /**
     * A rate for each branch, indexed by node number, that a proposal can change one at a time.
     */
    private static class NodeRates extends AbstractBranchRateModel {

        private final Parameter rates;

        NodeRates(double[] values) {
            super("nodeRates");
            rates = new Parameter.Default(values);
            addVariable(rates);
        }

        void setRate(int node, double rate) {
            rates.setParameterValue(node, rate);
        }

        @Override
        public double getBranchRate(Tree tree, NodeRef node) {
            return rates.getParameterValue(node.getNumber());
        }

        @Override
        protected void handleModelChangedEvent(Model model, Object object, int index) {
        }

        @Override
        protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
            fireModelChanged(null, index);
        }

        @Override
        protected void storeState() {
        }

        @Override
        protected void restoreState() {
        }

        @Override
        protected void acceptState() {
        }

        @Override
        public Model getModel() {
            return this;
        }

        @Override
        public double getLogLikelihood() {
            return 0.0;
        }

        @Override
        public void makeDirty() {
        }
    }
}
