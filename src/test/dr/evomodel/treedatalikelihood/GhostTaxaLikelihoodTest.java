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
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
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

    private enum Arrangement { SINGLES, CHERRY, CLADE, ROOT_SISTER, ONE_REAL, MIXED }

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

                        for (int iteration = 0; iteration < 40; ++iteration) {
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

        private final Random proposals;

        Fixture(long seed, Arrangement arrangement, boolean leaveOut, PartialsRescalingScheme scheme,
                int categories, boolean useAmbiguities) {
            Random random = new Random(seed);
            proposals = new Random(seed + 1000);

            tree = new DefaultTreeModel(randomTree(random));

            Set<String> ghosts = chooseGhosts(arrangement, random);

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
            HomogeneousBranchModel branchModel = new HomogeneousBranchModel(hky, rootFrequencies);
            GammaSiteRateModel siteRateModel = new GammaSiteRateModel("siteModel", 0.5, categories);

            BeagleDataLikelihoodDelegate delegate = new BeagleDataLikelihoodDelegate(tree, patterns, branchModel,
                    siteRateModel, useAmbiguities, false, scheme, false, PreOrderSettings.getDefault(), null, null,
                    leaveOut);

            likelihood = new TreeDataLikelihood(delegate, tree, new DefaultBranchRateModel());
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
            if (proposals.nextInt(3) == 0) {
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

        private boolean isAncestor(NodeRef ancestor, NodeRef node) {
            for (NodeRef n = node; n != null; n = tree.getParent(n)) {
                if (n == ancestor) {
                    return true;
                }
            }
            return false;
        }
    }
}
