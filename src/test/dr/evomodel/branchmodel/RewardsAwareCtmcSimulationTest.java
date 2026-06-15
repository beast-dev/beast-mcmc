/*
 * RewardsAwareCtmcSimulationTest.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.evomodel.branchmodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.RewardsAwareCtmcSimulation;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchmodel.RewardsAwareCtmcSimulationParser;
import dr.evomodelxml.branchratemodel.RewardsAwareMixtureBranchRatesParser;
import dr.inference.model.Parameter;
import dr.xml.XMLObject;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import test.dr.math.MathTestCase;

import javax.xml.parsers.DocumentBuilderFactory;
import java.lang.reflect.Method;

/*
 * @author Filippo Monti
 */
public class RewardsAwareCtmcSimulationTest extends MathTestCase {

    public void testAbsorbingRootStateProducesAtomicBranchesAndRewards() {
        final TreeModel tree = createThreeTipTree();
        final double[] q = new double[]{0.0, 0.0, 0.0, 0.0};
        final double[] rewards = new double[]{0.25, 1.50};

        final RewardsAwareCtmcSimulation simulation =
                RewardsAwareCtmcSimulation.withFixedRootState(
                        "absorbingRewardSimulation",
                        2,
                        rewards,
                        tree,
                        q,
                        1
                );

        assertEquals(0, simulation.getSummary().totalJumps);
        assertEquals(tree.getNodeCount() - 1, simulation.getSummary().atomicBranchCount);
        assertEquals(0, simulation.getSummary().continuousBranchCount);
        assertEquals(tree.getExternalNodeCount(), simulation.getSummary().tipStateCounts[1]);

        final int[] indicators = simulation.getIndicators();
        final int[] atomIndices = simulation.getAtomIndices();
        final double[] branchRewardTotals = simulation.getBranchRewardTotals();
        final double[] branchRewardProportions = simulation.getBranchRewardProportions();

        for (int branch = 0; branch < indicators.length; branch++) {
            assertEquals(1, indicators[branch]);
            assertEquals(1, atomIndices[branch]);

            final int nodeNumber = simulation.getNodeNumberForBranchIndex(branch);
            final NodeRef node = tree.getNode(nodeNumber);
            assertEquals(tree.getBranchLength(node) * rewards[1], branchRewardTotals[branch], 1.0e-12);
            assertEquals(rewards[1], branchRewardProportions[branch], 1.0e-12);
        }
    }

    public void testParserFillsProvidedRewardAwareParameters() throws Exception {
        final TreeModel tree = createThreeTipTree();
        final int branchCount = tree.getNodeCount() - 1;
        final Parameter branchRewardTotals = new Parameter.Default("simulatedTotals", branchCount, -1.0);
        final Parameter ctsRewards = new Parameter.Default("simulatedCts", branchCount, -1.0);
        final Parameter indicators = new Parameter.Default("simulatedIndicators", branchCount, -1.0);
        final Parameter atomIndices = new Parameter.Default("simulatedAtoms", branchCount, -1.0);

        final Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().newDocument();
        final XMLObject xo = xmlObject(document,
                RewardsAwareCtmcSimulationParser.PARSER_NAME,
                RewardsAwareCtmcSimulationParser.STATE_COUNT, "2",
                RewardsAwareCtmcSimulationParser.ROOT_STATE, "1");

        addChild(xo, nativeObject(document, "treeModel", tree));
        addChild(xo, wrapper(document, RewardsAwareCtmcSimulationParser.Q_MATRIX,
                nativeObject(document, "parameter", new Parameter.Default("q", new double[]{0.0, 0.0, 0.0, 0.0}))));
        addChild(xo, wrapper(document, RewardsAwareCtmcSimulationParser.REWARD_RATES,
                nativeObject(document, "parameter", new Parameter.Default("rewardRates", new double[]{0.25, 1.50}))));
        addChild(xo, wrapper(document, RewardsAwareCtmcSimulationParser.BRANCH_REWARD_TOTALS,
                nativeObject(document, "parameter", branchRewardTotals)));
        addChild(xo, wrapper(document, RewardsAwareMixtureBranchRatesParser.CTS,
                nativeObject(document, "parameter", ctsRewards)));
        addChild(xo, wrapper(document, RewardsAwareMixtureBranchRatesParser.INDICATOR,
                nativeObject(document, "parameter", indicators)));
        addChild(xo, wrapper(document, RewardsAwareMixtureBranchRatesParser.ATOMS_INDICES,
                nativeObject(document, "parameter", atomIndices)));

        final Object parsed = new RewardsAwareCtmcSimulationParser().parseXMLObject(xo);
        assertTrue(parsed instanceof RewardsAwareCtmcSimulation);

        for (int branch = 0; branch < branchCount; branch++) {
            assertTrue(branchRewardTotals.getParameterValue(branch) >= 0.0);
            assertEquals(1.50, ctsRewards.getParameterValue(branch), 1.0e-12);
            assertEquals(1.0, indicators.getParameterValue(branch), 1.0e-12);
            assertEquals(1.0, atomIndices.getParameterValue(branch), 1.0e-12);
        }
    }

    private static TreeModel createThreeTipTree() {
        final SimpleNode left = new SimpleNode();
        left.setTaxon(new Taxon("a"));
        left.setHeight(0.0);

        final SimpleNode right = new SimpleNode();
        right.setTaxon(new Taxon("b"));
        right.setHeight(0.0);

        final SimpleNode internal = new SimpleNode();
        internal.setHeight(1.0);
        internal.addChild(left);
        internal.addChild(right);

        final SimpleNode outgroup = new SimpleNode();
        outgroup.setTaxon(new Taxon("c"));
        outgroup.setHeight(0.0);

        final SimpleNode root = new SimpleNode();
        root.setHeight(2.0);
        root.addChild(internal);
        root.addChild(outgroup);

        final Tree tree = new SimpleTree(root);
        return new DefaultTreeModel("rewardAwareCtmcSimulationTree", tree);
    }

    private static XMLObject xmlObject(final Document document, final String name, final String... attributes) {
        final Element element = document.createElement(name);
        for (int i = 0; i < attributes.length; i += 2) {
            element.setAttribute(attributes[i], attributes[i + 1]);
        }
        return new XMLObject(element, null);
    }

    private static XMLObject nativeObject(final Document document, final String name, final Object object) {
        final XMLObject xo = xmlObject(document, name);
        xo.setNativeObject(object);
        return xo;
    }

    private static XMLObject wrapper(final Document document, final String name, final XMLObject child) throws Exception {
        final XMLObject xo = xmlObject(document, name);
        addChild(xo, child);
        return xo;
    }

    private static void addChild(final XMLObject parent, final XMLObject child) throws Exception {
        final Method addChild = XMLObject.class.getDeclaredMethod("addChild", Object.class);
        addChild.setAccessible(true);
        addChild.invoke(parent, child);
    }
}
