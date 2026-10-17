/*
 * FastTreeModelTest.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package test.dr.evomodel.tree;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.operators.ExchangeOperator;
import dr.evomodel.operators.FixedHeightSubtreePruneRegraftOperator;
import dr.evomodel.operators.SubtreeLeapOperator;
import dr.evomodel.operators.WilsonBalding;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.FastTreeModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.inference.operators.*;
import dr.math.MathUtils;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;

/**
 * Runs the same moves on a DefaultTreeModel and a FastTreeModel built from the same tree, in lock
 * step, and checks that they stay identical: the tree, the values and bounds seen through every
 * height parameter, the events fired, and the trees after store/accept/restore.
 */
public class FastTreeModelTest extends TestCase {

    private static final String NEWICK = "((((A:1.0,B:1.2):1.1,(C:1.5,D:1.4):0.5):1.0,(E:2.0,F:2.3):1.0):1.0," +
            "(((G:1.0,H:1.0):1.0,I:2.5):1.0,(J:0.5,K:0.8):2.5):1.5);";

    public FastTreeModelTest(String name) {
        super(name);
    }

    /**
     * One tree model with its height parameters, operators and a record of the events it fires.
     */
    private static class Setup {
        final TreeModel tree;
        final Parameter rootHeight;
        final Parameter internalHeights;
        final Parameter allInternalHeights;
        final Parameter leafHeightA;
        final Parameter leafHeightJ;
        final SimpleMCMCOperator[] operators;
        final List<String> events = new ArrayList<String>();

        Setup(boolean fast) throws Exception {
            Tree start = new NewickImporter(NEWICK).importTree(null);
            if (fast) {
                FastTreeModel model = new FastTreeModel("treeModel", start, true, false);
                tree = model;
                rootHeight = model.getRootHeightParameter();
                internalHeights = model.createNodeHeightsParameter(false, true, false);
                allInternalHeights = model.createNodeHeightsParameter(true, true, false);
                leafHeightA = model.getLeafHeightParameter(model.getExternalNode(model.getTaxonIndex("A")));
                leafHeightJ = model.getLeafHeightParameter(model.getExternalNode(model.getTaxonIndex("J")));
            } else {
                DefaultTreeModel model = new DefaultTreeModel("treeModel", start, true, false);
                tree = model;
                rootHeight = model.getRootHeightParameter();
                internalHeights = model.createNodeHeightsParameter(false, true, false);
                allInternalHeights = model.createNodeHeightsParameter(true, true, false);
                leafHeightA = model.getLeafHeightParameter(model.getExternalNode(model.getTaxonIndex("A")));
                leafHeightJ = model.getLeafHeightParameter(model.getExternalNode(model.getTaxonIndex("J")));
            }
            // an extra bound on a tip height (as the parser adds for an uncertain tip date)
            leafHeightJ.addBounds(new Parameter.DefaultBounds(0.6, 0.0, 1));

            Parameter rate = new Parameter.Default(1.0);
            rate.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

            operators = new SimpleMCMCOperator[]{
                    new SubtreeLeapOperator(tree, 1.0, 3.0, SubtreeLeapOperator.DistanceKernelType.NORMAL,
                            AdaptationMode.ADAPTATION_OFF, 0.234),
                    new FixedHeightSubtreePruneRegraftOperator(tree, 1.0),
                    new ExchangeOperator(ExchangeOperator.NARROW, tree, 1.0),
                    new ExchangeOperator(ExchangeOperator.WIDE, tree, 1.0),
                    new WilsonBalding(tree, 1.0),
                    new UniformOperator(internalHeights, 1.0),
                    new ScaleOperator(rootHeight, 0.75, AdaptationMode.ADAPTATION_OFF, 1.0),
                    new UpDownOperator(new Scalable[]{new Scalable.Default(allInternalHeights)},
                            new Scalable[]{new Scalable.Default(rate)}, 0.75, 1.0, AdaptationMode.ADAPTATION_OFF),
                    new RandomWalkOperator(leafHeightA, 0.2, RandomWalkOperator.BoundaryCondition.reflecting,
                            1.0, AdaptationMode.ADAPTATION_OFF),
                    new UniformOperator(leafHeightJ, 1.0)
            };

            // record the events, as a tree likelihood or prior would see them
            tree.addModelListener(new ModelListener() {
                public void modelChangedEvent(Model model, Object object, int index) {
                    if (object instanceof TreeChangedEvent) {
                        TreeChangedEvent event = (TreeChangedEvent) object;
                        events.add("tree:" + (event.getNode() == null ? "-" : event.getNode().getNumber()) + ":" +
                                event.isNodeChanged() + event.isTreeChanged() + event.isHeightChanged() +
                                event.isNodeParameterChanged() + event.isOnlyHeightChanged() + ":" + event.getIndex() +
                                ":" + (event.getParameter() == null ? "-" : event.getParameter().getId()));
                    } else if (object instanceof Parameter) {
                        events.add("parameter:" + ((Parameter) object).getId() + ":" + index);
                    } else {
                        events.add("other:" + object + ":" + index);
                    }
                }

                public void modelRestored(Model model) {
                }
            });
            addRecorder(rootHeight, "rootHeight");
            addRecorder(internalHeights, "internalHeights");
            addRecorder(allInternalHeights, "allInternalHeights");
        }

        private void addRecorder(Parameter parameter, final String name) {
            parameter.addParameterListener(new VariableListener() {
                public void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
                    events.add(name + ":" + index + ":" + type);
                }
            });
        }

        void storeAll() {
            tree.storeModelState();
            for (Parameter p : new Parameter[]{rootHeight, internalHeights, allInternalHeights, leafHeightA, leafHeightJ}) {
                p.storeParameterValues();
            }
        }

        void restoreAll() {
            tree.restoreModelState();
            for (Parameter p : new Parameter[]{rootHeight, internalHeights, allInternalHeights, leafHeightA, leafHeightJ}) {
                p.restoreParameterValues();
            }
        }

        void acceptAll() {
            tree.acceptModelState();
            for (Parameter p : new Parameter[]{rootHeight, internalHeights, allInternalHeights, leafHeightA, leafHeightJ}) {
                p.acceptParameterValues();
            }
        }
    }

    public void testLockStepWithDefaultTreeModel() throws Exception {
        Setup expected = new Setup(false);
        Setup actual = new Setup(true);
        assertSame(expected, actual, "start");

        int rootChanges = 0;
        for (int iteration = 0; iteration < 20000; iteration++) {
            expected.storeAll();
            actual.storeAll();
            NodeRef root = actual.tree.getRoot();

            // usually one move, sometimes several, each with the same random numbers for both trees
            int moves = MathUtils.nextInt(4) == 0 ? 2 + MathUtils.nextInt(3) : 1;
            for (int m = 0; m < moves; m++) {
                long seed = 1000L * iteration + m + 1;
                MathUtils.setSeed(seed);
                int op = MathUtils.nextInt(expected.operators.length);
                String e = operate(expected.operators[op]);
                MathUtils.setSeed(seed);
                MathUtils.nextInt(actual.operators.length);
                String a = operate(actual.operators[op]);
                assertEquals("operator " + op + " at " + iteration, e, a);
                assertSame(expected, actual, "after operator " + op + " at " + iteration);
            }
            if (actual.tree.getRoot() != root) {
                rootChanges++;
            }

            // accept half the time, but (in place of a prior) not if the tree has grown too tall for
            // the moves to change the root
            MathUtils.setSeed(iteration + 1);
            if (MathUtils.nextBoolean() && expected.tree.getNodeHeight(expected.tree.getRoot()) < 8.0) {
                expected.acceptAll();
                actual.acceptAll();
            } else {
                expected.restoreAll();
                actual.restoreAll();
                assertSame(expected, actual, "after restore at " + iteration);
            }
        }
        assertEquals(expected.events, actual.events);
        // make sure the moves did change the tree, including the root (and so swap the root height slot)
        assertTrue("root changes " + rootChanges, rootChanges > 100);
        assertTrue("events " + actual.events.size(), actual.events.size() > 100000);
    }

    /**
     * Makes a move and returns its outcome (the Hastings ratio, or the exception or assertion failure).
     */
    private static String operate(SimpleMCMCOperator operator) {
        try {
            return Double.toString(operator.doOperation());
        } catch (RuntimeException e) {
            return e.getClass().getName();
        } catch (AssertionError e) {
            // e.g. the tree validity check at the end of an edit (with assertions enabled)
            return e.getClass().getName();
        }
    }

    private static void assertSame(Setup expected, Setup actual, String message) {
        Tree e = expected.tree;
        Tree a = actual.tree;
        assertEquals(message, e.getRoot().getNumber(), a.getRoot().getNumber());
        for (int i = 0; i < e.getNodeCount(); i++) {
            NodeRef ne = e.getNode(i);
            NodeRef na = a.getNode(i);
            String m = message + ", node " + i;
            assertEquals(m, e.getNodeHeight(ne), a.getNodeHeight(na), 0.0);
            assertEquals(m, number(e.getParent(ne)), number(a.getParent(na)));
            assertEquals(m, e.getChildCount(ne), a.getChildCount(na));
            assertEquals(m, number(e.getChild(ne, 0)), number(a.getChild(na, 0)));
            assertEquals(m, number(e.getChild(ne, 1)), number(a.getChild(na, 1)));
            assertEquals(m, e.isRoot(ne), a.isRoot(na));
            assertEquals(m, e.isExternal(ne), a.isExternal(na));
            assertEquals(m, ((TreeModel) e).getNodeTaxon(ne), ((TreeModel) a).getNodeTaxon(na));
        }
        assertParameter(expected.rootHeight, actual.rootHeight, message + ", root height");
        assertParameter(expected.internalHeights, actual.internalHeights, message + ", internal heights");
        assertParameter(expected.allInternalHeights, actual.allInternalHeights, message + ", all internal heights");
        assertParameter(expected.leafHeightA, actual.leafHeightA, message + ", leaf height A");
        assertParameter(expected.leafHeightJ, actual.leafHeightJ, message + ", leaf height J");
        assertEquals(message, expected.events.size(), actual.events.size());
    }

    private static void assertParameter(Parameter e, Parameter a, String message) {
        assertEquals(message, e.getDimension(), a.getDimension());
        assertEquals(message, e.getParameterName(), a.getParameterName());
        for (int i = 0; i < e.getDimension(); i++) {
            assertEquals(message + " [" + i + "]", e.getParameterValue(i), a.getParameterValue(i), 0.0);
            assertEquals(message + " lower [" + i + "]", e.getBounds().getLowerLimit(i), a.getBounds().getLowerLimit(i), 0.0);
            assertEquals(message + " upper [" + i + "]", e.getBounds().getUpperLimit(i), a.getBounds().getUpperLimit(i), 0.0);
        }
    }

    private static int number(NodeRef node) {
        return node == null ? -1 : node.getNumber();
    }

    /**
     * Node rates and node traits are not supported.
     */
    public void testNodeRatesAndTraitsUnsupported() throws Exception {
        FastTreeModel tree = new FastTreeModel("treeModel", new NewickImporter(NEWICK).importTree(null), true, false);
        NodeRef node = tree.getNode(0);
        try {
            tree.createNodeRatesParameter(null, true, true, true);
            fail();
        } catch (UnsupportedOperationException e) {
            // expected
        }
        try {
            tree.setNodeRate(node, 2.0);
            fail();
        } catch (UnsupportedOperationException e) {
            // expected
        }
        try {
            tree.getNodeTrait(node, "trait");
            fail();
        } catch (UnsupportedOperationException e) {
            // expected
        }
        // as DefaultTreeModel without node rates: generic tree copying asks for the rate
        assertEquals(1.0, tree.getNodeRate(node), 0.0);
    }
}
