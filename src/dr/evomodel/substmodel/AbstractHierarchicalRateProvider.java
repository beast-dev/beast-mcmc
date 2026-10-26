/*
 * AbstractHierarchicalRateProvider.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Base class of the rate providers that compute the per-clade rates (lambdas) of a
 * HierarchicalSubstitutionModel from familiar parameters. The lambda vector is cached and recomputed
 * once after any change of a parameter, the frequencies or a sub-model.
 *
 * @author Xiang Ji
 */
public abstract class AbstractHierarchicalRateProvider extends AbstractModel
        implements HierarchicalSubstitutionModel.HierarchicalRateProvider {

    public static final String NUCLEOTIDE_NEWICK = "((A:1,G:1):1,(C:1,T:1):1);";

    protected final TreeModel tree;
    protected final FrequencyModel frequency;
    protected final DataType dataType;
    protected final int stateCount;
    protected final int[] stateOfNode;
    protected final int[] nodeOfState;
    protected final int[] parentOf;
    protected final int[] internalPreOrder;
    protected final int[][] statesBelow;
    protected final String[] cladeName;
    protected final double[] piNode;

    private final double[] lambda;
    private boolean known = false;
    private double normalization = Double.NaN;
    private int computeCount = 0;

    protected AbstractHierarchicalRateProvider(String name, TreeModel tree, FrequencyModel frequency) {
        super(name);
        this.tree = tree;
        this.frequency = frequency;
        this.dataType = frequency.getDataType();
        this.stateCount = tree.getExternalNodeCount();
        if (dataType.getStateCount() != stateCount) {
            throw new IllegalArgumentException("Character tree has " + stateCount + " tips but the data type has "
                    + dataType.getStateCount() + " states");
        }
        final int nodeCount = tree.getNodeCount();
        this.stateOfNode = new int[stateCount];
        this.nodeOfState = new int[stateCount];
        Arrays.fill(nodeOfState, -1);
        for (int i = 0; i < stateCount; i++) {
            final int state = dataType.getState(tree.getTaxonId(i));
            if (state < 0 || state >= stateCount || nodeOfState[state] != -1) {
                throw new IllegalArgumentException("Character tree tip " + tree.getTaxonId(i)
                        + " is not a distinct state of " + dataType.getDescription());
            }
            stateOfNode[i] = state;
            nodeOfState[state] = i;
        }

        this.parentOf = new int[nodeCount];
        Arrays.fill(parentOf, -1);
        this.statesBelow = new int[nodeCount][];
        this.cladeName = new String[nodeCount];
        final List<Integer> post = new ArrayList<Integer>();
        postOrder(tree, tree.getRoot(), post);
        for (int n : post) {
            final NodeRef node = tree.getNode(n);
            if (tree.isExternal(node)) {
                statesBelow[n] = new int[]{stateOfNode[n]};
            } else {
                final List<Integer> below = new ArrayList<Integer>();
                for (int c = 0; c < tree.getChildCount(node); c++) {
                    final int child = tree.getChild(node, c).getNumber();
                    parentOf[child] = n;
                    for (int state : statesBelow[child]) {
                        below.add(state);
                    }
                }
                statesBelow[n] = toArray(below);
            }
            cladeName[n] = sortedCodes(statesBelow[n]);
        }
        final List<Integer> pre = new ArrayList<Integer>();
        for (int i = post.size() - 1; i >= 0; i--) {
            if (!tree.isExternal(tree.getNode(post.get(i)))) {
                pre.add(post.get(i));
            }
        }
        this.internalPreOrder = toArray(pre);
        this.piNode = new double[nodeCount];
        this.lambda = new double[nodeCount];

        addModel(frequency);
    }

    public static TreeModel characterTree(String newick) {
        final Tree imported;
        try {
            imported = new NewickImporter(newick).importTree(null);
        } catch (Exception e) {
            throw new IllegalArgumentException("cannot read the character tree " + newick + ": " + e.getMessage());
        }
        return new DefaultTreeModel("characterTree", imported, true, true);
    }

    /**
     * Fill lambda by node number from the current parameters (piNode is up to date); return the normalisation
     * constant.
     */
    protected abstract double computeLambdas(double[] lambdaByNode);

    private void ensure() {
        if (!known) {
            fillPiNode();
            Arrays.fill(lambda, 0.0);
            normalization = computeLambdas(lambda);
            computeCount++;
            known = true;
        }
    }

    private void fillPiNode() {
        for (int n = 0; n < piNode.length; n++) {
            double sum = 0.0;
            for (int state : statesBelow[n]) {
                sum += frequency.getFrequency(state);
            }
            piNode[n] = sum;
        }
    }

    @Override
    public final double getRate(TreeModel tree, int node) {
        ensure();
        return lambda[node];
    }

    @Override
    public final TreeModel getTreeModel() {
        return tree;
    }

    @Override
    public final FrequencyModel getFrequencyModel() {
        return frequency;
    }

    @Override
    public final void addToModelListener(AbstractModel listener) {
        listener.addModel(this);
    }

    public final double getNormalization() {
        ensure();
        return normalization;
    }

    public final double getLambda(String sortedTipCodes) {
        ensure();
        return lambda[nodeOfClade(sortedTipCodes)];
    }

    public final String[] getCladeNames() {
        final String[] names = new String[internalPreOrder.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = cladeName[internalPreOrder[i]];
        }
        return names;
    }

    public final int getComputeCount() {
        return computeCount;
    }

    protected final int nodeOfClade(String sortedTipCodes) {
        for (int n : internalPreOrder) {
            if (cladeName[n].equals(sortedTipCodes)) {
                return n;
            }
        }
        throw new IllegalArgumentException("clade {" + braces(sortedTipCodes) + "} is not a clade of the character tree; its clades are "
                + describeClades());
    }

    protected final String describeClades() {
        final StringBuilder sb = new StringBuilder();
        for (int n : internalPreOrder) {
            if (sb.length() > 0) {
                sb.append(", ");
            }
            sb.append('{').append(braces(cladeName[n])).append('}');
        }
        return sb.toString();
    }

    private static String braces(String codes) {
        final StringBuilder sb = new StringBuilder();
        for (int i = 0; i < codes.length(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            sb.append(codes.charAt(i));
        }
        return sb.toString();
    }

    private String sortedCodes(int[] states) {
        final String[] codes = new String[states.length];
        for (int i = 0; i < states.length; i++) {
            codes[i] = dataType.getCode(states[i]);
        }
        Arrays.sort(codes);
        final StringBuilder sb = new StringBuilder();
        for (String code : codes) {
            sb.append(code);
        }
        return sb.toString();
    }

    private static void postOrder(Tree tree, NodeRef node, List<Integer> order) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            postOrder(tree, tree.getChild(node, i), order);
        }
        order.add(node.getNumber());
    }

    protected static int[] toArray(List<Integer> list) {
        final int[] out = new int[list.size()];
        for (int i = 0; i < out.length; i++) {
            out[i] = list.get(i);
        }
        return out;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        known = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        known = false;
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
        known = false;
    }

    @Override
    protected void acceptState() {
    }
}
