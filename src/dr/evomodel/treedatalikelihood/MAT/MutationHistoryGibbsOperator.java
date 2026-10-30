/*
 * MutationHistoryGibbsOperator.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

import java.util.Arrays;

/**
 * A Gibbs update of a MAT's mutational history, without changing the tree: the states of a random internal
 * node and of up to size - 1 of its ancestors are drawn from their full conditional given all the other
 * states (see {@link MutationHistoryResampler#resampleInPlace}).
 * <p>
 * The subtree moves resample the history only around the nodes they move, so a node's states otherwise
 * change only when it is moved. This operator lets every part of the history mix, including moving a
 * mutation across several branches at once. The local block needs it: without it, excess mutations in the
 * history are almost never removed, and the tree parameters follow them. A weight of about 10% of the
 * subtree moves' weight is enough.
 */
public class MutationHistoryGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String OPERATOR_NAME = "mutationHistoryGibbs";

    private final MutationHistoryResampler resampler;
    private final Tree tree;
    private final int size;

    /**
     * @param size the largest number of nodes in a block (a node and its ancestors)
     */
    public MutationHistoryGibbsOperator(MutationHistoryResampler resampler, int size, double weight) {
        if (size < 1) {
            throw new IllegalArgumentException("The block size must be at least 1");
        }
        this.resampler = resampler;
        resampler.setGibbsOperator();
        this.tree = resampler.getDelegate().getTree();
        this.size = size;
        setWeight(weight);
    }

    @Override
    public double doOperation() {
        resampler.resampleInPlace(chooseBlock());
        return 0.0;
    }

    /**
     * @return a random internal node and up to size - 1 of its ancestors
     */
    int[] chooseBlock() {
        NodeRef node = tree.getInternalNode(MathUtils.nextInt(tree.getInternalNodeCount()));
        int[] block = new int[size];
        int count = 0;
        for (NodeRef n = node; n != null && count < size; n = tree.getParent(n)) {
            block[count++] = n.getNumber();
        }
        return count == size ? block : Arrays.copyOf(block, count);
    }

    @Override
    protected String getBaseOperatorName() {
        return OPERATOR_NAME + "(" + tree.getId() + ")";
    }
}
