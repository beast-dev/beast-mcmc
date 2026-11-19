/*
 * GMRFUtils.java
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

package dr.evomodel.coalescent.gmrf;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.math.MathUtils;
import no.uib.cipr.matrix.UpperTriangBandMatrix;

/**
 * Functions shared by the Gaussian Markov random field (and Gaussian process) coalescent models.
 */
public class GMRFUtils {

    // the truncated value that these models have always used
    public static final double LOG_TWO_TIMES_PI = 1.837877;

    private GMRFUtils() {
    }

    /**
     * @return the log determinant of a Cholesky factor, leaving out any (near) zero elements of the diagonal
     */
    public static double logGeneralizedDeterminant(UpperTriangBandMatrix Matrix) {
        double returnValue = 0;

        for (int i = 0; i < Matrix.numColumns(); i++) {
            if (Matrix.get(i, i) > 0.0000001) {
                returnValue += Math.log(Matrix.get(i, i));
            }
        }

        return returnValue;
    }

    /**
     * Draws a new height for each internal node, between its parent and its older child.
     */
    public static void checkTree(TreeModel treeModel) {

        // todo Should only be run if there exists a zero-length interval

//        TreeModel treeModel = (TreeModel) tree;
        for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {
            NodeRef node = treeModel.getInternalNode(i);
            if (node != treeModel.getRoot()) {
                double parentHeight = treeModel.getNodeHeight(treeModel.getParent(node));
                double childHeight0 = treeModel.getNodeHeight(treeModel.getChild(node, 0));
                double childHeight1 = treeModel.getNodeHeight(treeModel.getChild(node, 1));
                double maxChild = childHeight0;
                if (childHeight1 > maxChild)
                    maxChild = childHeight1;
                double newHeight = maxChild + MathUtils.nextDouble() * (parentHeight - maxChild);
                treeModel.setNodeHeight(node, newHeight);
            }
        }
        treeModel.pushTreeChangedEvent();

    }
}
