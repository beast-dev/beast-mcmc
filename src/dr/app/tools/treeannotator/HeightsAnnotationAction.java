/*
 * HeightsAnnotationAction.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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

package dr.app.tools.treeannotator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.geo.contouring.ContourMaker;
import dr.geo.contouring.ContourPath;
import dr.geo.contouring.ContourWithSynder;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.*;

public class HeightsAnnotationAction implements CladeAction {
    private final TreeAnnotator.HeightsSummary heightsOption;
    private final double[] hpdIntervals;
    private final int hpdLimit;
    private final double[] kdeIntervals;
    private final int kdeCount;
    private final int kdeLimit;
    private final double posteriorLimit;
    private final double countLimit;
    private final static boolean PROCESS_BIVARIATE_ATTRIBUTES = true;

    HeightsAnnotationAction(final TreeAnnotator.HeightsSummary heightsOption,
                            final double[] hpdIntervals,
                            final int hpdLimit,
                            final double[] kdeIntervals,
                            final int kdeCount,
                            final int kdeLimit,
                            final double posteriorLimit,
                            final int countLimit) {
        this.heightsOption = heightsOption;
        this.posteriorLimit = posteriorLimit;
        this.countLimit = countLimit;
        this.hpdIntervals = hpdIntervals;
        this.hpdLimit = hpdLimit;
        this.kdeIntervals = kdeIntervals;
        this.kdeCount = kdeCount;
        this.kdeLimit = kdeLimit;
    }

    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        assert tree instanceof MutableTree;
        annotateNode((MutableTree)tree, node, clade);
    }

    @Override
    public boolean expectAllClades() {
        return false;
    }

    private void annotateNode(MutableTree tree, NodeRef node, Clade clade) {
        boolean filter = false;
        assert clade != null;

        if (!tree.isExternal(node)) {
            final double posterior = clade.getCredibility();
            tree.setNodeAttribute(node, "posterior", posterior);
            if (posterior < posteriorLimit || clade.getCount() < countLimit) {
                filter = true;
            }
        }

        if (!filter) {
            tree.setNodeAttribute(node, "height_mean", ((BiClade) clade).getMeanHeight());
            tree.setNodeAttribute(node, "height_median", ((BiClade) clade).getMedianHeight());
            if (((BiClade) clade).getHeightHPDs() != null){
                tree.setNodeAttribute(node, "height_95%_HPD", ((BiClade) clade).getHeightHPDs());
            }
            if (((BiClade) clade).getHeightRange() != null){
                tree.setNodeAttribute(node, "height_range", ((BiClade) clade).getHeightRange());
            }
        }
        if (heightsOption == TreeAnnotator.HeightsSummary.MEAN_HEIGHTS) {
            tree.setNodeHeight(node, ((BiClade) clade).getMeanHeight());
        } else if (heightsOption == TreeAnnotator.HeightsSummary.MEDIAN_HEIGHTS) {
            tree.setNodeHeight(node, ((BiClade) clade).getMedianHeight());
        } else {
            // keep the existing height
        }
//        assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 0))) >= 0.0;
//        assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 1))) >= 0.0;
    }

}
