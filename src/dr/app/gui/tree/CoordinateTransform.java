/*
 * CoordinateTransform.java
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

package dr.app.gui.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 * Date: Dec 3, 2004
 * Time: 9:59:49 AM
 *
 */
public interface CoordinateTransform {

    /**
     * @param height in the time units of the tree.
     * @return the Graphics2D x coordinate corresponding to the given height in the tree.
     */
    double xCoordinate(double height);

    /**
     * @param tree the tree
     * @param node the node whose y coordinate is in question.
     * @return the Graphics2D y coordinate corresponding to the given node.
     */
    double yCoordinate(Tree tree, NodeRef node);
}
