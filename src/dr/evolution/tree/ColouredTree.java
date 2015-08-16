/*
 * ColouredTree.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree;

import dr.evolution.datatype.DataType;

import java.util.List;

/**
 * @author Alexei Drummond
 *
 * @version $Id: ColouredTree.java,v 1.7 2005/04/28 09:29:07 alexei Exp $
 */
public interface ColouredTree extends Tree {

    /**
     * @param node
     * @return the colour of the given node.
     */
    int getColour(NodeRef node);

    /**
     * This method will return the colour on the ancestral branch if the given time
     * is older than the parent of the given node.
     * @param node
     * @return the colour of the branch above the given node at the given (absolute) time.
     */
    int getColour(NodeRef node, double time);


    /**
     * @return the dataType that describes the colours.
     */
    DataType getColourDataType();

    /**
     * @param node
     * @return an list of colour change objects representing the colour changes on the branch above this node
     */
    List getColourChanges(NodeRef node);

    public class Utils {

        /**
         * @param tree the tree to count the changes on
         * @return the total number of changes on the given tree
         */
        public static int getChangeCount(ColouredTree tree) {

            int changeCount = 0;
            int nodeCount = tree.getNodeCount();
            for (int i = 0; i < nodeCount; i++) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    changeCount += tree.getColourChanges(tree.getNode(i)).size();
                }
            }
            return changeCount;
        }
    }
}
