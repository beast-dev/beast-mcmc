/*
 * CorporealTreeModel.java
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

package dr.evomodel.bigfasttree.ghosttree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.bigfasttree.BigFastTreeModel;
import dr.evomodel.bigfasttree.ghosttree.GhostTreeModel;

/**
 * A tree model that does not allow for the usual edits so that it's unlikely to be
 * edited outside of the ghost tree it shadows.
 *
 * @author JT McCrone
 */
public class CorporealTreeModel extends BigFastTreeModel {
    public static final String CORPOREAL_TREE_MODEL = "corporealTreeModel";

    private final GhostTreeModel ghostTreeModel;
    public CorporealTreeModel(String name, Tree tree,GhostTreeModel ghostTreeModel) {
        super(name, tree);
        this.ghostTreeModel = ghostTreeModel;
        addModel(ghostTreeModel);
    }
    public CorporealTreeModel(Tree tree, GhostTreeModel ghostTreeModel) {
        this(CORPOREAL_TREE_MODEL, tree,ghostTreeModel);
    }

    public GhostTreeModel getGhostTreeModel(){
        return this.ghostTreeModel;
    }
    // *****************************************************************
    // Interface MutableTree
    // *****************************************************************

    //The tree should only be edited by the ghosttree model. The usual edit methods throw errors
    /**
     * Set a new node as root node.
     */
    public void setRoot(NodeRef newRoot) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void makeRoot(NodeRef newRoot) {
        super.setRoot(newRoot);
    }

    public void addChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void adoptChild(NodeRef p, NodeRef c) {
        super.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void disownChild(NodeRef p, NodeRef c) {
        super.removeChild(p, c);
    }

    public boolean beginTreeEdit() {
        return super.beginTreeEdit();
    }

    public void endTreeEdit() {
        // and cleanup
        super.endTreeEdit();
    }

    public void setNodeHeight(NodeRef n, double height) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setNodeHeightQuietly(NodeRef n, double height) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }
    protected void adjustNodeHeight(NodeRef n, double height) {
        super.setNodeHeight(n, height);
    }

    protected void adjustNodeHeightQuietly(NodeRef n, double height) {
        super.setNodeHeightQuietly(n, height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setNodeTrait(NodeRef n, String name, double value) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new UnsupportedOperationException("Function not available in CorporealTreeModel");
    }


}
