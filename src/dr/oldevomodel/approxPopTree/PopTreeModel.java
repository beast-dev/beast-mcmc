/*
 * PopTreeModel.java
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

package dr.oldevomodel.approxPopTree;

import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.Map;

/**
 * Package: PopTree
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: Feb 1, 2010
 *         Time: 11:39:59 AM
 */
public class PopTreeModel extends FlexibleTree {
    protected Map<NodeRef, LinkedList<NodeRef>> populations;
    protected double time;

    public PopTreeModel(TreeModel treeModel, double time) throws InvalidTreeException {
        super(treeModel);
        this.time = time;
        this.populations = new HashMap<NodeRef, LinkedList<NodeRef>>();
        enforcePopulations(time);
    }

    public void enforcePopulations(double time) {
        enforcePopulations(time, this.getRoot());
    }

    public void enforcePopulations(double time, NodeRef node) {
        beginTreeEdit();
        reducePopulations(time, node);
        endTreeEdit();
        this.adoptTreeModelOrdering();
    }

    protected void reducePopulations(double time, NodeRef node) {
        if (this.getNodeHeight(node) <= time) {
            if (!this.isExternal(node))
                this.setNodeHeight(node, time);
            mergeSubtreePopulation(node);

        } else {
            for (int i = 0; i < this.getChildCount(node); ++i) {
                NodeRef child = this.getChild(node, i);
                reducePopulations(time, child);
            }
        }
    }

    public void splitPopulation(NodeRef node) {
        splitPopulation(node, .5);
    }

    public void splitPopulation(NodeRef node, double frac) {
        if (populations.containsKey(node) && populations.get(node).size() >= 1) {
            //
        }
    }

    public void mergePopulations(NodeRef parent) {
        // assumes parent is an internal node which will become a new population node with all the crap below merged into this population
        setNodeHeight(parent, time);
        enforcePopulations(time, parent);
    }

    protected void mergeSubtreePopulation(NodeRef parent) {
        /*
        merge the subtree below node into a new population
         */
        LinkedList<NodeRef> sequenceNodes;
        sequenceNodes = new LinkedList<NodeRef>();
        if (this.getNodeTaxon(parent) == null) {
            this.setNodeTaxon(parent, new Taxon("popNode" + populations.size()));
        }

        if (this.isExternal(parent)) {

        }

        sequenceNodes.addFirst(parent); //push(parent);
        while (true) {
            NodeRef currentNode = sequenceNodes.removeFirst();
            if (this.isExternal(currentNode)) {
                sequenceNodes.addFirst(currentNode); //push(currentNode);
                //this.setNodeAttribute(parent,"sequenceNodes",sequenceNodes);
                populations.put(parent, sequenceNodes);
                return;
            }
            int numChildren = this.getChildCount(currentNode);
            for (int i = 0; getChildCount(currentNode) > 0;) {
                NodeRef child = this.getChild(currentNode, i);
                sequenceNodes.addFirst(child); //push(child);
                this.removeChild(currentNode, child);
            }
        }
    }

    public int getPopulationNodeCount() {
        return populations.size();
    }

    public String toString() {
        return super.toString() + formatPopulationNodes();
    }

    public String formatPopulationNodes() {
        String retval = "\n";
        retval += getExternalNodeCount() + "\n";
        for (NodeRef node : populations.keySet()) {
            retval += this.getNodeTaxon(node).getId() + "(" + this.isExternal(node) + ")::";
            for (NodeRef seqNode : populations.get(node)) {
                retval += this.getNodeTaxon(seqNode).getId() + ",";
            }
            retval += "\n";
        }
        return retval;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            double timeCutoff = xo.getDoubleAttribute(POP_HEIGHT_CUTOFF);
            PopTreeModel popTreeModel = null;
            try {
                popTreeModel = new PopTreeModel(treeModel, timeCutoff);
            } catch (InvalidTreeException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
            return popTreeModel;
        }

        public XMLSyntaxRule[] getSyntaxRules() {

            return rules;  //AUTOGENERATED METHOD IMPLEMENTATION
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeModel.class),
                AttributeRule.newDoubleRule(POP_HEIGHT_CUTOFF, false),
        };

        public String getParserDescription() {
            return "Creates a Population Tree Model with specified divergence cut-off for population nodes";
        }

        public Class getReturnType() {
            return PopTreeModel.class;
        }

        public String getParserName() {
            return POP_TREE_MODEL;
        }

    };

    public static final String POP_TREE_MODEL = "popTreeModel";
    public static final String POP_HEIGHT_CUTOFF = "populationNodeHeight";
    public static final String TREE_MODEL = "treeModel";
}
