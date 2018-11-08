/*
 * RateScaleOperator.java
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

package dr.evomodel.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.RateScaleOperatorParser;
import dr.inference.operators.AbstractAdaptableOperator;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.OperatorUtils;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A special operator for scaling rates in a subtree.
 *
 * @author Michael Defoin Platel
 */
public class RateScaleOperator extends AbstractAdaptableOperator {

    private TreeModel tree;

    private boolean noRoot;

    public RateScaleOperator(TreeModel tree, double scale, boolean noRoot, AdaptationMode mode) {

        super(mode);

        this.scaleFactor = scale;
        this.tree = tree;

        this.noRoot = noRoot;
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        int index;
        if (noRoot) {
            do {
                index = MathUtils.nextInt(tree.getNodeCount());
            } while (index == tree.getRoot().getNumber());
        } else {
            index = MathUtils.nextInt(tree.getNodeCount());
        }

        List<NodeRef> listNode = new ArrayList<NodeRef>();
        getSubtree(listNode, tree.getNode(index));

        double oldValue, newValue;
        double logq = 0;
        for (NodeRef node : listNode) {

            oldValue = tree.getNodeRate(node);
            newValue = oldValue * scale;

            tree.setNodeRate(node, newValue);
        }

        //  According to the hastings ratio defined in the original scale Operator
        logq = (listNode.size() - 2) * Math.log(scale);

        return logq;
    }

    void getSubtree(List<NodeRef> listNode, NodeRef parent) {

        listNode.add(parent);
        int nbChildren = tree.getChildCount(parent);
        for (int c = 0; c < nbChildren; c++) {
            getSubtree(listNode, tree.getChild(parent, c));
        }
    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "rateScale";
    }

    public double getAdaptableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setAdaptableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public String toString() {
        return RateScaleOperatorParser.SCALE_OPERATOR + "(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}
