/*
 * SubtreeLeapOperator.java
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
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.NodeHeightOperatorParser;
import dr.inference.operators.AdaptationMode;
import dr.inference.operators.RandomWalkOperator;
import dr.math.MathUtils;
import dr.util.Transform;

/**
 * Implements moves that changes node heights - Scale version
 *
 * Created to allow these operations on TreeModels that can't expose node heights as parameters for efficiency
 * reasons such as BigFastTreeModel.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public class ScaleNodeHeightOperator extends AbstractAdaptableTreeOperator {


    private double tunableParameter;

    private final TreeModel tree;
    private final NodeHeightOperatorParser.OperatorType operatorType;

    /**
     * Constructor
     *
     * @param tree   the tree
     * @param weight the weight
     * @param tunableParameter   size of move for types that need it
     * @param targetAcceptance the desired acceptance probability
     * @param operatorType the type of move to make
     * @param mode   coercion mode
     */
    public ScaleNodeHeightOperator(TreeModel tree, double weight, double tunableParameter, NodeHeightOperatorParser.OperatorType operatorType, AdaptationMode mode, double targetAcceptance) {
        super(mode, targetAcceptance);

        if (operatorType != NodeHeightOperatorParser.OperatorType.SCALEALL && operatorType != NodeHeightOperatorParser.OperatorType.SCALEROOT) {
            throw new UnsupportedOperationException("ScaleNodeHeightOperator can only perform SCALEALL or SCALEROOT operator types");
        }
        this.tree = tree;
        setWeight(weight);
        this.tunableParameter = tunableParameter;
        this.operatorType = operatorType;
    }

    /**
     * Do a subtree leap move.
     *
     * @return the log-transformed hastings ratio
     */
    public double doOperation() {
        final NodeRef root = tree.getRoot();
        double logq;

        if (operatorType == NodeHeightOperatorParser.OperatorType.SCALEALL) {
            logq = doScaleAll();
        } else {
            //if (operatorType == OperatorType.SCALEROOT) {
            final double lowerHeight = Math.max(
                    tree.getNodeHeight(tree.getChild(root, 0)),
                    tree.getNodeHeight(tree.getChild(root, 1)));

            final double oldHeight = tree.getNodeHeight(root);
            logq = doScaleRoot(root, oldHeight, lowerHeight);
        }

        return logq;
    }

    private double doScaleAll() {
        final double scaleFactor = tunableParameter;
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        for (NodeRef node : tree.getNodes()) {
            if (!tree.isExternal(node)) {
                double h = tree.getNodeHeight(node);
                // set quietly so a single update message is sent at the end
                tree.setNodeHeightQuietly(node, h * scale);
            }
        }

        if (!tree.isTreeValid()) {
            // node heights are no long valid (i.e. may be below tips) so force a move reject
            return Double.NEGATIVE_INFINITY;
        }

        tree.pushTreeChangedEvent(TreeChangedEvent.create(false, true));

        return (tree.getInternalNodeCount() - 2) * Math.log(scale);
    }

    private double doScaleRoot(NodeRef node, double oldValue, double lower) {
        final double scaleFactor = tunableParameter;
        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        double h = oldValue - lower;
        tree.setNodeHeight(node, h * scale + lower);
        return -Math.log(scale);
    }

    @Override
    protected void setAdaptableParameterValue(double value) {
            tunableParameter = 1.0 / (Math.exp(value) + 1.0);
    }

    @Override
    protected double getAdaptableParameterValue() {
            return Math.log(1.0 / tunableParameter - 1.0);
    }

    @Override
    public double getRawParameter() {
        return tunableParameter;
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public String getOperatorName() {
        switch (operatorType) {
            case SCALEROOT:
                return "scale(" + tree.getId() + " root height)";
            case SCALEALL:
                return "scaleAll(" + tree.getId() + " internal nodes)";
            default:
                throw new IllegalArgumentException("Unsupported OperatorType");
        }
    }
}
