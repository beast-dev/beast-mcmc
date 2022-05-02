/*
 * NodeHeightScaleOperator.java
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
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inferencexml.operators.ScaleOperatorParser;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.*;

/**
 * An operator for scaling 
 *
 * @author Andrew Rambaut
 */
// Cleaning out untouched stuff. Can be resurrected if needed
@Deprecated
public class NodeHeightScaleOperator extends AbstractAdaptableOperator {

    public static final String NODE_HEIGHT_SCALE_OPERATOR = "nodeHeightScaleOperator";
    public static final String SCALE_FACTOR = ScaleOperatorParser.SCALE_FACTOR;
    public static final String SCALE_ALL = ScaleOperatorParser.SCALE_ALL;

    private final TreeModel tree;
    private final boolean scaleAll;
    Set<Double> tipDates = new TreeSet<Double>();
    private final Map<Double, Integer> intervals = new TreeMap<Double, Integer>();
    private Parameter nodeHeightParameter;

    public NodeHeightScaleOperator(TreeModel tree, double scale, boolean scaleAll, AdaptationMode mode) {

        super(mode);

        this.scaleFactor = scale;
        this.scaleAll = scaleAll;
        this.tree = tree;
//        nodeHeightParameter = tree.get

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            double h = tree.getNodeHeight(tree.getExternalNode(i));
            tipDates.add(h);
        }
    }


    /**
     * scale the rates of a subtree and return the hastings ratio.
     */
    public final double doOperation() {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        List<NodeRef> listNode = new ArrayList<NodeRef>();

        double logq = 0;
        if (scaleAll) {

        } else {
            int r = MathUtils.nextInt(tree.getInternalNodeCount());
            NodeRef node = tree.getNode(r);

            double oldValue, newValue;

            oldValue = tree.getNodeHeight(node);
            newValue = oldValue * scale;

            tree.setNodeHeight(node, newValue);

            logq = getJacobian(node);
        }

        //  According to the hastings ratio defined in the original scale Operator
        logq = (listNode.size() - 2) * Math.log(scale);

        return logq;
    }

    private double getJacobian(NodeRef node) {

        intervals.clear();
        for (Double date : tipDates) {
            intervals.put(date, 0);
        }

        traverse(tree, node);

        double logq = 0.0;
        for (Double date : intervals.keySet()) {
            double s = tree.getNodeHeight(tree.getRoot()) - date;
            int k = intervals.get(date);

            logq += logFactorial(k) - (double) k * Math.log(s);
        }

        return logq;
    }

    private double logFactorial(int n) {
        if (n == 0) {
            return 0;
        }

        double rValue = 0;

        for (int i = n; i > 0; i--) {
            rValue += Math.log(i);
        }
        return rValue;
    }


    private Double traverse(Tree tree, NodeRef node) {
        Double date;
        if (tree.isExternal(node)) {
            date = tree.getNodeHeight(node);
            if (!intervals.keySet().contains(date)) {
                throw new RuntimeException("Tip date not found");
            }

        } else {
            Double date1 = traverse(tree, tree.getChild(node, 0));
            Double date2 = traverse(tree, tree.getChild(node, 1));
            date = (date1 > date2 ? date1 : date2);
            if (!tree.isRoot(node)) {
                intervals.put(date, intervals.get(date) + 1);
            }
        }

        return date;
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
        return "nodeHeightScale";
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setAdaptableParameterValue(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    @Override
    public String getAdaptableParameterName() {
        return "scaleFactor";
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return NODE_HEIGHT_SCALE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AdaptationMode mode = AdaptationMode.parseMode(xo);

            final double weight = xo.getDoubleAttribute(WEIGHT);
            final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
            final boolean scaleAll = xo.getBooleanAttribute(SCALE_ALL);

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            NodeHeightScaleOperator operator = new NodeHeightScaleOperator(treeModel, scaleFactor, scaleAll, mode);
            operator.setWeight(weight);
            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a nodeHeightScale operator on a given tree.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(SCALE_ALL, true),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(TreeModel.class),
        };

    };

    public String toString() {
        return "nodeHeightScaleOperator(" + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private double scaleFactor = 0.5;
}