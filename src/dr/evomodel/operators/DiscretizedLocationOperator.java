/*
 * DiscretizedLocationOperator.java
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

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.inference.operators.AbstractCoercableOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.MCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.awt.geom.Point2D;
import java.util.*;
import java.util.logging.Logger;


/**
 * @author Marc Suchard
 */
public class DiscretizedLocationOperator extends AbstractCoercableOperator {

    public static final String GIBBS_OPERATOR = "discretizedLocationOperator";
    public static final String INTERNAL_ONLY = "onlyInternalNodes";
    public static final String DISK = "neighborhoodSize";
    public static final String RANDOMIZE = "randomize";

    public DiscretizedLocationOperator(AbstractMultivariateTraitLikelihood traitModel, boolean onlyInternalNodes, int disk, CoercionMode mode) {
        super(mode);
        this.treeModel = traitModel.getTreeModel();
        this.traitName = traitModel.getTraitName();
        this.onlyInternalNodes = onlyInternalNodes;
        allLocations = makeLocationList();
        nearestNeighborMap = makeNearestNeighborMap();
        this.disk = disk;
        this.autoOptimize = convertToAutoOptimizeValue(disk);

        if (disk > allLocations.size() - 2)
            throw new RuntimeException("Neighborhood size is too large");

        printInfo();
    }

    private Map<Point2D, List<WeightedPoint2D>> makeNearestNeighborMap() {
        Map<Point2D, List<WeightedPoint2D>> map = new HashMap<Point2D, List<WeightedPoint2D>>();

        for (Point2D location : allLocations) {

            List<WeightedPoint2D> weightedNeighbors = new ArrayList<WeightedPoint2D>();
            for (Point2D neighbor : allLocations) {
                double distance = location.distance(neighbor);
                if (distance > 0)
                    weightedNeighbors.add(
                            new WeightedPoint2D(neighbor.getX(), neighbor.getY(), distance)
                    );
            }
            Collections.sort(weightedNeighbors);

            map.put(location, weightedNeighbors);
        }

//        for (Point2D location : map.keySet()) {
//            List<WeightedPoint2D> neighbors = map.get(location);
//            System.err.println("Location: "+location+"\n");
//            System.err.println("\t");
//            int count = 0;
//            for (WeightedPoint2D neighbor : neighbors) {
//                count++;
//                if (count < 3)
//                    System.err.println(" "+neighbor);
//            }
//            System.err.println("\n");
//        }

        return map;
    }


    private void recursivelySetTrait(NodeRef node, double[] trait, NodeRef fromNode) {

        treeModel.setMultivariateTrait(node, traitName, trait);

        for (int i = 0; i < treeModel.getChildCount(node); i++) {
            NodeRef child = treeModel.getChild(node, i);
            if (child != fromNode && treeModel.getBranchLength(child) == 0) {
//                System.err.println("recursing down");
                recursivelySetTrait(child, trait, node);
            }
        }
        if (!treeModel.isRoot(node) && treeModel.getBranchLength(node) == 0) {
//            System.err.println("recursing up");
            recursivelySetTrait(treeModel.getParent(node), trait, node);
        }
    }


    public void randomizeNodes() {

        List<Point2D> listLocations = new ArrayList<Point2D>();
        listLocations.addAll(allLocations);

        for (int i = 0; i < treeModel.getInternalNodeCount(); i++) {

            NodeRef node = treeModel.getInternalNode(i);
            double[] trait = treeModel.getMultivariateNodeTrait(node, traitName);

            Point2D newPt = listLocations.get(MathUtils.nextInt(listLocations.size()));

            trait[0] = newPt.getX();
            trait[1] = newPt.getY();

            recursivelySetTrait(node, trait, null);

//            treeModel.setMultivariateTrait(node, traitName, trait);            
        }
        System.err.println("Done with randomization");
//        System.exit(-1);
    }

    private void printInfo() {
        StringBuffer sb = new StringBuffer();
        sb.append("\nCreating a discretized location sampler:\n");
        sb.append("\tTip count: " + treeModel.getExternalNodeCount() + "\n");
        sb.append("\tUnique locations: " + allLocations.size() + "\n");
        sb.append("\tNeighborhood size: " + disk + "\n");
        Logger.getLogger("dr.evomodel.operators").info(sb.toString());
    }

    private Set<Point2D> makeLocationList() {

        Set<Point2D> uniquePoints = new HashSet<Point2D>();

        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            double[] leafTrait = treeModel.getMultivariateNodeTrait(node, traitName);

            Point2D.Double point = new Point2D.Double(leafTrait[0], leafTrait[1]);
            if (!uniquePoints.contains(point)) {
                uniquePoints.add(point);
                savedPt = point;
            }
        }

        return uniquePoints;
    }

    public double doOperation() {

        NodeRef node;

        if (onlyInternalNodes)
            node = treeModel.getInternalNode(MathUtils.nextInt(
                    treeModel.getInternalNodeCount()));
        else
            node = treeModel.getNode(MathUtils.nextInt(
                    treeModel.getNodeCount()));

        double[] trait = treeModel.getMultivariateNodeTrait(node, traitName);
        Point2D currentPt = new Point2D.Double(trait[0], trait[1]);

        List<WeightedPoint2D> neighbors = nearestNeighborMap.get(currentPt);

        if (neighbors == null)
            throw new RuntimeException("Node location outside allowable values: " + currentPt);

//        Point2D newPt = neighbors.get(MathUtils.nextInt(disk));
        Point2D newPt = neighbors.get(MathUtils.nextInt(convertFromAutoOptimizeValue(autoOptimize)));

        trait[0] = newPt.getX();
        trait[1] = newPt.getY();

//        treeModel.setMultivariateTrait(node, traitName, trait);
        recursivelySetTrait(node, trait, null);

        return 0;

    }

    private int convertFromAutoOptimizeValue(double value) {
        return 1 + (int) Math.exp(autoOptimize);
    }

    private double convertToAutoOptimizeValue(int value) {
        return Math.log(value - 1);
    }

    public double getCoercableParameter() {
        return autoOptimize;
    }

    public void setCoercableParameter(double value) {
        autoOptimize = value;
    }

    public double getRawParameter() {
        return convertFromAutoOptimizeValue(autoOptimize);
    }

//    public double getScaleFactor() {
//        return scaleFactor;
//    }

    public double getTargetAcceptanceProbability() {
        return 0.50;
    }

    public final String getPerformanceSuggestion() {

//        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//        double targetProb = getTargetAcceptanceProbability();
//        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
//        double sf = OperatorUtils.optimizeScaleFactor(autoOptimize, prob, targetProb);
//        if (prob < getMinimumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else if (prob > getMaximumGoodAcceptanceLevel()) {
//            return "Try setting scaleFactor to about " + formatter.format(sf);
//        } else return "";
        return "I have no idea.";
    }


    public String getOperatorName() {
        return GIBBS_OPERATOR;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GIBBS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            boolean onlyInternalNodes = xo.getAttribute(INTERNAL_ONLY, true);
            int disk = xo.getAttribute(DISK, 4);


            AbstractMultivariateTraitLikelihood traitModel = (AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);

            DiscretizedLocationOperator operator = new DiscretizedLocationOperator(traitModel, onlyInternalNodes, disk, mode);
            operator.setWeight(weight);

            // Get node priors
//            for (int i = 0; i < xo.getChildCount(); i++) {
//                if (xo.getChild(i) instanceof MultivariateDistributionLikelihood) {
//                    MultivariateDistribution dist = ((MultivariateDistributionLikelihood) xo.getChild(i)).getDistribution();
//                    if (dist instanceof GeoSpatialDistribution) {
//                        GeoSpatialDistribution prior = (GeoSpatialDistribution) dist;
//                        String nodeLabel = prior.getLabel();
//                        TreeModel treeModel = traitModel.getTreeModel();
//
//                        // Get taxon node from tree
//                        int index = treeModel.getTaxonIndex(nodeLabel);
//                        if (index == -1) {
//                            throw new XMLParseException("taxon '" + nodeLabel + "' not found for geoSpatialDistribution element in traitGibbsOperator element");
//                        }
//                        NodeRef node = treeModel.getExternalNode(index);
////                        operator.setTaxonPrior(node, prior);
//                        System.err.println("Adding truncated prior for " + node);
//                    }
//                }
//            }

            boolean randomize = xo.getAttribute(RANDOMIZE, false);

            if (randomize)
                operator.randomizeNodes();

            return operator;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate Gibbs operator on traits for possible all nodes.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(INTERNAL_ONLY, true),
                new ElementRule(AbstractMultivariateTraitLikelihood.class),
                AttributeRule.newIntegerRule(DISK, true),
                AttributeRule.newBooleanRule(RANDOMIZE, true),
        };

    };

    public class WeightedPoint2D extends Point2D.Double implements Comparable {

        public double weight;

        public WeightedPoint2D(double x, double y, double weight) {
            super(x, y);
            this.weight = weight;
        }

        public int compareTo(Object o) {
            WeightedPoint2D pt = (WeightedPoint2D) o;
            if (weight > pt.weight)
                return 1;
            if (weight < pt.weight)
                return -1;
            return 0;
        }

        public String toString() {
            return super.toString() + "(" + weight + ")";
        }
    }

    private Map<Point2D, List<WeightedPoint2D>> nearestNeighborMap;
    private Set<Point2D> allLocations;
    private final MutableTreeModel treeModel;
    private String traitName;
    private double autoOptimize;
    private boolean onlyInternalNodes = true;
    private int disk = 4;

    private Point2D savedPt;

}
