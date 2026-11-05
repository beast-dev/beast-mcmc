/*
 * NodeHeightTransform.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.continuous.hmc.NodeHeightTransformParser;
import dr.inference.model.Parameter;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.Transform;
import dr.xml.Reportable;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class NodeHeightTransform extends Transform.MultivariateTransform implements Reportable, Citable {

    private AbstractNodeHeightTransformDelegate nodeHeightTransformDelegate;
    private TreeModel tree;
    private Parameter nodeHeights;

    public NodeHeightTransform(Parameter ratios,
                               TreeModel tree,
                               BranchRateModel branchrateModel,
                               boolean withRoot) {
        super(tree.getInternalNodeCount(), withRoot ? tree.getInternalNodeCount() : tree.getInternalNodeCount() - 1);
        this.tree = tree;
        if (withRoot) {
            this.nodeHeights = new NodeHeightProxyParameter("nodeHeightToFullRatiosProxyParameter", tree, true);
            this.nodeHeightTransformDelegate = new NodeHeightToRatiosFullTransformDelegate(tree, nodeHeights, ratios, branchrateModel);
        } else {
            this.nodeHeights = new NodeHeightProxyParameter("nodeHeightToRatiosProxyParameter", tree, false);
            this.nodeHeightTransformDelegate = new NodeHeightToRatiosTransformDelegate(tree, nodeHeights, ratios, branchrateModel);
        }
    }

    @Deprecated
    public NodeHeightTransform(Parameter nodeHeights,
                               TreeModel tree,
                               GMRFSkyrideLikelihood skyrideLikelihood) {
        super(nodeHeights.getDimension());
        this.tree = tree;
        this.nodeHeightTransformDelegate = new NodeHeightToCoalescentIntervalsDelegate(tree, nodeHeights, skyrideLikelihood);
    }

    public AbstractNodeHeightTransformDelegate getNodeHeightTransformDelegate() {
        return nodeHeightTransformDelegate;
    }

    public Parameter getNodeHeights() {
        return nodeHeightTransformDelegate.getNodeHeights();
    }

    public Parameter getParameter() {
        return nodeHeightTransformDelegate.getParameter();
    }

    @Override
    public double transform(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    protected double[] transform(double[] values) {
        return nodeHeightTransformDelegate.transform(values);
    }

    @Override
    public double inverse(double value) {
        throw new RuntimeException("Should not be called.");
    }

    @Override
    protected double[] inverse(double[] values) {
        return nodeHeightTransformDelegate.inverse(values);
    }

    @Override
    public boolean isInInteriorDomain(double[] values) {
        return true; //TODO Is it necessary to check entry values ?
    }

    @Override
    public double[] inverse(double[] values, int from, int to, double sum) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradient(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[] gradientInverse(double[] values, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public String getTransformName() {
        return NodeHeightTransformParser.NAME;
    }

    @Override
    protected double getLogJacobian(double[] values) {
        return nodeHeightTransformDelegate.getLogJacobian(values);
    }

    @Override
    protected double[] updateGradientLogDensity(double[] gradient, double[] value) {
        return nodeHeightTransformDelegate.updateGradientLogDensity(gradient, value);
    }

    @Override
    public double[] updateGradientUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        return nodeHeightTransformDelegate.updateGradientUnWeightedLogDensity(gradient, value, from, to);
    }

    @Override
    public double[] updateGradientInverseUnWeightedLogDensity(double[] gradient, double[] value, int from, int to) {
        return nodeHeightTransformDelegate.updateGradientUnWeightedLogDensity(gradient, inverse(value, from, to), from, to);
    }

    @Override
    public String getReport() {
        return nodeHeightTransformDelegate.getReport();
    }

    @Override
    protected double[] getGradientLogJacobianInverse(double[] values) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double[][] computeJacobianMatrixInverse(double[] values) {
        throw new RuntimeException("Not yet implemented!");
    }

    // the ratio transform is cited; the coalescent-interval transform is not
    @Override
    public Citation.Category getCategory() {
        return Citation.Category.FRAMEWORK;
    }

    @Override
    public String getDescription() {
        return "Using the ratio transformation of internal node heights";
    }

    @Override
    public List<Citation> getCitations() {
        return nodeHeightTransformDelegate instanceof NodeHeightToRatiosTransformDelegate ?
                Collections.singletonList(RATIO_CITATION) : Collections.<Citation>emptyList();
    }

    public static final Citation RATIO_CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author("AA", "Fisher"),
                    new Author("S", "Su"),
                    new Author("JL", "Thorne"),
                    new Author("B", "Potter"),
                    new Author("P", "Lemey"),
                    new Author("G", "Baele"),
                    new Author("MA", "Suchard"),
            },
            "Scalable Bayesian divergence time estimation with ratio transformations",
            2023,
            "Systematic Biology",
            72,
            1136,
            1153,
            "10.1093/sysbio/syad039");

    /**
     * The ratio transform composed with the map of the ratios to the real line, cited as the ratio transform.
     */
    public static final class RealLine extends Transform.ComposeMultivariable implements Citable {

        private final NodeHeightTransform ratios;

        public RealLine(Transform.MultivariableTransform outer, NodeHeightTransform ratios) {
            super(outer, ratios);
            this.ratios = ratios;
        }

        @Override
        public Citation.Category getCategory() {
            return ratios.getCategory();
        }

        @Override
        public String getDescription() {
            return ratios.getDescription();
        }

        @Override
        public List<Citation> getCitations() {
            return ratios.getCitations();
        }
    }
}
