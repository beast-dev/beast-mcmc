/*
 * RandomLocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.ScaledByTreeTimeBranchRateModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * @author Marc A. Suchard
 * @author Alexei J. Drummond
 * @author Alexander Fisher
 */

public class ScaledByTreeTimeBranchRateModel extends AbstractBranchRateModel implements DifferentiableBranchRates, Citable {

    private final TreeModel treeModel;
    private final BranchRateModel branchRateModel;
    private final Parameter meanRateParameter;

    private boolean scaleFactorKnown;
    private boolean storedScaleFactorKnown;

    private double scaleFactor;
    private double storedScaleFactor;

    private double branchTotal;
    private double storedBranchTotal;
    private double timeTotal;
    private double storedTimeTotal;

    public ScaledByTreeTimeBranchRateModel(TreeModel treeModel,
                                           BranchRateModel branchRateModel,
                                           Parameter meanRateParameter) {

        super(ScaledByTreeTimeBranchRateModelParser.TREE_TIME_BRANCH_RATES);

        this.treeModel = treeModel;
        this.branchRateModel = branchRateModel;
        this.meanRateParameter = meanRateParameter;

        addModel(treeModel);
        addModel(branchRateModel);

        if (meanRateParameter != null) {
            addVariable(meanRateParameter);
        }

        scaleFactorKnown = false;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        scaleFactorKnown = false;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        scaleFactorKnown = false;
        fireModelChanged();
    }

    protected void storeState() {
        storedScaleFactor = scaleFactor;
        storedScaleFactorKnown = scaleFactorKnown;

        storedBranchTotal = branchTotal;
        storedTimeTotal = timeTotal;
    }

    protected void restoreState() {
        scaleFactor = storedScaleFactor;
        scaleFactorKnown = storedScaleFactorKnown;

        branchTotal = storedBranchTotal;
        timeTotal = storedTimeTotal;
    }

    protected void acceptState() { }

    @Override
    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }

        return ((DifferentiableBranchRates)branchRateModel).getBranchRateDifferential(tree, node);
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }

        return ((DifferentiableBranchRates)branchRateModel).getBranchRateSecondDifferential(tree, node);
    }

    @Override
    public Parameter getRateParameter() {
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }
        return ((DifferentiableBranchRates)branchRateModel).getRateParameter();
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        if (!(branchRateModel instanceof DifferentiableBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }
        return ((DifferentiableBranchRates)branchRateModel).getParameterIndexFromNode(node);
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {

        if (!scaleFactorKnown) {
            calculateScaleFactor();
            scaleFactorKnown = true;
        }

//        // TODO Profile and possibly optimize
//        //   for example: with index -> (length, rate) map
//
//        double[] result = new double[treeModel.getNodeCount() - 1];
//
//        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
//            final NodeRef nodeI = treeModel.getNode(i);
//            if (!treeModel.isRoot(nodeI)) {
//                final int indexOne = getParameterIndexFromNode(nodeI);
//
//                final double crossTermNodeI = scaleFactor * scaleFactor / timeTotal * treeModel.getBranchLength(nodeI);
//
//                for (int j = 0; j < treeModel.getNodeCount(); ++j) {
//                    final NodeRef nodeJ = treeModel.getNode(j);
//                    if (!treeModel.isRoot(nodeJ)) {
//                        final int indexTwo = getParameterIndexFromNode(nodeJ);
//
//                        result[indexOne] -=  crossTermNodeI * branchRateModel.getBranchRate(treeModel, nodeJ) * gradient[indexTwo];
//
//                        if (indexOne == indexTwo) {
//                            result[indexOne] += scaleFactor * gradient[indexTwo];
//                        }
//
//                    }
//                }
//            }
//        }
//
//        return result;

        final int dim = gradient.length;
        final DenseMatrix64F Jacobian = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F gradVector = new DenseMatrix64F(dim, 1);
        final DenseMatrix64F scaledGradient = new DenseMatrix64F(dim, 1);

        // compute Jacobian matrix
        double tempTotal;
        int rootNodeIndex = treeModel.getRoot().getNumber(); // to ignore rootNode

        for (int row = 0; row < rootNodeIndex; ++row) {
            NodeRef nodeJ = treeModel.getNode(row);
            for (int col = 0; col < rootNodeIndex; ++col) {
                NodeRef nodeI = treeModel.getNode(col);
                tempTotal = getTempTotal(nodeI, nodeJ);

                Jacobian.unsafe_set(row, col, tempTotal);
            }
        }

        if(rootNodeIndex < dim) {
            for (int row = rootNodeIndex + 1; row < dim + 1; ++row) {
                NodeRef nodeJ = treeModel.getNode(row);
                for (int col = 0; col < rootNodeIndex; ++col) {
                    NodeRef nodeI = treeModel.getNode(col);
                    tempTotal = getTempTotal(nodeI, nodeJ);

                    Jacobian.unsafe_set(row - 1 , col - 1, tempTotal);
                }
            }
        }

        for(int row = 0; row<dim; row++) {
            // add to diagonals & setup vector for multiplication
            tempTotal = Jacobian.unsafe_get(row, row);
            Jacobian.unsafe_set(row, row, tempTotal + scaleFactor);
            gradVector.set(row, 0, gradient[row]);
        }


        CommonOps.mult(Jacobian, gradVector, scaledGradient);

        for (int i = 0; i < dim; ++i){
            gradient[i] = scaledGradient.unsafe_get(0,i);
        }

        return gradient;
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert tree == treeModel;

        if (!scaleFactorKnown) {
            calculateScaleFactor();
            scaleFactorKnown = true;
        }

        return scaleFactor * branchRateModel.getBranchRate(tree, node);
    }

    private void calculateScaleFactor() {

        double timeTotal = 0.0;
        double branchTotal = 0.0;

        for (int i = 0; i < treeModel.getNodeCount(); i++) {
            NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {

                double branchTime = treeModel.getBranchLength(node);

                double branchLength = branchTime * branchRateModel.getBranchRate(treeModel, node);

                timeTotal += branchTime;
                branchTotal += branchLength;
            }
        }

        double scaleFactor = timeTotal / branchTotal;

        if (meanRateParameter != null) {
            scaleFactor *= meanRateParameter.getParameterValue(0);
        }

        this.scaleFactor = scaleFactor;
        this.branchTotal = branchTotal;
        this.timeTotal = timeTotal;
    }

    private double getTempTotal(NodeRef nodeI, NodeRef nodeJ){
        double total = -branchRateModel.getBranchRate(treeModel, nodeI) * treeModel.getBranchLength(nodeJ) * scaleFactor;
        total /= branchTotal;
        return total;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        String description =
                (branchRateModel instanceof Citable) ?
                        ((Citable) branchRateModel).getDescription() :
                        "Unknown clock model";

        description += " with scaling-by-tree-time";
        return description;
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> list = 
                (branchRateModel instanceof Citable) ?
                        new ArrayList<>(((Citable) branchRateModel).getCitations()) :
                        new ArrayList<>();
        list.add(RandomLocalClockModel.CITATION);
        return list;
    }
}