/*
 * DiscreteTraitBranchSubstitutionParameterGradient.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.*;
import dr.evomodel.branchmodel.BranchSpecificSubstitutionParameterBranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

import static dr.math.MachineAccuracy.SQRT_EPSILON;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class BranchSubstitutionParameterGradient
        implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    protected final TreeDataLikelihood treeDataLikelihood;
    protected final TreeTrait treeTraitProvider;
    protected final Tree tree;
    protected final boolean useHessian;

    protected final CompoundParameter branchParameter;
    protected final DifferentiableBranchRates branchRateModel;

    protected final Double nullableTolerance;
    private static final boolean DEBUG = true;

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    protected long getGradientLogDensityCount = 0;
    private final double smallGradientThreshold = 0.5;

    private BranchDifferentialMassProvider save;

    public BranchSubstitutionParameterGradient(String traitName,
                                               TreeDataLikelihood treeDataLikelihood,
                                               BeagleDataLikelihoodDelegate likelihoodDelegate,
                                               CompoundParameter branchParameter,
                                               DifferentiableBranchRates branchRateModel,
                                               Double tolerance,
                                               boolean useHessian,
                                               int dim) {
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.branchParameter = branchParameter;
        this.branchRateModel = branchRateModel;
        this.useHessian = useHessian;
        this.nullableTolerance = tolerance;

        String name = BranchSubstitutionParameterDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {

            BranchSpecificSubstitutionParameterBranchModel branchModel = (BranchSpecificSubstitutionParameterBranchModel) likelihoodDelegate.getBranchModel();

            List<DifferentialMassProvider> differentialMassProviderList = new ArrayList<>();
            for (int i = 0; i < tree.getNodeCount(); ++i) {
                NodeRef node = tree.getNode(i);
                if (!tree.isRoot(node)) {
                    DifferentiableSubstitutionModel substitutionModel = (DifferentiableSubstitutionModel) branchModel.getSubstitutionModel(node);
                    Parameter parameter = branchParameter.getParameter(node.getNumber());
                    DifferentialMassProvider.DifferentialWrapper.WrtParameter wrtParameter = substitutionModel.factory(parameter, dim);
                    differentialMassProviderList.add(new DifferentialMassProvider.DifferentialWrapper(substitutionModel, wrtParameter));
                }
            }

            BranchDifferentialMassProvider branchDifferentialMassProvider =
                    new BranchDifferentialMassProvider(branchRateModel, differentialMassProviderList);

            this.save = branchDifferentialMassProvider;

            ProcessSimulationDelegate gradientDelegate = new BranchSubstitutionParameterDelegate(traitName,
                    treeDataLikelihood.getTree(),
                    likelihoodDelegate,
                    treeDataLikelihood.getBranchRateModel(),
                    branchDifferentialMassProvider);
            TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
            treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
        }

        treeTraitProvider = treeDataLikelihood.getTreeTrait(name);
        assert (treeTraitProvider != null);

        int nTraits = treeDataLikelihood.getDataLikelihoodDelegate().getTraitCount();
        if (nTraits != 1) {
            throw new RuntimeException("Not yet implemented for >1 traits");
        }
    }


    @Override
    public double[] getDiagonalHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return branchRateModel.getRateParameter();
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] result = new double[getDimension()];

        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = branchRateModel.getParameterIndexFromNode(node);
                result[destinationIndex] = gradient[destinationIndex] *
                        branchRateModel.getBranchRateDifferential(tree, node);
            }
            // TODO Handle root node at most point
        }

        if (COUNT_TOTAL_OPERATIONS) {
            ++getGradientLogDensityCount;
        }

        return result;
    }

    protected double getChainGradient(Tree tree, NodeRef node) {
        if (branchRateModel instanceof ArbitraryBranchRates) {
            final double raw = getParameter().getParameterValue(branchRateModel.getParameterIndexFromNode(node));
            return branchRateModel.getTransform().differential(raw, tree, node);
        } else {
            return 1.0;
        }
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "BranchSubstitutionParameterGradientReport");
    }



    protected boolean valuesAreSufficientlyLarge(double[] vector) {
        for (double x : vector) {
            if (Math.abs(x) < SQRT_EPSILON * 1.2) {
                return false;
            }
        }

        return true;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, null);

//        BranchSpecificSubstitutionParameterBranchModel branchModel = (BranchSpecificSubstitutionParameterBranchModel)
//                ((BeagleDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate()).getBranchModel();
//
//                DifferentiableSubstitutionModel substitutionModel = (DifferentiableSubstitutionModel) branchModel.getSubstitutionModel(tree.getNode(0));
//
////                substitutionModel.getInfinitesimalDifferentialMatrix()
//        double[] differential = save.getDifferentialMassMatrixForBranch(tree.getNode(0), tree.getBranchLength(tree.getNode(0)));
//        double[] generator = new double[differential.length];
//        substitutionModel.getInfinitesimalMatrix(generator);
//        int len = substitutionModel.getFrequencyModel().getDataType().getStateCount();
//
//        checkCommutability(new WrappedMatrix.Raw(generator, 0, len, len), new WrappedMatrix.Raw(differential, 0, len, len));
//
//        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, nullableTolerance, smallGradientThreshold);
    }
}
