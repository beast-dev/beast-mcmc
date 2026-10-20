/*
 * HomogeneousSubstitutionParameterGradient.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodel.substmodel.DifferentialMassProvider.Mode;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class HomogeneousSubstitutionParameterGradient implements GradientWrtParameterProvider, Reportable {
    //TODO: create an AbstractClass to unify this and the BranchSubstitutionParameterGradient

    private final Parameter parameter;
    private final TreeDataLikelihood treeDataLikelihood;
    private final List<TreeTrait> treeTraitProviderList = new ArrayList<>();
    private final Tree tree;
    private final Mode mode;
    private final Integer dim;
    private final Double tolerance;

    public HomogeneousSubstitutionParameterGradient(String traitName,
                                                    TreeDataLikelihood treeDataLikelihood,
                                                    Parameter parameter,
                                                    BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                    Mode mode,
                                                    Double tolerance) {
        this(traitName, treeDataLikelihood, parameter, likelihoodDelegate, null, mode, tolerance);
    }

    /**
     * @param tolerance report check tolerance; null only prints the comparison with the numerical gradient
     */
    public HomogeneousSubstitutionParameterGradient(String traitName,
                                                    TreeDataLikelihood treeDataLikelihood,
                                                    Parameter parameter,
                                                    BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                    Integer dim,
                                                    Mode mode,
                                                    Double tolerance) {
        this.parameter = parameter;
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.mode = mode;
        this.dim = dim;
        this.tolerance = tolerance;

        final String name = BranchSubstitutionParameterDelegate.getName(traitName);
        if (treeDataLikelihood.getTreeTrait(name) != null) {
            throw new IllegalArgumentException("Trait name '" + traitName + "' is already used on "
                    + treeDataLikelihood.getId() + "; give each substitution-parameter gradient a unique traitName");
        }

        if (likelihoodDelegate.getBranchModel().getSubstitutionModels().size() != 1) {
            throw new RuntimeException("Homogeneous process should contain only one substitution model!");
        }
        SubstitutionModel model = likelihoodDelegate.getBranchModel().getSubstitutionModels().get(0);
        if (isFrequencyParameter(model, parameter)) {
            throw new IllegalArgumentException("Parameter " + parameter.getId() + " holds the stationary frequencies of "
                    + model.getId() + "; this gradient omits the root term and the simplex constraint."
                    + " Use the K-dimensional frequency gradient instead");
        }
        DifferentiableSubstitutionModel substitutionModel = (DifferentiableSubstitutionModel) model;

        if (dim != null) {
            if (dim < 0 || dim >= parameter.getDimension()) {
                throw new IllegalArgumentException("Dimension index " + dim + " is out of range [0, "
                        + parameter.getDimension() + ")");
            }
            buildTreeTraitProviderList(dim, substitutionModel, likelihoodDelegate, traitName, name);
        } else {
            for (int i = 0; i < parameter.getDimension(); i++) {
                buildTreeTraitProviderList(i, substitutionModel, likelihoodDelegate, traitName, name);
            }
        }
    }

    private static boolean isFrequencyParameter(SubstitutionModel model, Parameter parameter) {
        FrequencyModel frequencyModel = model.getFrequencyModel();
        if (frequencyModel == null) {
            return false;
        }
        try {
            return frequencyModel.getFrequencyParameter() == parameter;
        } catch (RuntimeException e) {
            // derived frequency models (e.g. Markov-modulated) have no single frequency parameter
            return false;
        }
    }

    public void buildTreeTraitProviderList(int i, DifferentiableSubstitutionModel substitutionModel,
                                           BeagleDataLikelihoodDelegate likelihoodDelegate,
                                           String traitName, String name){
        List<DifferentialMassProvider> differentialMassProviderList = new ArrayList<>();
        DifferentialMassProvider.DifferentialWrapper.WrtParameter wrtParameter = substitutionModel.factory(parameter, i);
        DifferentialMassProvider differentialMassProvider = new DifferentialMassProvider.DifferentialWrapper(
                substitutionModel, wrtParameter, mode);
        differentialMassProviderList.add(differentialMassProvider);
        BranchDifferentialMassProvider branchDifferentialMassProvider =
                new BranchDifferentialMassProvider(null, differentialMassProviderList);

        ProcessSimulationDelegate gradientDelegate = new BranchSubstitutionParameterDelegate(traitName,
                treeDataLikelihood.getTree(),
                likelihoodDelegate,
                treeDataLikelihood.getBranchRateModel(),
                branchDifferentialMassProvider);

        TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, gradientDelegate);
        treeDataLikelihood.addTraits(traitProvider.getTreeTraits());
        treeTraitProviderList.add(treeDataLikelihood.getTreeTrait(name));
    }



    @Override
    public Likelihood getLikelihood() {
        return treeDataLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        if (dim == null) {
            return parameter.getDimension();
        } else {
            return 1;
        }
    }

    @Override
    public double[] getGradientLogDensity() {
        if (treeTraitProviderList.size() == 1) {
            double result = 0.0;
            double[] gradient = (double[]) treeTraitProviderList.get(0).getTrait(tree, null);
            for (int i = 0; i < gradient.length; i++) {
                result += gradient[i]; // sum over all branches
            }
            return new double[]{result};
        } else {
            double[] result = new double[treeTraitProviderList.size()];
            for (int i = 0; i < treeTraitProviderList.size(); i++) {
                double[] gradient = (double[]) treeTraitProviderList.get(i).getTrait(tree, null);
                for (int j = 0; j < gradient.length; j++) {
                    result[i] += gradient[j]; // sum over all branches
                }
            }
            return result;
        }

    }

    @Override
    public String getReport() {
        return mode.getReport() + " " + GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, tolerance);
    }
}
