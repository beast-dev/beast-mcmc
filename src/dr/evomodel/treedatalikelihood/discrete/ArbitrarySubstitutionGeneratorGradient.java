/*
 * HomogeneousSubstitutionParameterGradient.java
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

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.coalescent.OldGMRFSkyrideLikelihood;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider;
import dr.evomodel.substmodel.OldGLMSubstitutionModel;
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

/**
 * @author Marc A. Suchard
 */
public class ArbitrarySubstitutionGeneratorGradient implements GradientWrtParameterProvider, Reportable {
    //TODO: create an AbstractClass to unify this and the BranchSubstitutionParameterGradient

    private final Parameter parameter;
    private final TreeDataLikelihood treeDataLikelihood;
    private final TreeTrait treeTraitProvider;
    private final Tree tree;
    private final OldGLMSubstitutionModel substitutionModel;
    private final int stateCount;

    public ArbitrarySubstitutionGeneratorGradient(String traitName,
                                                  TreeDataLikelihood treeDataLikelihood,
//                                                  Parameter parameter,
                                                  BeagleDataLikelihoodDelegate likelihoodDelegate,
                                                  OldGLMSubstitutionModel substitutionModel) {
//        this.parameter = parameter;
        this.treeDataLikelihood = treeDataLikelihood;
        this.tree = treeDataLikelihood.getTree();
        this.substitutionModel = substitutionModel;
        this.parameter = null;
        this.stateCount = substitutionModel.getDataType().getStateCount();

        final String name = BranchSubstitutionParameterDelegate.getName(traitName);
        TreeTrait test = treeDataLikelihood.getTreeTrait(name);

        if (test == null) {
            
            DifferentialMassProvider differentialMassProvider = time -> {

                double[] differential = new double[stateCount * stateCount];
                double[] generator = new double[stateCount * stateCount];

                substitutionModel.getInfinitesimalMatrix(generator);
                differential[0] = generator[0];
                differential[1] = generator[1];
                differential[2] = generator[2]; // Comment out for irreversible case
                differential[3] = generator[3]; // Comment out for irreversible case

                return differential;
            };

            List<DifferentialMassProvider> differentialMassProviderList = new ArrayList<DifferentialMassProvider>();
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
        }

        this.treeTraitProvider = treeDataLikelihood.getTreeTrait(name);
        assert(treeTraitProvider != null);
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
        return 1;
    }

    @Override
    public double[] getGradientLogDensity() {
        double result = 0.0;

        double[] gradient = (double[]) treeTraitProvider.getTrait(tree, null);

        for (int i = 0; i < gradient.length; i++) {
            result += gradient[i];
        }

        return new double[]{result};
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, null);
    }

    private final double tolerance = 1E-2;
}
