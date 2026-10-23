/*
 * SiteRateModelGradient.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.siteratemodel.DifferentiableSiteRateModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * Analytic gradient of log L with respect to a parameter of the discretized site-rate model (gamma shape, relative
 * rate or proportion of invariant sites), through the category rates r_c and proportions w_c:
 * <p>
 * d log L / d theta = sum_c (d log L / d r_c)(d r_c / d theta) + sum_c (d log L / d w_c)(d w_c / d theta).
 * <p>
 * The rate term is sum_e length_e g_e, where g_e comes from one pre-order pass with (d r_c / d theta) Q~ as the
 * per-category derivative matrices ({@link SiteRateDerivativeDelegate}). The weight term uses
 * d log L / d w_c = sum_s w_s L_{s,c} / L_s from the root partials. The category Jacobian comes from the site-rate
 * model itself ({@link DifferentiableSiteRateModel}), so it follows that model's normalization.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class SiteRateModelGradient implements GradientWrtParameterProvider, Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final DifferentiableSiteRateModel siteRateModel;
    private final Parameter parameter;
    private final Tree tree;
    private final TreeTrait rateTrait;
    private final Double tolerance;

    public SiteRateModelGradient(TreeDataLikelihood treeDataLikelihood, Parameter parameter, Double tolerance) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.parameter = parameter;
        this.tolerance = tolerance;
        this.tree = treeDataLikelihood.getTree();

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(treeDataLikelihood.getId() + " must use a single BEAGLE delegate");
        }
        this.likelihoodDelegate = (BeagleDataLikelihoodDelegate) delegate;

        String reason = whyIneligible(likelihoodDelegate);
        if (reason != null) {
            throw new IllegalArgumentException(reason);
        }

        final SiteRateModel model = likelihoodDelegate.getSiteRateModel();
        if (!(model instanceof DifferentiableSiteRateModel)
                || !((DifferentiableSiteRateModel) model).isDifferentiableWrtParameter(parameter)) {
            throw new IllegalArgumentException(parameter.getId() + " is not a differentiable parameter of the site model"
                    + " of " + treeDataLikelihood.getId() + " (supported: gamma shape, proportion of invariant sites, relative rate)");
        }
        if (parameter.getDimension() != 1) {
            throw new IllegalArgumentException(parameter.getId() + " must be a scalar");
        }
        this.siteRateModel = (DifferentiableSiteRateModel) model;

        // The trait is keyed by the parameter object, since ids may be absent or shared
        final String name = "siteRateModelGradient:" + parameter.getId()
                + "@" + Integer.toHexString(System.identityHashCode(parameter));
        final String traitName = SiteRateDerivativeDelegate.gradientTraitName(name);
        if (treeDataLikelihood.getTreeTrait(traitName) == null) {
            SiteRateDerivativeDelegate rateDelegate = new SiteRateDerivativeDelegate(name, tree, likelihoodDelegate,
                    treeDataLikelihood.getBranchRateModel(), siteRateModel, parameter);
            treeDataLikelihood.addTraits(new ProcessSimulation(treeDataLikelihood, rateDelegate).getTreeTraits());
        }
        this.rateTrait = treeDataLikelihood.getTreeTrait(traitName);
    }

    /**
     * @return null if the gradient applies to this likelihood, otherwise the reason it does not
     */
    static String whyIneligible(BeagleDataLikelihoodDelegate likelihoodDelegate) {
        if (!likelihoodDelegate.isUsePreOrder() || !likelihoodDelegate.getPreOrderSettings().isBranchRateDerivative()) {
            return "site-model gradients need usePreOrder=\"true\" and branchRateDerivative=\"true\"";
        }
        if (likelihoodDelegate.getPreOrderSettings().useHierarchicalRepresentation()) {
            return "site-model gradients do not support the hierarchical representation";
        }
        BranchModel branchModel = likelihoodDelegate.getBranchModel();
        if (branchModel.getSubstitutionModels().size() != 1 || branchModel.requiresMatrixConvolution()) {
            return "site-model gradients need a single substitution model without matrix convolution";
        }
        return null;
    }

    /**
     * @return the part of the gradient that flows through the category rates
     */
    public double getRateTerm() {
        treeDataLikelihood.getLogLikelihood();
        final double[] dRates = new double[siteRateModel.getCategoryCount()];
        siteRateModel.getCategoryDifferentials(parameter, dRates, new double[dRates.length]);
        if (isZero(dRates)) {
            return 0.0;
        }

        final double[] perBranch = (double[]) rateTrait.getTrait(tree, null); // d log L / d length_e, non-root order
        final BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        double sum = 0.0;
        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final double length = branchRateModel.getBranchRate(tree, node) * tree.getBranchLength(node);
                sum += perBranch[v] * length;
                v++;
            }
        }
        return sum;
    }

    /**
     * @return the part of the gradient that flows through the category proportions
     */
    public double getWeightTerm() {
        treeDataLikelihood.getLogLikelihood();
        final int categoryCount = siteRateModel.getCategoryCount();
        final double[] dProportions = new double[categoryCount];
        siteRateModel.getCategoryDifferentials(parameter, new double[categoryCount], dProportions);
        if (isZero(dProportions)) {
            return 0.0;
        }

        final double[] weightDerivatives = new double[categoryCount];
        likelihoodDelegate.getCategoryWeightDerivatives(tree.getRoot().getNumber(), weightDerivatives);
        double sum = 0.0;
        for (int c = 0; c < categoryCount; ++c) {
            sum += weightDerivatives[c] * dProportions[c];
        }
        return sum;
    }

    @Override
    public double[] getGradientLogDensity() {
        return new double[]{getRateTerm() + getWeightTerm()};
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

    // GammaFunction.incompleteGammaP (AS32, accuracy 1e-8) gives the gamma quantiles ~1e-9 jitter in the shape,
    // so the report's finite differences step above it
    @Override
    public double getNumericGradientStepSize() {
        return siteRateModel.hasQuantileNoise(parameter) ? QUANTILE_NOISE_STEP_SIZE
                : GradientWrtParameterProvider.super.getNumericGradientStepSize();
    }

    private static final double QUANTILE_NOISE_STEP_SIZE = 1E-4;

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder("siteRateModelGradient." + parameter.getId() + "\n");
        sb.append("rate term: ").append(getRateTerm()).append("\n");
        sb.append("weight term: ").append(getWeightTerm()).append("\n");
        sb.append(GradientWrtParameterProvider.getReportAndCheckForError(this,
                parameter.getBounds().getLowerLimit(0), parameter.getBounds().getUpperLimit(0), tolerance));
        return sb.toString();
    }

    private static boolean isZero(double[] values) {
        for (double value : values) {
            if (value != 0.0) {
                return false;
            }
        }
        return true;
    }
}
