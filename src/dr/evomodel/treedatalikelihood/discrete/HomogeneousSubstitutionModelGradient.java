/*
 * HomogeneousSubstitutionModelGradient.java
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

import beagle.BeagleException;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.DifferentiableSubstitutionModel;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.WrtFrequencyParameter;
import dr.evomodel.substmodel.WrtRateParameter;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.ModelDependency;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Analytic gradient of log L with respect to all dimensions of one homogeneous substitution-model parameter.
 * For stationary frequencies pi it is g_k = B_k + R_k: the branch term through dQ~/dpi_k plus the root term
 * sum_s w_s p_root,s,k / L_s. Only the tangent part of g is meaningful on the simplex.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class HomogeneousSubstitutionModelGradient implements GradientWrtParameterProvider, Reportable {

    public enum Backend {
        AUTO("auto"),
        ADJOINT("adjoint"),
        PER_DIMENSION("perDimension");

        private final String label;

        Backend(String label) {
            this.label = label;
        }

        public String getLabel() {
            return label;
        }

        public static Backend parse(String label) {
            for (Backend backend : values()) {
                if (backend.label.equalsIgnoreCase(label)) {
                    return backend;
                }
            }
            throw new IllegalArgumentException("Unknown backend '" + label + "'; use auto, adjoint or perDimension");
        }
    }

    private static final Logger LOGGER = Logger.getLogger("dr.evomodel.treedatalikelihood.discrete");

    private final TreeDataLikelihood treeDataLikelihood;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final DifferentiableSubstitutionModel substitutionModel;
    private final Parameter parameter;
    private final List<WrtParameter> wrts;
    private final boolean includeBranch;
    private final boolean includeRoot;
    private final boolean isFrequency;
    private final Double tolerance;

    private BranchTermEngine engine;
    private BranchTermEngine perDimension;
    private boolean adjointDisabled = false;
    private boolean fallbackLogged = false;

    public HomogeneousSubstitutionModelGradient(TreeDataLikelihood treeDataLikelihood,
                                     Parameter parameter,
                                     Backend backend,
                                     Double tolerance) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.parameter = parameter;
        this.tolerance = tolerance;
        this.likelihoodDelegate = getBeagleDelegate(treeDataLikelihood);

        final BranchModel branchModel = likelihoodDelegate.getBranchModel();
        if (branchModel.getSubstitutionModels().size() != 1 || branchModel.requiresMatrixConvolution()) {
            throw new IllegalArgumentException("A substitution-model gradient needs a single substitution model"
                    + " without matrix convolution");
        }
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(0);
        if (!(model instanceof DifferentiableSubstitutionModel)) {
            throw new IllegalArgumentException(model.getId() + " is not differentiable");
        }
        this.substitutionModel = (DifferentiableSubstitutionModel) model;

        final FrequencyModel branchFrequencies = model.getFrequencyModel();
        final FrequencyModel rootFrequencies = branchModel.getRootFrequencyModel();
        final boolean branchHolds = frequencyParameterOf(branchFrequencies) == parameter;
        this.includeRoot = frequencyParameterOf(rootFrequencies) == parameter;
        this.isFrequency = branchHolds || includeRoot;

        this.wrts = new ArrayList<>();
        if (isFrequency) {
            requirePlainFrequencyModel(branchHolds ? branchFrequencies : null);
            requirePlainFrequencyModel(includeRoot ? rootFrequencies : null);
            auditHolders(branchHolds ? branchFrequencies : null, includeRoot ? rootFrequencies : null);
            if (branchHolds) {
                for (int k = 0; k < parameter.getDimension(); ++k) {
                    WrtParameter wrt = substitutionModel.factory(parameter, k);
                    if (!(wrt instanceof WrtFrequencyParameter)) {
                        throw new IllegalArgumentException(model.getId() + " has no frequency differential");
                    }
                    wrts.add(wrt);
                }
            }
        } else {
            for (int k = 0; k < parameter.getDimension(); ++k) {
                wrts.add(substitutionModel.factory(parameter, k));
            }
        }
        this.includeBranch = !wrts.isEmpty();

        if (includeBranch) {
            this.engine = resolveEngine(backend);
        }
    }

    /**
     * @return true if parameter is the stationary or root frequency parameter of the likelihood's single model
     */
    public static boolean isFrequencyParameter(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            return false;
        }
        BranchModel branchModel = ((BeagleDataLikelihoodDelegate) delegate).getBranchModel();
        if (frequencyParameterOf(branchModel.getRootFrequencyModel()) == parameter) {
            return true;
        }
        for (SubstitutionModel model : branchModel.getSubstitutionModels()) {
            if (frequencyParameterOf(model.getFrequencyModel()) == parameter) {
                return true;
            }
        }
        return false;
    }

    /**
     * @return true if parameter sets a relative rate of the likelihood's single substitution model through the generic
     * rate differential (e.g. TN93 kappa1/kappa2, GTR rates), which only the adjoint engine needs
     */
    public static boolean isRateParameter(TreeDataLikelihood treeDataLikelihood, Parameter parameter) {
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            return false;
        }
        List<SubstitutionModel> models = ((BeagleDataLikelihoodDelegate) delegate).getBranchModel().getSubstitutionModels();
        if (models.size() != 1 || !(models.get(0) instanceof DifferentiableSubstitutionModel)) {
            return false;
        }
        try {
            return ((DifferentiableSubstitutionModel) models.get(0)).factory(parameter, 0) instanceof WrtRateParameter;
        } catch (RuntimeException e) {
            return false;
        }
    }

    private static Parameter frequencyParameterOf(FrequencyModel frequencyModel) {
        if (frequencyModel == null) {
            return null;
        }
        try {
            return frequencyModel.getFrequencyParameter();
        } catch (RuntimeException e) {
            // derived frequency models (e.g. Markov-modulated) have no single frequency parameter
            return null;
        }
    }

    private static BeagleDataLikelihoodDelegate getBeagleDelegate(TreeDataLikelihood treeDataLikelihood) {
        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(treeDataLikelihood.getId() + " must use a single BEAGLE delegate");
        }
        return (BeagleDataLikelihoodDelegate) delegate;
    }

    private static void requirePlainFrequencyModel(FrequencyModel frequencyModel) {
        if (frequencyModel != null && frequencyModel.getClass() != FrequencyModel.class) {
            throw new UnsupportedOperationException("Frequency gradients are not supported for "
                    + frequencyModel.getClass().getSimpleName());
        }
    }

    private void auditHolders(FrequencyModel branchFrequencies, FrequencyModel rootFrequencies) {
        for (Model holder : ModelDependency.findHolders(treeDataLikelihood, parameter)) {
            if (holder != branchFrequencies && holder != rootFrequencies) {
                throw new IllegalArgumentException(parameter.getId() + " also enters " + holder.getId()
                        + "; its gradient would be missing that dependence");
            }
        }
    }

    private BranchTermEngine resolveEngine(Backend backend) {
        final String adjointReason = AdjointEngine.whyIneligible(likelihoodDelegate);
        switch (backend) {
            case ADJOINT:
                if (adjointReason != null) {
                    throw new IllegalArgumentException("backend=\"adjoint\" is unavailable: " + adjointReason);
                }
                return newAdjointEngine();
            case PER_DIMENSION:
                return newPerDimensionEngine();
            case AUTO:
            default:
                if (adjointReason == null) {
                    return newAdjointEngine();
                }
                LOGGER.info("Gradient wrt " + parameter.getId() + " uses the per-dimension path: " + adjointReason);
                return newPerDimensionEngine();
        }
    }

    private BranchTermEngine newAdjointEngine() {
        return new AdjointEngine(treeDataLikelihood, likelihoodDelegate, substitutionModel, wrts, stateCount());
    }

    private BranchTermEngine newPerDimensionEngine() {
        // The engine reuses traits by name, so the name must identify the parameter object: ids may be absent or shared
        return new PerDimensionEngine(treeDataLikelihood, likelihoodDelegate, substitutionModel, wrts,
                "homogeneousSubstitutionModelGradient:" + parameter.getId()
                        + "@" + Integer.toHexString(System.identityHashCode(parameter)));
    }

    private int stateCount() {
        return ((SubstitutionModel) substitutionModel).getDataType().getStateCount();
    }

    private double[] getBranchGradient() {
        if (!(engine instanceof AdjointEngine)) {
            return engine.getBranchGradient();
        }
        if (adjointDisabled) {
            return getPerDimension().getBranchGradient();
        }
        try {
            return engine.getBranchGradient();
        } catch (UnsatisfiedLinkError | BeagleException e) {
            // the BEAGLE library cannot run the adjoint: stay on the per-dimension path for this session
            logFallback(e.toString());
            adjointDisabled = true;
        } catch (AdjointEngine.NonFiniteAdjointException e) {
            // only this evaluation falls back
            logFallback(e.getMessage());
        }
        return getPerDimension().getBranchGradient();
    }

    private BranchTermEngine getPerDimension() {
        if (perDimension == null) {
            perDimension = newPerDimensionEngine();
        }
        return perDimension;
    }

    private void logFallback(String reason) {
        if (!fallbackLogged) {
            System.err.println("Gradient wrt " + parameter.getId() + " falls back to the per-dimension path: "
                    + reason);
            fallbackLogged = true;
        }
    }

    public double[] getBranchTerm() {
        treeDataLikelihood.getLogLikelihood();
        return includeBranch ? getBranchGradient() : new double[parameter.getDimension()];
    }

    public double[] getRootTerm() {
        final double[] root = new double[parameter.getDimension()];
        if (includeRoot) {
            treeDataLikelihood.getLogLikelihood();
            likelihoodDelegate.getRootFrequencyDerivatives(treeDataLikelihood.getTree().getRoot().getNumber(), root);
        }
        return root;
    }

    public String getEngineName() {
        if (!includeBranch) {
            return "root only";
        }
        return adjointDisabled ? "perDimension (adjoint disabled)" : engine.getName();
    }

    @Override
    public double[] getGradientLogDensity() {
        final double[] gradient = getBranchTerm();
        final double[] root = getRootTerm();
        for (int k = 0; k < gradient.length; ++k) {
            gradient[k] += root[k];
        }
        return gradient;
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
        return parameter.getDimension();
    }

    @Override
    public double getSimplexSum() {
        return isFrequency ? 1.0 : Double.NaN;
    }

    @Override
    public String getReport() {
        final double[] branch = getBranchTerm();
        final double[] root = getRootTerm();
        final double[] pi = parameter.getParameterValues();

        double piRoot = 0.0;
        double piGradient = 0.0;
        for (int k = 0; k < pi.length; ++k) {
            piRoot += pi[k] * root[k];
            piGradient += pi[k] * (branch[k] + root[k]);
        }

        StringBuilder sb = new StringBuilder("homogeneousSubstitutionModelGradient." + parameter.getId()
                + " (engine: " + getEngineName() + "; terms:" + (includeBranch ? " branch" : "")
                + (includeRoot ? " root" : "") + ")\n");
        sb.append("branch: ").append(new dr.math.matrixAlgebra.Vector(branch)).append("\n");
        sb.append("root: ").append(new dr.math.matrixAlgebra.Vector(root)).append("\n");
        if (isFrequency) {
            sb.append("pi.root: ").append(piRoot).append("\n");
            sb.append("pi.g: ").append(piGradient).append("\n");
        }
        sb.append(GradientWrtParameterProvider.getReportAndCheckForError(this,
                Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, tolerance));
        return sb.toString();
    }
}
