/*
 * MascotGradient.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DifferentiableBranchRates;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * HMC-facing gradient provider for one part (migration rates, population
 * sizes, or the clock rate) of a MASCOT likelihood's parameters. MascotCore
 * itself still returns one combined flat gradient per evaluation (see
 * {@link MascotLikelihood#getGradientLogDensity()}); this class exposes just
 * the requested part as an independently-addressable
 * {@link GradientWrtParameterProvider}, so migration rates, population
 * sizes, and the clock rate can each have their own HMC operator/step size,
 * prior, and BEAUti panel, matching how BASTA's XML already separates
 * migration/popSizes.
 * <p>
 * {@code Part.CLOCK_RATE} requires the likelihood's {@link BranchRateModel}
 * to implement {@link DifferentiableBranchRates}: MascotCore itself never
 * assumes anything beyond a plain per-lineage {@code double[]}, and this
 * class's {@code CLOCK_RATE} reduction (build a per-branch array, then call
 * {@link DifferentiableBranchRates#updateGradientLogDensity}) is the same
 * generic mechanism the rest of BEAST-X already uses for e.g.
 * {@code ArbitraryBranchRates}, {@code LocalClockModel}, etc -- so any
 * current or future {@code DifferentiableBranchRates} implementation works
 * here with no further changes, not just {@code StrictClockBranchRates}.
 */
public final class MascotGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    public enum Part {
        MIGRATION, POPULATION_SIZE, CLOCK_RATE
    }

    private final MascotLikelihood likelihood;
    private final Part part;
    private double numericGradientStepSize = GradientWrtParameterProvider.super.getNumericGradientStepSize();

    public MascotGradient(MascotLikelihood likelihood, Part part) {
        this.likelihood = likelihood;
        this.part = part;
        if (part == Part.MIGRATION) {
            String error = likelihood.getMigrationGradientCompatibilityError();
            if (error != null) {
                throw new IllegalArgumentException(error);
            }
        }
        if (part == Part.CLOCK_RATE && !(likelihood.getBranchRateModel() instanceof DifferentiableBranchRates)) {
            BranchRateModel branchRateModel = likelihood.getBranchRateModel();
            throw new IllegalArgumentException("mascotGradient part=\"clockRate\" requires a branchRateModel " +
                    "that implements DifferentiableBranchRates, got: " +
                    (branchRateModel == null ? "null (no branchRateModel supplied)" :
                            branchRateModel.getClass().getName()));
        }
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        switch (part) {
            case MIGRATION:
                return likelihood.getMigrationRates();
            case POPULATION_SIZE:
                return likelihood.getPopSizes();
            case CLOCK_RATE:
                return differentiableBranchRates().getRateParameter();
            default:
                throw new IllegalStateException("unknown part: " + part);
        }
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        switch (part) {
            case MIGRATION:
                return likelihood.getMigrationGradientLogDensity();
            case POPULATION_SIZE:
                return likelihood.getPopSizeGradientLogDensity();
            case CLOCK_RATE:
                return clockGradientLogDensity();
            default:
                throw new IllegalStateException("unknown part: " + part);
        }
    }

    private DifferentiableBranchRates differentiableBranchRates() {
        return (DifferentiableBranchRates) likelihood.getBranchRateModel();
    }

    /**
     * Builds the per-branch {@code double[nodeCount-1]} array MascotCore's raw,
     * lineage-id-indexed {@code clockGradient} does not directly have (raw is
     * indexed by node number over the full {@code [0, nodeCount)} range,
     * including the unused root slot; {@code getParameterIndexFromNode} may
     * also remap indices for clock models with a non-identity mapping), then
     * delegates the actual clock-parameter reduction to the branch-rate model
     * itself -- exactly the pattern used by
     * {@code DiscreteTraitBranchRateGradient}/{@code
     * BranchRateGradientForDiscreteTrait} elsewhere in BEAST-X.
     */
    private double[] clockGradientLogDensity() {
        DifferentiableBranchRates branchRates = differentiableBranchRates();
        Tree tree = likelihood.getTreeModel();
        double[] raw = likelihood.getClockGradientLogDensity();
        double[] result = new double[tree.getNodeCount() - 1];
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                result[branchRates.getParameterIndexFromNode(node)] = raw[node.getNumber()];
            }
        }
        return branchRates.updateGradientLogDensity(result, null, 0, result.length);
    }

    @Override
    public double getNumericGradientStepSize() {
        return numericGradientStepSize;
    }

    @Override
    public void setNumericGradientStepSize(double ratio) {
        this.numericGradientStepSize = ratio;
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(
                this,
                Double.NEGATIVE_INFINITY,
                Double.POSITIVE_INFINITY,
                1.0e-3,
                1.0e-8
        );
    }

    @Override
    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }
}
