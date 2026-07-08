/*
 * MascotGradient.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * HMC-facing gradient provider for the flat log-MASCOT parameter vector.
 */
public final class MascotGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    public static final String MASCOT_GRADIENT = "mascotGradient";

    private final MascotLikelihood likelihood;
    private double numericGradientStepSize = GradientWrtParameterProvider.super.getNumericGradientStepSize();

    public MascotGradient(MascotLikelihood likelihood) {
        this.likelihood = likelihood;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return likelihood.getTheta();
    }

    @Override
    public int getDimension() {
        return getParameter().getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return likelihood.getGradientLogDensity();
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
