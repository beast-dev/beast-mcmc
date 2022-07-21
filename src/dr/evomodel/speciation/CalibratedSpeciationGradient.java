/*
 * CalibratedSpeciationGradient.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.speciation;

import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CalibratedSpeciationGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final SpeciationLikelihoodGradient speciationLikelihoodGradient;
    private final SpeciationLikelihood unCalibratedSpeciationLikelihood;
    private final CalibratedSpeciationLikelihood calibratedSpeciationLikelihood;
    private final TreeModel tree;

    public CalibratedSpeciationGradient(SpeciationLikelihoodGradient speciationLikelihoodGradient,
                                        CalibratedSpeciationLikelihood calibratedSpeciationLikelihood) {
        this.speciationLikelihoodGradient = speciationLikelihoodGradient;
        this.unCalibratedSpeciationLikelihood = (SpeciationLikelihood) speciationLikelihoodGradient.getLikelihood();
        this.calibratedSpeciationLikelihood = calibratedSpeciationLikelihood;
        this.tree = speciationLikelihoodGradient.getTree();
    }

    @Override
    public Likelihood getLikelihood() {
        return calibratedSpeciationLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return speciationLikelihoodGradient.getParameter();
    }

    @Override
    public int getDimension() {
        return speciationLikelihoodGradient.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        double[] gradient = speciationLikelihoodGradient.getGradientLogDensity();
        for (CalibratedSpeciationLikelihood.CalibrationLikelihood calibrationLikelihood : calibratedSpeciationLikelihood.getCalibrationLikelihoods()) {
            final int nodeIndex = calibrationLikelihood.getMrcaNodeNumber() - tree.getExternalNodeCount();
            double[] calibrationGradient = calibrationLikelihood.getGradientLogDensity();
            assert(calibrationGradient.length == 1);
            gradient[nodeIndex] += calibrationGradient[0];
        }
        return gradient;
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "CalibratedSpeciationGradient report");
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-2);
    }
}
