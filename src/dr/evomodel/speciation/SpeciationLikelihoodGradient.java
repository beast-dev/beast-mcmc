/*
 * UltrametricSpeciationGradient.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
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
public class SpeciationLikelihoodGradient implements GradientWrtParameterProvider, Reportable, Loggable {

    private final SpeciationModel speciationModel;
    private final SpeciationLikelihood likelihood;
    private final Parameter nodeHeightParameter;
    private final TreeModel tree;

    public SpeciationLikelihoodGradient(SpeciationLikelihood likelihood,
                                        TreeModel tree) {

        this.likelihood = likelihood;
        this.speciationModel = likelihood.speciationModel;
        this.tree = tree;
        this.nodeHeightParameter = new NodeHeightProxyParameter("internalNodeParameter", tree, true);

    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return nodeHeightParameter;
    }

    @Override
    public int getDimension() {
        return nodeHeightParameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradient = new double[tree.getInternalNodeCount()];

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            gradient[i] = speciationModel.getNodeGradient(tree, tree.getNode(i + tree.getExternalNodeCount()));
        }

        return gradient;
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "SpeciationLikelihoodGradient check");
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-3);
    }
}
