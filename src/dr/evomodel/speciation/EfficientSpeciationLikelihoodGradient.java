/*
 * EfficientSpeciationLikelihoodGradient.java
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

import dr.evolution.tree.TreeTrait;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.xml.Reportable;

/**
 * @author Andy Magee
 * @author Yucai Shao
 * @author Marc Suchard
 */
public class EfficientSpeciationLikelihoodGradient extends AbstractModel
        implements GradientWrtParameterProvider, Reportable, Loggable {

    static final String GRADIENT_KEY = "speciationGradient";

    private final EfficientSpeciationLikelihood likelihood;
    private final Parameter parameter;
    private final SpeciationLikelihoodGradient.WrtParameter wrtParameter;
    private final TreeModel tree;
    private final SpeciationModel speciationModel;
    private final BigFastTreeIntervals treeIntervals;
    private final SpeciationModelGradientProvider provider;

    private final TreeTrait gradientProvider;

    private boolean gradientKnown;
    private double[] gradient;
    private double[] storedGradient;

    public EfficientSpeciationLikelihoodGradient(EfficientSpeciationLikelihood likelihood,
                                                 SpeciationLikelihoodGradient.WrtParameter wrtParameter) {
        super("efficientSpeciationLikelihoodGradient");

        this.likelihood = likelihood;
        this.wrtParameter = wrtParameter;

        this.tree = likelihood.getTreeModel();

        this.speciationModel = likelihood.getSpeciationModel();
        this.treeIntervals = likelihood.getTreeIntervals();
        this.provider = likelihood.getGradientProvider();
        this.parameter = wrtParameter.getParameter(provider, tree);

        likelihood.addModel(this);

        if (wrtParameter == SpeciationLikelihoodGradient.WrtParameter.NODE_HEIGHT) {
            speciationModel.addModelListener(this);
            treeIntervals.addModelListener(this);
        }

        gradientKnown = false;

        this.gradientProvider = getGradientDelegateSingleton(likelihood);
    }

    private TreeTrait getGradientDelegateSingleton(EfficientSpeciationLikelihood likelihood) {
        TreeTrait singleton = likelihood.getTreeTrait(GRADIENT_KEY);
        if (singleton == null) {
            CachedGradientDelegate delegate = new CachedGradientDelegate(likelihood);
            addModel(delegate);
            singleton = delegate;
            likelihood.addTrait(singleton);
        }
        return singleton;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
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
    public double[] getGradientLogDensity() {
        if (wrtParameter == SpeciationLikelihoodGradient.WrtParameter.NODE_HEIGHT) {
            if (!gradientKnown) {
                gradient = wrtParameter.getGradientLogDensity(provider, tree); // TODO Harmonize
                gradientKnown = true;
            }
            return gradient;
        } else {
            return wrtParameter.filter((double[]) gradientProvider.getTrait(null, null));
        }
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "SpeciationLikelihoodGradient check");
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0.0, Double.POSITIVE_INFINITY, 1E-3);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == speciationModel || model == treeIntervals) {
            gradientKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown model: " + model.getId());
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new IllegalArgumentException("Unknown variable: " + variable.getId());
    }

    @Override
    protected void storeState() {
        if (wrtParameter == SpeciationLikelihoodGradient.WrtParameter.NODE_HEIGHT) {
            if (gradient != null) {
                if (storedGradient == null) {
                    storedGradient = new double[gradient.length];
                }
                System.arraycopy(gradient, 0, storedGradient, 0, gradient.length);
            }
        }
    }

    @Override
    protected void restoreState() {
        if (wrtParameter == SpeciationLikelihoodGradient.WrtParameter.NODE_HEIGHT) {
            double[] tmp = gradient;
            gradient = storedGradient;
            storedGradient = tmp;
        }
    }

    @Override
    protected void acceptState() {
        // Do nothing
    }
}
