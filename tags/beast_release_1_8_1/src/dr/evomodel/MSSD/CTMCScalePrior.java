/*
 * CTMCScalePrior.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.MSSD;

import dr.app.beagle.evomodel.substmodel.SubstitutionModel;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.GammaFunction;

/**
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * @author Marc A. Suchard
 *         <p/>
 *         Date: Aug 22, 2008
 *         Time: 3:26:57 PM
 */
public class CTMCScalePrior extends AbstractModelLikelihood {
    final private Parameter ctmcScale;
    final private TreeModel treeModel;
    private double treeLength;
    private boolean treeLengthKnown;

    final private boolean reciprocal;
    final private SubstitutionModel substitutionModel;
    final private boolean trial;

    private static final double logGammaOneHalf = GammaFunction.lnGamma(0.5);

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel) {
        this(name, ctmcScale, treeModel, false);
    }

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel, boolean reciprocal) {
        this(name, ctmcScale, treeModel, reciprocal, null);
    }

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel, boolean reciprocal,
                          SubstitutionModel substitutionModel) {
        this(name, ctmcScale, treeModel, reciprocal, substitutionModel, false);
    }

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel, boolean reciprocal,
                          SubstitutionModel substitutionModel, boolean trial) {
        super(name);
        this.ctmcScale = ctmcScale;
        this.treeModel = treeModel;
        addModel(treeModel);
        treeLengthKnown = false;
        this.reciprocal = reciprocal;
        this.substitutionModel = substitutionModel;
        this.trial = trial;
    }

    private void updateTreeLength() {
        treeLength = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == treeModel) {
            treeLengthKnown = false;
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    }

    protected void storeState() {
    }

    protected void restoreState() {
        treeLengthKnown = false;
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    private double calculateTrialLikelihood() {
        double totalTreeTime = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());

        double[] eigenValues = substitutionModel.getEigenDecomposition().getEigenValues();
        // Find second largest
        double lambda2 = Double.NEGATIVE_INFINITY;
        for (double l : eigenValues) {
            if (l > lambda2 && l < 0.0) {
                lambda2 = l;
            }
        }
        lambda2 = -lambda2;

        double logNormalization = 0.5 * Math.log(lambda2) - logGammaOneHalf;

        double logLike = 0;
        for (int i = 0; i < ctmcScale.getDimension(); ++i) {
            double ab = ctmcScale.getParameterValue(i) * totalTreeTime;
            logLike += logNormalization - 0.5 * Math.log(ab) - ab * lambda2;
        }
        return logLike;
    }

    public double getLogLikelihood() {

//        if (!treeLengthKnown) {
//            updateTreeLength();
//            treeLengthKnown = true;
//        }
//        double totalTreeTime = treeLength;

        if (trial) return calculateTrialLikelihood();

        double totalTreeTime = Tree.Utils.getTreeLength(treeModel, treeModel.getRoot());
        if (reciprocal) {
            totalTreeTime = 1.0 / totalTreeTime;
        }
        if (substitutionModel != null) {
            double[] eigenValues = substitutionModel.getEigenDecomposition().getEigenValues();
            // Find second largest
            double lambda2 = Double.NEGATIVE_INFINITY;
            for (double l : eigenValues) {
                if (l > lambda2 && l < 0.0) {
                    lambda2 = l;
                }
            }
            totalTreeTime *= -lambda2; // TODO Should this be /=?
        }
        double logNormalization = 0.5 * Math.log(totalTreeTime) - logGammaOneHalf;
        double logLike = 0;
        for (int i = 0; i < ctmcScale.getDimension(); ++i) {
            double ab = ctmcScale.getParameterValue(i);
            logLike += logNormalization - 0.5 * Math.log(ab) - ab * totalTreeTime; // TODO Change to treeLength and confirm results
        }
        return logLike;
    }

    public void makeDirty() {
        treeLengthKnown = false;
    }
}
