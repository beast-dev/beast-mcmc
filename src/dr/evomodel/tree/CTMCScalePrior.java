/*
 * CTMCScalePrior.java
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

package dr.evomodel.tree;

import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.math.GammaFunction;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 * @author Marc A. Suchard
 *         <p/>
 *         Date: Aug 22, 2008
 *         Time: 3:26:57 PM
 */
public class CTMCScalePrior extends AbstractModelLikelihood implements GradientWrtParameterProvider, Citable {
    final private Parameter ctmcScale;
    final private TreeModel treeModel;
    private Set<Taxon> taxa = null;
    private double treeLength;
    private double storedTreeLength;
    private boolean treeLengthKnown;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown;

    final private boolean reciprocal;
    final private SubstitutionModel substitutionModel;
    final private boolean trial;

    private static final double shape = 0.5;

    private static final double logGammaOneHalf = GammaFunction.lnGamma(shape);

    public CTMCScalePrior(String name, Parameter ctmcScale, TreeModel treeModel, TaxonList taxonList, boolean reciprocal,
                          SubstitutionModel substitutionModel, boolean trial) {
        super(name);
        this.ctmcScale = ctmcScale;
        this.treeModel = treeModel;

        addVariable(ctmcScale);

        if (taxonList != null) {
            this.taxa = new HashSet<>();
            for (Taxon taxon : taxonList) {
                this.taxa.add(taxon);
            }
        }
        addModel(treeModel);

        treeLengthKnown = false;
        likelihoodKnown = false;

        this.reciprocal = reciprocal;
        this.substitutionModel = substitutionModel;
        if (substitutionModel != null) {
            addModel(substitutionModel);
        }
        this.trial = trial;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        treeLengthKnown = false;
        likelihoodKnown = false;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        treeLengthKnown = false;
        likelihoodKnown = false;
    }

    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedTreeLength = treeLength;
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        treeLength = storedTreeLength;
        treeLengthKnown = true;
        likelihoodKnown = true;
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    private double calculateTrialLikelihood() {
        double totalTreeTime = getTreeLength();

        double lambda2 = -secondLargestEigenvalue();

        double logNormalization = shape * Math.log(lambda2) - logGammaOneHalf;

        double logLike = 0;
        for (int i = 0; i < ctmcScale.getDimension(); ++i) {
            double ab = ctmcScale.getParameterValue(i) * totalTreeTime;
            logLike += logNormalization - shape * Math.log(ab) - ab * lambda2;
        }
        return logLike;
    }

    private double secondLargestEigenvalue() {
        double[] eigenValues = substitutionModel.getEigenDecomposition().getEigenValues();
        // Find second largest
        double lambda2 = Double.NEGATIVE_INFINITY;
        for (double l : eigenValues) {
            if (l > lambda2 && l < 0.0) {
                lambda2 = l;
            }
        }
        return lambda2;
    }

    private double scaledTotalTreeTime() {
        double totalTreeTime = getTreeLength();
        if (reciprocal) {
            totalTreeTime = 1.0 / totalTreeTime;
        }

        if (substitutionModel != null) {
            totalTreeTime *= -secondLargestEigenvalue(); // TODO Should this be /=?
        }
        return totalTreeTime;
    }
    
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }

        return logLikelihood;
    }


    private double calculateLogLikelihood() {

        if (trial) return calculateTrialLikelihood();

        double totalTreeTime = scaledTotalTreeTime();

        double logNormalization = shape * Math.log(totalTreeTime) - logGammaOneHalf;
        double logLike = 0;
        for (int i = 0; i < ctmcScale.getDimension(); ++i) {
            double ab = ctmcScale.getParameterValue(i);
            logLike += logNormalization - shape * Math.log(ab) - ab * totalTreeTime; // TODO Change to treeLength and confirm results
        }
        return logLike;
    }

    private double getTreeLength() {
        if (!treeLengthKnown) {
            if (taxa == null) {
                treeLength = TreeUtils.getTreeLength(treeModel, treeModel.getRoot());
            } else {
                treeLength = TreeUtils.getSubTreeLength(treeModel, taxa);
            }
            treeLengthKnown = true;
        }

        return treeLength;
    }

    public void makeDirty() {
        treeLengthKnown = false;
        likelihoodKnown = false;
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "CTMC Scale Reference Prior model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("MAR", "Ferreira"),
                    new Author("MA", "Suchard")
            },
            "Bayesian analysis of elapsed times in continuous-time Markov chains",
            2008,
            "Canadian Journal of Statistics",
            36,
            355, 368,
            Citation.Status.PUBLISHED
    );

    @Override
    public Likelihood getLikelihood() {
        return this;
    }

    @Override
    public Parameter getParameter() {
        return ctmcScale;
    }

    @Override
    public int getDimension() {
        return ctmcScale.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradLogLike = new double[ctmcScale.getDimension()];

        double totalTreeTime = scaledTotalTreeTime();

        for (int i = 0; i < ctmcScale.getDimension(); ++i) {
            double ab = ctmcScale.getParameterValue(i);
            gradLogLike[i] = -shape / ab - totalTreeTime;
        }
        return gradLogLike;
    }
}
