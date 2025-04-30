/*
 * TaxonEffectGradient.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.TaxonEffectTraitDataModel;
import dr.evomodel.treedatalikelihood.preorder.TipGradientViaFullConditionalDelegate;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.util.StopWatch;
import dr.util.TaskPool;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class TaxonEffectGradient implements GradientWrtParameterProvider, Reportable {

    private final List<Partition> partitions;
    private Likelihood compoundLikelihood;
    private final Parameter effects;
    private final int dimTrait;
    private final int nTraits;
    private final int nTaxa;

    @SuppressWarnings("unused")
    public TaxonEffectGradient(List<Partition> partitions,
                               TaskPool taskPool,
                               ThreadUseProvider threadUseProvider) {

        this.partitions = partitions;
        this.effects = partitions.get(0).model.getEffects();

        this.dimTrait = partitions.get(0).delegate.getTraitDim();
        this.nTraits = partitions.get(0).delegate.getTraitCount();
        this.nTaxa = partitions.get(0).tree.getExternalNodeCount();

        if (TIMING) {
            int length = 5;
            stopWatches = new StopWatch[length];
            for (int i = 0; i < length; ++i) {
                stopWatches[i] = new StopWatch();
            }
        }
    }

    @Override
    public Likelihood getLikelihood() {
        if (compoundLikelihood == null) {

            List<Likelihood> likelihoods = new ArrayList<>();
            for (Partition p : partitions) {
                likelihoods.add(p.likelihood);
            }

            compoundLikelihood = new CompoundLikelihood(likelihoods);
        }
        return compoundLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return effects;
    }

    @Override
    public int getDimension() {
        return nTaxa * dimTrait * nTraits;
    }

    @Override
    public double[] getGradientLogDensity() {

        if (TIMING) {
            stopWatches[0].start();
        }

        double[] gradient = new double[getDimension()];

        for (Partition p : partitions) {
            for (int taxon = 0; taxon < nTaxa; ++taxon) {
                Tree tree = p.tree;
                TreeTrait trait = p.treeTraitProvider;
                double[] taxonGradient = (double[]) trait.getTrait(tree, tree.getExternalNode(taxon));
                final int effectIndex = p.model.getMap().getEffectIndex(taxon);
                final int sign = p.model.getMap().getSign(taxon);
                final int offsetOutput = effectIndex * dimTrait * nTraits;
                for (int i = 0; i < dimTrait; ++i) {
                    gradient[offsetOutput + i] -= sign * taxonGradient[i];
                }
            }
        }

        if (TIMING) {
            stopWatches[0].stop();
        }

        return gradient;
    }

    public enum ThreadUseProvider {
        PARALLEL {
            @Override
            boolean usePool() {
                return true;
            }
        },

        SERIAL {
            @Override
            boolean usePool() {
                return false;
            }
        };

        @SuppressWarnings("unused")
        abstract boolean usePool();
    }

    @Override
    public String getReport() {

        String report = "";

        if (TIMING) {
            report += timingInfo();
        }

        report += GradientWrtParameterProvider.getReportAndCheckForError(
                this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY,
                1E-3);

        if (TIMING) {
            report += timingInfo();
        }

        return report;
    }

    private String timingInfo() {
        StringBuilder sb = new StringBuilder("\nTiming in TaxonEffectGradient\n");
        for (StopWatch stopWatch : stopWatches) {
            sb.append("\t").append(stopWatch.toString()).append("\n");
            stopWatch.reset();
        }
        return sb.toString();
    }

    protected StopWatch[] stopWatches;
    protected static final boolean TIMING = true;

    public static class Partition {
        public final Tree tree;
        public final TreeDataLikelihood likelihood;
        public final ContinuousDataLikelihoodDelegate delegate;
        public final TaxonEffectTraitDataModel model;
        public final TreeTrait treeTraitProvider;

        public Partition(TreeDataLikelihood likelihood,
                         ContinuousDataLikelihoodDelegate delegate,
                         TaxonEffectTraitDataModel model) {
            this.tree = likelihood.getTree();
            this.likelihood = likelihood;
            this.delegate = delegate;
            this.model = model;

            String traitName = delegate.getModelName();
            int dimTrait = delegate.getTraitDim();

            String name = TipGradientViaFullConditionalDelegate.getName(traitName);
            TreeTrait test = likelihood.getTreeTrait(name);

            if (test == null) {
                delegate.addFullConditionalGradientTrait(traitName, 0, dimTrait);
            }

            this.treeTraitProvider = likelihood.getTreeTrait(name);
        }
    }
}
