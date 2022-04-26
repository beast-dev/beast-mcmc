/*
 * CalibratedSpeciationLikelihood.java
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

import dr.evolution.tree.TreeUtils;
import dr.evomodel.tree.TMRCAStatistic;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.distributions.Distribution;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CalibratedSpeciationLikelihood extends AbstractModelLikelihood implements Citable {

    private final SpeciationLikelihood speciationLikelihood;
    private final TreeModel tree;
    private final List<CalibrationLikelihood> calibrationLikelihoods;

    public CalibratedSpeciationLikelihood(String name,
                                          SpeciationLikelihood speciationLikelihood,
                                          TreeModel tree,
                                          List<CalibrationLikelihood> calibrationLikelihoods) {
        super(name);
        this.speciationLikelihood = speciationLikelihood;
        this.tree = tree;
        this.calibrationLikelihoods = calibrationLikelihoods;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        double lnL = speciationLikelihood.getLogLikelihood();
        for (CalibrationLikelihood calibrationLikelihood : calibrationLikelihoods) {
            lnL += calibrationLikelihood.getLogLikelihood();
        }
        return lnL;
    }

    public List<CalibrationLikelihood> getCalibrationLikelihoods() {
        return calibrationLikelihoods;
    }

    @Override
    public void makeDirty() {
        speciationLikelihood.makeDirty();
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Calibration prior distributions";
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> calibrationCitations = new ArrayList<>();
        for (CalibrationLikelihood calibrationLikelihood : calibrationLikelihoods) {
            if (calibrationLikelihood.getDistribution() instanceof Citable) {
                Citation citation = ((Citable) calibrationLikelihood.getDistribution()).getCitations().get(0);
                if (!calibrationCitations.contains(citation)) {
                    calibrationCitations.add(((Citable) calibrationLikelihood.getDistribution()).getCitations().get(0));
                }
            }
        }
        return calibrationCitations;
    }

    public static class CalibrationLikelihood {

        private final TMRCAStatistic tmrcaStatistic;
        private final Distribution distribution;
        private final int mrcaNodeNumber;

        public CalibrationLikelihood(TMRCAStatistic tmrcaStatistic,
                                     Distribution distribution) {
            this.tmrcaStatistic = tmrcaStatistic;
            this.distribution = distribution;
            this.mrcaNodeNumber = TreeUtils.getCommonAncestorNode(tmrcaStatistic.getTree(), tmrcaStatistic.getLeafSet()).getNumber();
        }

        public int getMrcaNodeNumber() {
            return mrcaNodeNumber;
        }

        public Distribution getDistribution() {
            return distribution;
        }

        public double getLogLikelihood() {
            final double nodeHeight = getNodeHeight();
            return distribution.logPdf(nodeHeight);
        }

        private final double getNodeHeight() {
            return tmrcaStatistic.getTree().getNodeHeight(tmrcaStatistic.getTree().getNode(mrcaNodeNumber));
        }

        public double[] getGradientLogDensity() {
            return ((GradientProvider) distribution).getGradientLogDensity(getNodeHeight());
        }
    }
}
