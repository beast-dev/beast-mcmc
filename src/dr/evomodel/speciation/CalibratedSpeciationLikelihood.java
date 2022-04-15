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
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.distributions.Distribution;

import java.util.List;

/**
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class CalibratedSpeciationLikelihood extends AbstractModelLikelihood {

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

    @Override
    public void makeDirty() {

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

        public double getLogLikelihood() {
            final double nodeHeight = tmrcaStatistic.getTree().getNodeHeight(tmrcaStatistic.getTree().getNode(mrcaNodeNumber));
            return distribution.logPdf(nodeHeight);
        }

    }

}
