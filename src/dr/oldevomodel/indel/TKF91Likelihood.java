/*
 * TKF91Likelihood.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.oldevomodel.indel;

import dr.evolution.alignment.Alignment;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodelxml.indel.TKF91LikelihoodParser;
import dr.inference.model.*;


/**
 * Calculates the likelihood of a set of continuous attributes on a tree.
 *
 * @author Alexei Drummond
 * @version $Id: TKF91Likelihood.java,v 1.17 2005/06/20 07:24:25 rambaut Exp $
 */
public class TKF91Likelihood extends AbstractModelLikelihood {

    public TKF91Likelihood(TreeModel treeModel, Alignment alignment, GammaSiteModel siteModel, TKF91Model tkfModel) {

        super(TKF91LikelihoodParser.TKF91_LIKELIHOOD);

        if (siteModel.getAlphaParameter() != null)
            throw new IllegalArgumentException("TKF91 model cannot handle gamma-distributed rates");
        if (siteModel.getPInvParameter() != null)
            throw new IllegalArgumentException("TKF91 model cannot handle invariant sites");

        addModel(siteModel);
        addModel(tkfModel);
        addModel(treeModel);

        this.treeModel = treeModel;
        this.alignment = alignment;
        this.siteModel = siteModel;
        this.tkfModel = tkfModel;

        recursion = new HomologyRecursion();
        recursion.init(
                treeModel,
                alignment,
                siteModel.getSubstitutionModel(),
                siteModel.getMutationRateParameter().getParameterValue(0),
                tkfModel.getLengthDistributionValue(),
                tkfModel.getDeathRate(1));

        addStatistic(new AlignmentLengthStatistic());
    }

    public class AlignmentLengthStatistic extends Statistic.Abstract {

        public double getStatisticValue(int dim) {
            return alignment.getSiteCount();
        }

        public int getDimension() {
            return 1;
        }

        public String getStatisticName() {
            return "alignmentLength";
        }
    }

    public void acceptState() {
        //throw new RuntimeException("Not implemented!");
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        //throw new RuntimeException("Not implemented!");
    }

    public void handleModelChangedEvent(Model m, Object o, int i) {
        //throw new RuntimeException("Not implemented!");
    }

    /**
     * Sets the alignment of this tkf91likelihood
     */
    public void setAlignment(Alignment a) {
        alignment = a;
        //System.out.println("Set new alignment");
    }

    protected void storeState() {
        //super.storeState();
        //storedAlignment = alignment;
        //System.out.println("Stored alignment");

    }

    protected void restoreState() {
        //super.restoreState();
        //alignment = storedAlignment;
        //System.out.println("restored alignment");
    }

    public Alignment getAlignment() {
        return alignment;
    }

    public GammaSiteModel getSiteModel() {
        return siteModel;
    }

    public TreeModel getTreeModel() {
        return treeModel;
    }

    public void makeDirty() { // this is always dirty
    }

    public Model getModel() {
        return this;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    public double getLogLikelihood() {

        recursion.init(
                treeModel,
                alignment,
                siteModel.getSubstitutionModel(),
                siteModel.getMutationRateParameter().getParameterValue(0),
                tkfModel.getLengthDistributionValue(),
                tkfModel.getDeathRate(1));

        double logL = recursion.recursion();

        //System.out.println("logL = " + logL);
        return logL;
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }


    public String toString() {
        if (id != null) {
            return id;
        }
        return super.toString();
    }

    protected boolean getLikelihoodKnown() {
        return false;
    }

    private final TreeModel treeModel;
    private Alignment alignment;
    //private Alignment storedAlignment;
    private final GammaSiteModel siteModel;
    private final TKF91Model tkfModel;
    private HomologyRecursion recursion = null;
}