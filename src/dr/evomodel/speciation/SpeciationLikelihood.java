/*
 * SpeciationLikelihood.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

import dr.evolution.tree.Tree;
import dr.evolution.util.Units;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * A likelihood function for speciation processes. Takes a tree and a speciation model.
 *
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @version $Id: SpeciationLikelihood.java,v 1.10 2005/05/18 09:51:11 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class SpeciationLikelihood extends AbstractModel implements Likelihood, Units {

    // PUBLIC STUFF

    public static final String SPECIATION_LIKELIHOOD = "speciationLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "speciesTree";

    public SpeciationLikelihood(Tree tree, SpeciationModel speciationModel) {
        this(SPECIATION_LIKELIHOOD, tree, speciationModel);
    }

    public SpeciationLikelihood(String name, Tree tree, SpeciationModel speciationModel) {

        super(name);

        this.tree = tree;
        this.speciationModel = speciationModel;
        if (tree instanceof TreeModel) {
            addModel((TreeModel)tree);
        }
        if (speciationModel != null) {
            addModel(speciationModel);
        }
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    // **************************************************************
    // ParameterListener IMPLEMENTATION
    // **************************************************************

    protected final void handleParameterChangedEvent(Parameter parameter, int index) { } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: in this case the intervals
     */
    protected final void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: that is the intervals of the tree.
     */
    protected final void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() { } // nothing to do

    /**
     * Adopt the state of the model from source.
     * @param source the source model
     */
    protected final void adoptState(Model source) {
        // all we need to do is force a recalculation of intervals
        makeDirty();
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() { return this; }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public final void makeDirty() {
        likelihoodKnown = false;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     * @return the log likelihood
     */
    public double calculateLogLikelihood() {

        double logL = 0.0;

        for (int j = 0; j < tree.getInternalNodeCount(); j++) {
            logL += speciationModel.logNodeProbability( tree, tree.getInternalNode(j));
        }

        if (speciationModel.includeExternalNodesInLikelihoodCalculation()) {
            for (int j = 0; j < tree.getExternalNodeCount(); j++) {
                logL += speciationModel.logNodeProbability( tree, tree.getExternalNode(j));
            }
        }

        return logL;
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {
        return new dr.inference.loggers.LogColumn[] {
                new LikelihoodColumn(getId())
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) { super(label); }
        public double getDoubleValue() { return getLogLikelihood(); }
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u)
    {
        speciationModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits()
    {
        return speciationModel.getUnits();
    }

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return SPECIATION_LIKELIHOOD; }

        public Object parseXMLObject(XMLObject xo) {

            XMLObject cxo = (XMLObject)xo.getChild(MODEL);
            SpeciationModel specModel = (SpeciationModel)cxo.getChild(SpeciationModel.class);

            cxo = (XMLObject)xo.getChild(TREE);
            TreeModel treeModel = (TreeModel)cxo.getChild(TreeModel.class);

            return new SpeciationLikelihood(treeModel, specModel);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the speciation.";
        }

        public Class getReturnType() { return Likelihood.class; }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(MODEL, new XMLSyntaxRule[] {
                        new ElementRule(SpeciationModel.class)
                }),
                new ElementRule(TREE, new XMLSyntaxRule[] {
                        new ElementRule(TreeModel.class)
                }),
        };
    };



    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /** The speciation model. */
    SpeciationModel speciationModel = null;

    /** The tree. */
    Tree tree = null;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}