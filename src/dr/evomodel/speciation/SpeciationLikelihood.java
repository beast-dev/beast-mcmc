/*
 * SpeciationLikelihood.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Set;

/**
 * A likelihood function for speciation processes. Takes a tree and a speciation model.
 * <p/>
 * Parts of this class were derived from C++ code provided by Oliver Pybus.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: SpeciationLikelihood.java,v 1.10 2005/05/18 09:51:11 rambaut Exp $
 */
public class SpeciationLikelihood extends AbstractModelLikelihood implements Units {

    // PUBLIC STUFF
    /**
     * @param tree            the tree
     * @param speciationModel the model of speciation
     * @param id              a unique identifier for this likelihood
     * @param exclude         taxa to exclude from this model
     */
    public SpeciationLikelihood(Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude, String id) {
        this(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, tree, speciationModel, exclude);
        setId(id);
    }

    public SpeciationLikelihood(Tree tree, SpeciationModel speciationModel, String id) {
        this(tree, speciationModel, null, id);
    }

    public SpeciationLikelihood(String name, Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude) {

        super(name);

        this.tree = tree;
        this.speciationModel = speciationModel;
        this.exclude = exclude;

        if (tree instanceof Model) {
            addModel((Model) tree);
        }
        if (speciationModel != null) {
            addModel(speciationModel);
        }
    }

    public SpeciationLikelihood(Tree tree, SpeciationModel specModel, String id, CalibrationPoints calib) {
        this(tree, specModel, id);
        this.calibration = calib;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected final void handleModelChangedEvent(Model model, Object object, int index) {
        likelihoodKnown = false;
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
    } // No parameters to respond to

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    /**
     * Stores the precalculated state: likelihood
     */
    protected final void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    /**
     * Restores the precalculated state: computed likelihood
     */
    protected final void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    protected final void acceptState() {
    } // nothing to do

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

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
     *
     * @return the log likelihood
     */
    private double calculateLogLikelihood() {
        if (exclude != null) {
            return speciationModel.calculateTreeLogLikelihood(tree, exclude);
        }

        if ( calibration != null ) {
            return speciationModel.calculateTreeLogLikelihood(tree, calibration);
        }

        return speciationModel.calculateTreeLogLikelihood(tree);
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    /**
     * @return the log columns.
     */
    public final dr.inference.loggers.LogColumn[] getColumns() {

        String columnName = getId();
        if (columnName == null) columnName = getModelName() + ".likelihood";

        return new dr.inference.loggers.LogColumn[]{
                new LikelihoodColumn(columnName)
        };
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u) {
        speciationModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits() {
        return speciationModel.getUnits();
    }

    @Override
    public String prettyName() {
        String s = speciationModel.getClass().getName();
        String[] parts = s.split("\\.");
        s = parts[parts.length - 1];
        if( speciationModel.getId() != null ) {
           s = s + '/' + speciationModel.getId();
        }
        s = s + '(' + tree.getId() + ')';
        return s;
    }
    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    /**
     * The speciation model.
     */
    SpeciationModel speciationModel = null;

    /**
     * The tree.
     */
    Tree tree = null;
    private final Set<Taxon> exclude;

    private CalibrationPoints calibration;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}