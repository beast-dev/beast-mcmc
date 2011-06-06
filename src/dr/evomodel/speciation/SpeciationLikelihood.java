/*
 * SpeciationLikelihood.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.speciation;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.evomodelxml.speciation.SpeciationLikelihoodParser;
import dr.inference.model.*;
import dr.math.distributions.Distribution;

import java.util.List;
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

    public SpeciationLikelihood(Tree tree, SpeciationModel speciationModel, String id,
                                List<Distribution> dists, List<Taxa> taxa, List<Boolean> forParent,
                                Statistic userPDF) {
        this(tree, speciationModel, id);
        this.distribution = dists.toArray(new Distribution[dists.size()]);
        this.taxa = new int[taxa.size()][];
        this.forParent = new boolean[taxa.size()];

        for(int k = 0; k < taxa.size(); ++k) {
            final Taxa tk = taxa.get(k);
            final int tkcount = tk.getTaxonCount();
            this.taxa[k] = new int[tkcount];
            for(int nt = 0; nt < tkcount; ++nt) {
                this.taxa[k][nt] = tree.getTaxonIndex(tk.getTaxon(nt));
            }
            this.forParent[k] = forParent.get(k);
        }
        this.calibrationLogPDF = userPDF;

        if( userPDF != null && taxa.size() == 2 ) {
            assert this.taxa[0].length < this.taxa[1].length;
            for( int t : this.taxa[0] ) {
                boolean found = false;
                for( int x : this.taxa[1] ) {
                    if( x == t ) {
                        found = true;
                        break;
                    }
                }
                assert found;
            }
        }
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

        if ( distribution != null ) {
            //return speciationModel.calculateTreeLogLikelihood(tree, taxa, distribution, coefficients);
            return speciationModel.calculateTreeLogLikelihood(tree, taxa, forParent, distribution, calibrationLogPDF);
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

    private  Distribution[] distribution = null;
    private  int[][] taxa = null;
    private  boolean[] forParent = null;

    private  Statistic calibrationLogPDF = null;
    //private double[] coefficients = null;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}