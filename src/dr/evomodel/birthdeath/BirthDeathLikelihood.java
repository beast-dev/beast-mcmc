/*
 * SpeciationLikelihood.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

package dr.evomodel.birthdeath;

import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.Reportable;

import java.util.Set;

/**
 * A likelihood function for birth-death processes. Takes a tree and a birth-death model.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class BirthDeathLikelihood extends AbstractModelLikelihood implements Reportable, Units {

    // PUBLIC STUFF
    /**
     * @param tree            the tree
     * @param birthDeathModel the model of birth-death
     * @param id              a unique identifier for this likelihood
     * @param exclude         taxa to exclude from this model
     */
    public BirthDeathLikelihood(Tree tree, BirthDeathModel birthDeathModel, Set<Taxon> exclude, String id) {
        this("birthDeathLikelihood", tree, birthDeathModel, exclude);
        setId(id);
    }

    public BirthDeathLikelihood(Tree tree, BirthDeathModel birthDeathModel, String id) {
        this(tree, birthDeathModel, null, id);
    }

    public BirthDeathLikelihood(String name, Tree tree, BirthDeathModel birthDeathModel, Set<Taxon> exclude) {

        super(name);

        this.tree = tree;
        this.birthDeathModel = birthDeathModel;
        this.exclude = exclude;

        if (tree instanceof Model) {
            addModel((Model) tree);
        }
        if (birthDeathModel != null) {
            addModel(birthDeathModel);
        }
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    // TODO Make final again after done with EfficientSpeciationLikelihood
    protected void handleModelChangedEvent(Model model, Object object, int index) {
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
    double calculateLogLikelihood() {

        if (exclude != null) {
            return birthDeathModel.calculateTreeLogLikelihood(tree, exclude);
        }

        return birthDeathModel.calculateTreeLogLikelihood(tree);
    }

    // Super-clean interface (just one intrusive function) and a better place, since `Likelihood`s have gradients (`Model`s do not).
    public BirthDeathModelGradientProvider getGradientProvider() {
        if (gradientProvider == null) {
            gradientProvider = birthDeathModel.getProvider();
        }
        return gradientProvider;
    }

    private BirthDeathModelGradientProvider gradientProvider = null;

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

    @Override
    public String getReport() {
        String report = "BirthDeathLikelihoodReport:\n" +
                "Model: " + birthDeathModel.getModelName() + "\n" +
                "ID: " + birthDeathModel.getId() + "\n" +
                "lnL: " + calculateLogLikelihood();
        return report;
    }

    private final class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
        public LikelihoodColumn(String label) {
            super(label);
        }

        public double getDoubleValue() {
            return getLogLikelihood();
        }
    }

    public BirthDeathModel getBirthDeathModel() {
        return birthDeathModel;
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are
     * measured in.
     */
    public final void setUnits(Type u) {
        birthDeathModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are
     * measured in.
     */
    public final Type getUnits() {
        return birthDeathModel.getUnits();
    }

    @Override
    public String prettyName() {
        String s = birthDeathModel.getClass().getName();
        String[] parts = s.split("\\.");
        s = parts[parts.length - 1];
        if( birthDeathModel.getId() != null ) {
           s = s + '/' + birthDeathModel.getId();
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
    final BirthDeathModel birthDeathModel;

    /**
     * The tree.
     */
    Tree tree;
    private final Set<Taxon> exclude;

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}