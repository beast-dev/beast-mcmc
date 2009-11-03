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
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.inference.model.*;
import dr.xml.*;

import java.util.HashSet;
import java.util.Set;
import java.util.logging.Logger;

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

    public static final String SPECIATION_LIKELIHOOD = "speciationLikelihood";
    public static final String MODEL = "model";
    public static final String TREE = "speciesTree";
    public static final String INCLUDE = "include";
    public static final String EXCLUDE = "exclude";

    /**
     * @param tree            the tree
     * @param speciationModel the model of speciation
     * @param id              a unique identifier for this likelihood
     * @param exclude         taxa to exclude from this model
     */
    public SpeciationLikelihood(Tree tree, SpeciationModel speciationModel, Set<Taxon> exclude, String id) {
        this(SPECIATION_LIKELIHOOD, tree, speciationModel, exclude);
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

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return SPECIATION_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) {

            XMLObject cxo = (XMLObject) xo.getChild(MODEL);
            SpeciationModel specModel = (SpeciationModel) cxo.getChild(UltrametricSpeciationModel.class);

            cxo = (XMLObject) xo.getChild(TREE);
            Tree tree = (Tree) cxo.getChild(Tree.class);

            Set<Taxon> excludeTaxa = null;

            if (xo.hasChildNamed(INCLUDE)) {
                excludeTaxa = new HashSet<Taxon>();
                for (int i = 0; i < tree.getTaxonCount(); i++) {
                    excludeTaxa.add(tree.getTaxon(i));
                }

                cxo = (XMLObject) xo.getChild(INCLUDE);
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    TaxonList taxonList = (TaxonList) cxo.getChild(i);
                    for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                        excludeTaxa.remove(taxonList.getTaxon(j));
                    }
                }
            }

            if (xo.hasChildNamed(EXCLUDE)) {
                excludeTaxa = new HashSet<Taxon>();
                cxo = (XMLObject) xo.getChild(EXCLUDE);
                for (int i = 0; i < cxo.getChildCount(); i++) {
                    TaxonList taxonList = (TaxonList) cxo.getChild(i);
                    for (int j = 0; j < taxonList.getTaxonCount(); j++) {
                        excludeTaxa.add(taxonList.getTaxon(j));
                    }
                }
            }
            if (excludeTaxa != null) {
                Logger.getLogger("dr.evomodel").info("Speciation model excluding " + excludeTaxa.size() + " taxa from prior - " +
                        (tree.getTaxonCount() - excludeTaxa.size()) + " taxa remaining.");
            }

            return new SpeciationLikelihood(tree, specModel, excludeTaxa, null);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the speciation.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MODEL, new XMLSyntaxRule[]{
                        new ElementRule(UltrametricSpeciationModel.class)
                }),
                new ElementRule(TREE, new XMLSyntaxRule[]{
                        new ElementRule(Tree.class)
                }),

                new ElementRule(INCLUDE, new XMLSyntaxRule[]{
                        new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
                }, "One or more subsets of taxa which should be included from calculate the likelihood (the remaining taxa are excluded)", true),

                new ElementRule(EXCLUDE, new XMLSyntaxRule[]{
                        new ElementRule(Taxa.class, 1, Integer.MAX_VALUE)
                }, "One or more subsets of taxa which should be excluded from calculate the likelihood (which is calculated on the remaining subtree)", true)
        };
    };

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

    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;
}