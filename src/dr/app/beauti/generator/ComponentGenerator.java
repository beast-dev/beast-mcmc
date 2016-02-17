/*
 * ComponentGenerator.java
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

package dr.app.beauti.generator;

import dr.app.beauti.util.XMLWriter;

/**
 * This interface is for generators of model components that need to insert XML
 * at various places throughout the document. The implementation needs to return
 * whether it wishes to make use of an insertion point and then generate code
 * for that insertion point. Most implementations will be derived from the abstract
 * BaseComponent.
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface ComponentGenerator {
    enum InsertionPoint {
        BEFORE_TAXA,                // at the beginning of the document
        IN_TAXON,                   // in each individual taxon
        AFTER_TAXA,                 // after all taxon sets have been defined
        AFTER_SEQUENCES,            // after all alignments have been defined
        AFTER_PATTERNS,             // after all patterns
        IN_TREE_MODEL,              // in the tree model(s)
        AFTER_TREE_MODEL,           // after the tree model
        AFTER_TREE_PRIOR,           // after the tree prior
        AFTER_SUBSTITUTION_MODEL,   // after all substitution models
        AFTER_SITE_MODEL,           // after all site models
        IN_TREE_LIKELIHOOD,         // in the tree likelihood(s)
        AFTER_TREE_LIKELIHOOD,      // after all tree likelihoods
        AFTER_TRAITS,				// after each traits mapping 
        BEFORE_OPERATORS,           // before the operator schedule
        IN_OPERATORS,               // in the operator schedule
        AFTER_OPERATORS,            // after the operator schedule
        IN_MCMC_PRIOR,              // in the prior section of the MCMC
        IN_MCMC_LIKELIHOOD,         // in the likelihood section of the MCMC
        IN_SCREEN_LOG,              // in the screen log
        AFTER_SCREEN_LOG,           // after the screen log
        IN_FILE_LOG_PARAMETERS,     // in the file log after the parameters have been logged
        IN_FILE_LOG_LIKELIHOODS,    // in the file log after the likelihoods have been logged
        AFTER_FILE_LOG,             // after the file log
        IN_TREES_LOG,               // in the trees log
        AFTER_TREES_LOG,            // after the trees log
        AFTER_MCMC                 // after the mcmc element
    }

    /**
     * Opportunity for a component to do some pre-generation checks of the options selected.
     * Throwing a GeneratorException will result in a dialog box being shown to the user and
     * generation being disallowed.
     * @throws Generator.GeneratorException
     */
    void checkOptions() throws Generator.GeneratorException;

    /**
     * Returns whether this component requires access to a particular insertion point
     * @param point the insertion point
     * @return whether it requires it
     */
    boolean usesInsertionPoint(InsertionPoint point);

    /**
     * Called to allow the component to generate at the particular insertion point. For
     * some insertion points (currently only 'IN_TAXON') the specific item is given.
     * @param generator the calling generator
     * @param point the insertion point
     * @param item a reference to the item being generated (or null if not applicable)
     * @param writer the XMLWriter
     */
    void generateAtInsertionPoint(Generator generator, InsertionPoint point, Object item, XMLWriter writer);

    /**
     * Called to allow the component to generate at the particular insertion point. For
     * some insertion points the specific item is given plus a prefix to allow for unique IDs.
     * @param generator the calling generator
     * @param point the insertion point
     * @param item a reference to the item being generated (or null if not applicable)
     * @param prefix a prefix string
     * @param writer the XMLWriter
     */
    void generateAtInsertionPoint(Generator generator, InsertionPoint point, Object item, String prefix, XMLWriter writer);

}
