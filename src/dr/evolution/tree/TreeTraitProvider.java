/*
 * TreeTraitProvider.java
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

package dr.evolution.tree;

import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * An interface for objects that can provide TreeTraits (i.e., information about the nodes and
 * branches of a tree).
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
public interface TreeTraitProvider {

    /**
     * Returns an array of all the available traits
     * @return the array
     */
    TreeTrait[] getTreeTraits();

    /**
     * Returns a trait that is stored using a specific key. This will often be the same
     * as the 'name' of the trait but may not be depending on the application.
     * @param key a unique key
     * @return the trait
     */
    TreeTrait getTreeTrait(String key);

    /**
     * This is a concrete helper class that can be used either as a mix in base to another
     * class or as a delegate. It is itself a TreeTraitProvider so can be instantiated and
     * passed as is.
     */
    public class Helper implements TreeTraitProvider {

        /**
         * Default constructor
         */
        public Helper() {
        }

        /**
         * Constructor taking a single initial trait
         * @param trait the TreeTrait
         */
        public Helper(TreeTrait trait) {
           addTrait(trait);
        }

        /**
         * Constructor taking a single initial trait
         * @param trait the TreeTrait
         */
        public Helper(String key, TreeTrait trait) {
           addTrait(key, trait);
        }

        /**
         * Constructor taking an array of initial traits
         * @param traits the array of TreeTraits
         */
        public Helper(TreeTrait[] traits) {
            addTraits(traits);
        }

        /**
         * Constructor taking a collection of initial traits
         * @param traits the collection of TreeTraits
         */
        public Helper(Collection<TreeTrait> traits) {
            addTraits(traits);
        }

        /**
         * Add a single trait to the list keyed by its name
         * @param trait the TreeTrait
         */
        public void addTrait(TreeTrait trait) {
            traits.put(trait.getTraitName(), trait);
        }

        /**
         * Add a single trait to the list using a given key
         * @param key the key
         * @param trait the TreeTrait
         */
        public void addTrait(String key, TreeTrait trait) {
            if (traits.containsKey(key)) {
                throw new RuntimeException("All traits must have unique names");
            }
            traits.put(key, trait);
        }

        /**
         * Add an array of traits to the list, keying each with its
         * respective name
         * @param traits the array of TreeTraits
         */
        public void addTraits(TreeTrait[] traits) {
            for (TreeTrait trait : traits) {
                this.traits.put(trait.getTraitName(), trait);
            }
        }

        /**
         * Add a collection of traits to the list, keying each with its
         * respective name
         * @param traits the collection of TreeTraits
         */
        public void addTraits(Collection<TreeTrait> traits) {
            for (TreeTrait trait : traits) {
                this.traits.put(trait.getTraitName(), trait);
            }
        }


        // Implementation of TreeTraitProvider interface

        public TreeTrait[] getTreeTraits() {
            TreeTrait[] traitArray = new TreeTrait[traits.values().size()];
            return traits.values().toArray(traitArray);
        }

        public TreeTrait getTreeTrait(String key) {
            return traits.get(key);
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            Helper helper = (Helper) o;

            if (traits != null ? !traits.equals(helper.traits) : helper.traits != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            return traits != null ? traits.hashCode() : 0;
        }

        // Private members

        private Map<String, TreeTrait> traits = new HashMap<String, TreeTrait>();
    }
}
