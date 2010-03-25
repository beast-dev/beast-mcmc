/*
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

package dr.app.beauti.options;

import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 * @version $Id$
 */
public abstract class TraitsOptions extends PartitionOptions {

    public static enum TraitType {
        DISCRETE,
        INTEGER,
        CONTINUOUS
    }

    public static enum Traits {
        TRAIT_SPECIES("species"),
        TRAIT_LOCATIONS("locations");

        Traits(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

    public static final List<TraitGuesser> traits = new ArrayList<TraitGuesser>();  // traits list
//    public Map<String , TraitGuesser> traits = new HashMap<String, TraitGuesser>(); // traits map
    
//    public TraitsOptions(BeautiOptions options) {
//    	this.options = options;
//        initTraitParametersAndOperators();
//    }

    protected final String PREFIX_;
//    private final TraitGuesser traitGuesser;

    public TraitsOptions(TraitGuesser traitGuesser) {
//        this.traitGuesser = traitGuesser;
        this.partitionName = traitGuesser.getTraitName();
        PREFIX_ = partitionName + ".";
        initTraitParametersAndOperators();
    }

    protected abstract void initTraitParametersAndOperators();

    protected abstract void selectParameters(List<Parameter> params);
    protected abstract void selectOperators(List<Operator> ops);

    public abstract boolean isSpecifiedTraitAnalysis(String traitName);

    /////////////////////////////////////////////////////////////
    public static boolean containTrait(String traitName) {
        for (TraitGuesser trait : traits) {
            if (trait.getTraitName().equalsIgnoreCase(traitName))
                return true;
        }
        return false;
    }

    public static TraitGuesser getTrait(String traitName) {
        for (TraitGuesser trait : traits) {
            if (trait.getTraitName().equalsIgnoreCase(traitName))
                return trait;
        }
        return null;
    }
    
    public static List<String> getStatesListOfTrait(Taxa taxonList, String traitName) {
        List<String> species = new ArrayList<String>();
        String sp;

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                sp = (String) taxon.getAttribute(traitName);

                if (sp == null) return null;

                if (!species.contains(sp)) {
                    species.add(sp);
                }
            }
            return species;
        } else {
            return null;
        }
    }

    public static List<TraitGuesser> getDiscreteTraitsExcludeSpecies() { // exclude species at moment
        List<TraitGuesser> discreteTraitsExcludeSpecies = new ArrayList<TraitGuesser>();
        for (TraitGuesser trait : traits) {
            if (  (!trait.getTraitName().equalsIgnoreCase(TraitsOptions.Traits.TRAIT_SPECIES.toString()))
                    && trait.getTraitType() == TraitsOptions.TraitType.DISCRETE) {
                discreteTraitsExcludeSpecies.add(trait);
            }
        }
        return discreteTraitsExcludeSpecies;
    }

    public static boolean hasDiscreteTraitsExcludeSpecies() { // exclude species at moment
        return getDiscreteTraitsExcludeSpecies() != null && getDiscreteTraitsExcludeSpecies().size() > 0;
    }

    public static boolean hasPhylogeographic() {
        return containTrait(TraitsOptions.Traits.TRAIT_LOCATIONS.toString());
    }

    public static String getPhylogeographicDescription() {
        return "Discrete phylogeographic inference in BEAST (PLoS Comput Biol. 2009 Sep;5(9):e1000520)";
    }

    public String getPrefix() {
        return PREFIX_;
    }
}