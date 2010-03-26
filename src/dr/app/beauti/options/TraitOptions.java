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
public abstract class TraitOptions extends PartitionModelOptions {

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

//    public static final List<TraitGuesser> traits = new ArrayList<TraitGuesser>();  // traits list
//    public Map<String , TraitGuesser> traits = new HashMap<String, TraitGuesser>(); // traits map
    
//    public TraitsOptions(BeautiOptions options) {
//    	this.options = options;
//        initTraitParametersAndOperators();
//    }

    protected final String PREFIX_;
//    private final TraitGuesser traitGuesser;

    public TraitOptions(TraitData partition) {
        this.partitionName = partition.getName();

        allPartitionData.clear();
        addPartitionData(partition);

        PREFIX_ = partitionName + ".";
        initTraitParametersAndOperators();
    }

    protected abstract void initTraitParametersAndOperators();

    protected abstract void selectParameters(List<Parameter> params);
    protected abstract void selectOperators(List<Operator> ops);

    public abstract boolean isSpecifiedTraitAnalysis(String traitName);

    /////////////////////////////////////////////////////////////

    public static List<String> getStatesListOfTrait(Taxa taxonList, String traitName) {
        List<String> states = new ArrayList<String>();
        String attr;

        if (taxonList != null) {
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                attr = (String) taxon.getAttribute(traitName);

                if (attr == null) return null;

                if (!states.contains(attr)) {
                    states.add(attr);
                }
            }
            return states;
        } else {
            return null;
        }
    }

    public String getPrefix() {
        return PREFIX_;
    }
}