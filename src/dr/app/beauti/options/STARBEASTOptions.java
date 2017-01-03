/*
 * STARBEASTOptions.java
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

package dr.app.beauti.options;


import dr.app.beauti.generator.Generator;
import dr.app.beauti.types.OperatorType;
import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.TreePriorType;
import dr.evomodelxml.operators.TreeNodeSlideParser;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.speciation.SpeciesTreeModelParser;
import dr.evomodelxml.speciation.YuleModelParser;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Walter Xie
 * @version $Id$
 */
@Deprecated
public class STARBEASTOptions extends ModelOptions {
    private static final long serialVersionUID = -668734758207432455L;

    public static final String TREE_FILE_NAME = "trees";

    public static final String POP_MEAN = "popMean";
    public static final String SPECIES_TREE_FILE_NAME = TraitData.TRAIT_SPECIES
    							+ "." + STARBEASTOptions.TREE_FILE_NAME; // species.trees

    public static final String CITATION = "<html>Joseph Heled and Alexei J. Drummond,<br>" +
                "Bayesian Inference of Species Trees from Multilocus Data,<br>" +
                "Molecular Biology and Evolution 2010 27(3):570-580</html>";

    public static final String EXAMPLE_FORMAT = "<html>A proper trait file is tab delimited. " +
            "The first row is always <font color=red>traits</font> followed by the keyword " +
            "(e.g. <font color=red>species</font> in *BEAST) in the second column and separated " +
            "by <font color=red>tab</font>. The rest rows are mapping taxa to species, which list " +
            "taxon name in the first column and species name in the second column separated by " +
            "<font color=red>tab</font>. For example: <br>" +
            "traits\tspecies<br>" +
            "taxon1\tspeciesA<br>" +
            "taxon2\tspeciesA<br>" +
            "taxon3\tspeciesB<br>" +
            "... ...<br>" +
            "Once mapping file is loaded, the trait named by keyword <font color=red>species</font> " +
            "is displayed in the main panel, and the message of using *BEAST is also displayed on " +
            "the bottom of main frame.<br>" +
            "For multi-alignment, the default of *BEAST is unlinking all models: substitution model, " +
            "clock model, and tree models.</html>";

    private final BeautiOptions options;

    public STARBEASTOptions(BeautiOptions options) {
        this.options = options;
        initModelParametersAndOpererators();
    }

    @Override
    public void initModelParametersAndOpererators() {
        double spWeights = 5.0;
        double spTuning = 0.9;

        createParameterOneOverXPrior(TraitData.TRAIT_SPECIES + "." + POP_MEAN, "Species tree: population hyper-parameter operator",
                PriorScaleType.TIME_SCALE, 1.0);
        // species tree Yule
        createParameterOneOverXPrior(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE,
                "Species tree: Yule process birth rate", PriorScaleType.BIRTH_RATE_SCALE, 1.0);

        // species tree Birth Death
        createParameterOneOverXPrior(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME,
                "Species tree: Birth Death model mean growth rate", PriorScaleType.BIRTH_RATE_SCALE, 1.0);
        createZeroOneParameterUniformPrior(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME,
                "Species tree: Birth Death model relative death rate", 0.5);

        createParameterOneOverXPrior(SpeciesTreeModelParser.SPECIES_TREE + "." + Generator.SPLIT_POPS, "Species tree: population size operator",
                PriorScaleType.TIME_SCALE, 1.0);

        createParameter(TraitData.TRAIT_SPECIES + "." + TreeNodeSlideParser.TREE_NODE_REHEIGHT, "Species tree: tree node operator");

        createScaleOperator(TraitData.TRAIT_SPECIES + "." + POP_MEAN, spTuning, spWeights);

        createScaleOperator(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE, demoTuning, demoWeights);
        createScaleOperator(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, demoTuning, demoWeights);
        createScaleOperator(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, demoTuning, demoWeights);

        createScaleOperator(SpeciesTreeModelParser.SPECIES_TREE + "." + Generator.SPLIT_POPS, 0.5, 94);

        createOperator(TraitData.TRAIT_SPECIES + "." + TreeNodeSlideParser.TREE_NODE_REHEIGHT, OperatorType.NODE_REHIGHT, demoTuning, 94);

        //TODO: more

    }


    @Override
    public List<Parameter> selectParameters(List<Parameter> params) {

        params.add(getParameter(TraitData.TRAIT_SPECIES + "." + POP_MEAN));

        if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
            params.add(getParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            params.add(getParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
        } else if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE
                || options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION) {
            params.add(getParameter(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE));
        }

//    	params.add(getParameter(SpeciesTreeModel.SPECIES_TREE + "." + Generator.SPLIT_POPS));

        return params;
    }

    @Override
    public List<Operator> selectOperators(List<Operator> ops) {
        ops.add(getOperator(TraitData.TRAIT_SPECIES + "." + POP_MEAN));

        if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_BIRTH_DEATH) {
            ops.add(getOperator(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            ops.add(getOperator(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));

//            ops.add(getOperator("upDownBirthDeathSpeciesTree"));
//            ops.add(getOperator("upDownBirthDeathSTPop"));
//
//            for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            	ops.add(getOperator(tree.getPrefix() + "upDownBirthDeathGeneTree"));
//            }
        } else if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE
                || options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION) {
            ops.add(getOperator(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE));

//            ops.add(getOperator("upDownYuleSpeciesTree"));
//            ops.add(getOperator("upDownYuleSTPop"));
//
//            for (PartitionTreeModel tree : getPartitionTreeModels()) {
//            	ops.add(getOperator(tree.getPrefix() + "upDownYuleGeneTree"));
//            }
        }

        ops.add(getOperator(SpeciesTreeModelParser.SPECIES_TREE + "." + Generator.SPLIT_POPS));

        ops.add(getOperator(TraitData.TRAIT_SPECIES + "." + TreeNodeSlideParser.TREE_NODE_REHEIGHT));

        return ops;
    }

    @Override
    public String getPrefix() {
        return "";
    }

    /////////////////////////////////////////////////////////////

    public List<String> getSpeciesList() {
        return new ArrayList<String>(TraitData.getStatesListOfTrait(options.taxonList, TraitData.TRAIT_SPECIES));
    }

    public int getEmptySpeciesIndex() {
        return TraitData.getEmptyStateIndex(options.taxonList, TraitData.TRAIT_SPECIES);
    }

    public String getDescription() {
        return "Species definition: binds taxa, species and gene trees";
    }

}
