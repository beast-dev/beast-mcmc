/*
 * STARBEASTGenerator.java
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

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.TraitData;
import dr.app.beauti.types.PopulationSizeModelType;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.PloidyType;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.speciation.CalibrationPoints;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.coalescent.OldCoalescentSimulatorParser;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.speciation.*;
import dr.evomodelxml.tree.TMRCAStatisticParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.evoxml.util.XMLUnits;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.model.ParameterParser;
import dr.inferencexml.distribution.DistributionModelParser;
import dr.inferencexml.distribution.MixedDistributionLikelihoodParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class STARBEASTGenerator extends Generator {

    public static final String ALL_SPECIES = "allSpecies";
    private int numOfSpecies; // used in private String getIndicatorsParaValue()

    public STARBEASTGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    /**
     * *BEAST block, write the species, species tree, species tree model, likelihood, etc.
     *
     * @param writer XMLWriter
     */
    public void writeSTARBEAST(XMLWriter writer) {
        writeSpeciesTreePrior(writer);
        writeSpeciesTreeLikelihood(writer);
        writeSpeciesTreeRootHeight(writer);
        writeGeneUnderSpecies(writer);
    }

    public void writeSpecies(XMLWriter writer) {
        String traitName = TraitData.TRAIT_SPECIES;
        writer.writeText("");
        writer.writeComment(options.starBEASTOptions.getDescription());

        writer.writeOpenTag(traitName, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitName)});
        //new Attribute.Default<String>("traitType", traitType)});
        writeMultiSpecies(options.taxonList, writer);
        writer.writeCloseTag(traitName);

        writer.writeText("");
        writer.writeComment("full species set for species tree root height");
        writer.writeOpenTag(TaxaParser.TAXA, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, ALL_SPECIES)});
        for (String eachSp : options.starBEASTOptions.getSpeciesList()) {
            writer.writeIDref(SpeciesBindingsSPinfoParser.SP, eachSp);
        }
        writer.writeCloseTag(TaxaParser.TAXA);
    }

    public void writeStartingTreeForCalibration(XMLWriter writer) {
        writer.writeComment("species starting tree for calibration");
        writer.writeText("");
        writer.writeOpenTag(OldCoalescentSimulatorParser.COALESCENT_TREE,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID, SP_START_TREE)}
        );

        Attribute[] taxaAttribute = {new Attribute.Default<String>(XMLParser.IDREF, ALL_SPECIES)};

        writer.writeOpenTag(OldCoalescentSimulatorParser.CONSTRAINED_TAXA);
        writer.writeTag(TaxaParser.TAXA, taxaAttribute, true);
        for (Taxa taxa : options.speciesSets) {
            Parameter statistic = options.getStatistic(taxa);

            Attribute mono = new Attribute.Default<Boolean>(
                    OldCoalescentSimulatorParser.IS_MONOPHYLETIC, options.speciesSetsMono.get(taxa));

            writer.writeOpenTag(OldCoalescentSimulatorParser.TMRCA_CONSTRAINT, mono);

            writer.writeIDref(TaxaParser.TAXA, taxa.getId());

            if (options.getPartitionTreePriors().get(0).getNodeHeightPrior() == TreePriorType.SPECIES_YULE_CALIBRATION
                    && statistic.priorType == PriorType.UNIFORM_PRIOR) {
                writeDistribution(statistic, false, writer);
            }

            writer.writeCloseTag(OldCoalescentSimulatorParser.TMRCA_CONSTRAINT);

        }
        writer.writeCloseTag(OldCoalescentSimulatorParser.CONSTRAINED_TAXA);


        writer.writeOpenTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL,
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, "spInitDemo"),
                        new Attribute.Default<String>("units", Units.Utils.getDefaultUnitName(options.units))
                });
        writer.writeOpenTag(ConstantPopulationModelParser.POPULATION_SIZE);

        double popSizeValue = options.getPartitionTreePriors().get(0).getParameter("constant.popSize").getInitial(); // "initial" is "value"

        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, "sp.popSize"),
                new Attribute.Default<Double>(ParameterParser.VALUE, popSizeValue)
        }, true);
        writer.writeCloseTag(ConstantPopulationModelParser.POPULATION_SIZE);
        writer.writeCloseTag(ConstantPopulationModelParser.CONSTANT_POPULATION_MODEL);

        writer.writeCloseTag(OldCoalescentSimulatorParser.COALESCENT_TREE);
    }


    /**
     * write tag <sp>
     *
     * @param taxonList TaxonList
     * @param writer    XMLWriter
     */
    private void writeMultiSpecies(TaxonList taxonList, XMLWriter writer) {
        List<String> species = options.starBEASTOptions.getSpeciesList();
        String sp;

        numOfSpecies = species.size(); // used in private String getIndicatorsParaValue()

        for (String eachSp : species) {
            writer.writeOpenTag(SpeciesBindingsSPinfoParser.SP, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, eachSp)});

            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = null;
                try {
                    taxon = taxonList.getTaxon(i);
                    sp = taxon.getAttribute(TraitData.TRAIT_SPECIES).toString();
                } catch (Exception e) {
                    throw new IllegalArgumentException("Cannot get value from Taxon " + taxon.getId());
                }

                if (sp.equals(eachSp)) {
                    writer.writeIDref(TaxonParser.TAXON, taxon.getId());
                }
            }
            writer.writeCloseTag(SpeciesBindingsSPinfoParser.SP);
        }

        writeGeneTrees(writer);
    }

    private void writeGeneTrees(XMLWriter writer) {
        writer.writeComment("Collection of Gene Trees");

        writer.writeOpenTag(SpeciesBindingsParser.GENE_TREES, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, SpeciesBindingsParser.GENE_TREES)});

        boolean isSameAllPloidyType = true;
        PloidyType checkSamePloidyType = options.getPartitionTreeModels().get(0).getPloidyType();
        for (PartitionTreeModel model : options.getPartitionTreeModels()) {
            if (checkSamePloidyType != model.getPloidyType()) {
                isSameAllPloidyType = false;
                break;
            }
        }

        if (isSameAllPloidyType) {
            // generate gene trees regarding each data partition
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                writer.writeIDref(TreeModel.TREE_MODEL, model.getPrefix() + TreeModel.TREE_MODEL);
            }
        } else {
            // give ploidy
            for (PartitionTreeModel model : options.getPartitionTreeModels()) {
                writer.writeOpenTag(SpeciesBindingsParser.GTREE, new Attribute[]{
                        new Attribute.Default<String>(SpeciesBindingsParser.PLOIDY, Double.toString(model.getPloidyType().getValue()))
                }
                );
                writer.writeIDref(TreeModel.TREE_MODEL, model.getPrefix() + TreeModel.TREE_MODEL);
                writer.writeCloseTag(SpeciesBindingsParser.GTREE);
            }
        }

        writer.writeCloseTag(SpeciesBindingsParser.GENE_TREES);
    }


    public void writeSpeciesTree(XMLWriter writer, boolean calibration) {
        writer.writeComment("Species Tree: Provides Per branch demographic function");

        List<Attribute> attributes = new ArrayList<Attribute>();

        attributes.add(new Attribute.Default<String>(XMLParser.ID, SP_TREE));
        // *BEAST always share same tree prior
        if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS_CONSTANT) {
            attributes.add(new Attribute.Default<String>(SpeciesTreeModelParser.CONST_ROOT_POPULATION, "true"));
        } else if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONSTANT) {
            attributes.add(new Attribute.Default<String>(SpeciesTreeModelParser.CONSTANT_POPULATION, "true"));
        }

        writer.writeOpenTag(SpeciesTreeModelParser.SPECIES_TREE, attributes);

        writer.writeIDref(TraitData.TRAIT_SPECIES, TraitData.TRAIT_SPECIES);

        // take sppSplitPopulations value from partionModel(?).constant.popSize
        // *BEAST always share same tree prior
        double popSizeValue = options.getPartitionTreePriors().get(0).getParameter("constant.popSize").getInitial(); // "initial" is "value"
        writer.writeOpenTag(SpeciesTreeModelParser.SPP_SPLIT_POPULATIONS, new Attribute[]{
                new Attribute.Default<Double>(ParameterParser.VALUE, popSizeValue)});

        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SpeciesTreeModelParser.SPECIES_TREE + "." + SPLIT_POPS)}, true);

        writer.writeCloseTag(SpeciesTreeModelParser.SPP_SPLIT_POPULATIONS);

        if (calibration) writer.writeIDref(OldCoalescentSimulatorParser.COALESCENT_TREE, SP_START_TREE);

        writer.writeCloseTag(SpeciesTreeModelParser.SPECIES_TREE);
    }

    private void writeSpeciesTreePrior(XMLWriter writer) {
        Parameter para;

        TreePriorType nodeHeightPrior = options.getPartitionTreePriors().get(0).getNodeHeightPrior();

        if (nodeHeightPrior == TreePriorType.SPECIES_BIRTH_DEATH) {
            writer.writeComment("Species tree prior: Birth Death Model");

            writer.writeOpenTag(BirthDeathModelParser.BIRTH_DEATH_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, BirthDeathModelParser.BIRTH_DEATH),
                    new Attribute.Default<String>(XMLUnits.UNITS, XMLUnits.SUBSTITUTIONS)});

            writer.writeOpenTag(BirthDeathModelParser.BIRTHDIFF_RATE);
            para = options.starBEASTOptions.getParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME);
//            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
//                    new Attribute.Default<String>(XMLParser.ID, TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME),
//                    new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(para.initial)),
//                    new Attribute.Default<String>(ParameterParser.LOWER, Double.toString(para.lower)),
//                    new Attribute.Default<String>(ParameterParser.UPPER, Double.toString(para.upper))}, true);
            writeParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, para, writer);
            writer.writeCloseTag(BirthDeathModelParser.BIRTHDIFF_RATE);

            writer.writeOpenTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);
            para = options.starBEASTOptions.getParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME);
//            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
//                    new Attribute.Default<String>(XMLParser.ID, TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME),
//                    new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(para.initial)),
//                    new Attribute.Default<String>(ParameterParser.LOWER, Double.toString(para.lower)),
//                    new Attribute.Default<String>(ParameterParser.UPPER, Double.toString(para.upper))}, true);
            writeParameter(TraitData.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, para, writer);
            writer.writeCloseTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);

            writer.writeCloseTag(BirthDeathModelParser.BIRTH_DEATH_MODEL);
        } else if (nodeHeightPrior == TreePriorType.SPECIES_YULE ||
                nodeHeightPrior == TreePriorType.SPECIES_YULE_CALIBRATION) {
            writer.writeComment("Species tree prior: Yule Model");

            writer.writeOpenTag(YuleModelParser.YULE_MODEL, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, YuleModelParser.YULE),
                    new Attribute.Default<String>(XMLUnits.UNITS, XMLUnits.SUBSTITUTIONS)});

            writer.writeOpenTag(YuleModelParser.BIRTH_RATE);
            para = options.starBEASTOptions.getParameter(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE);
//            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
//                    new Attribute.Default<String>(XMLParser.ID, TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE),
//                    new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(para.initial)),
//                    new Attribute.Default<String>(ParameterParser.LOWER, Double.toString(para.lower)),
//                    new Attribute.Default<String>(ParameterParser.UPPER, Double.toString(para.upper))}, true);
            writeParameter(TraitData.TRAIT_SPECIES + "." + YuleModelParser.YULE + "." + YuleModelParser.BIRTH_RATE, para, writer);
            writer.writeCloseTag(YuleModelParser.BIRTH_RATE);

            writer.writeCloseTag(YuleModelParser.YULE_MODEL);
        } else if (nodeHeightPrior == TreePriorType.SPECIES_YULE_CALIBRATION) {

        } else {
            throw new IllegalArgumentException("Get wrong species tree prior using *BEAST : " + nodeHeightPrior.toString());
        }

    }


    private void writeSpeciesTreeLikelihood(XMLWriter writer) {

        TreePriorType nodeHeightPrior = options.getPartitionTreePriors().get(0).getNodeHeightPrior();

        if (nodeHeightPrior == TreePriorType.SPECIES_BIRTH_DEATH) {
            writer.writeComment("Species Tree Likelihood: Birth Death Model");

            writer.writeOpenTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, SPECIATION_LIKE)});

            writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
            writer.writeIDref(BirthDeathModelParser.BIRTH_DEATH_MODEL, BirthDeathModelParser.BIRTH_DEATH);
            writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);

        } else if (nodeHeightPrior == TreePriorType.SPECIES_YULE ||
                nodeHeightPrior == TreePriorType.SPECIES_YULE_CALIBRATION) {
            writer.writeComment("Species Tree Likelihood: Yule Model");

            writer.writeOpenTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD, new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, SPECIATION_LIKE)});

            writer.writeOpenTag(SpeciationLikelihoodParser.MODEL);
            writer.writeIDref(YuleModelParser.YULE_MODEL, YuleModelParser.YULE);
            writer.writeCloseTag(SpeciationLikelihoodParser.MODEL);

            if (nodeHeightPrior == TreePriorType.SPECIES_YULE_CALIBRATION) {
                // should be only 1 calibrated node with monophyletic for species tree at moment
                if (options.speciesSets.size() == 1 && options.speciesSetsMono.size() == 1) {
                    Taxa t = options.speciesSets.get(0);
                    Parameter nodeCalib = options.getStatistic(t);

                    writer.writeOpenTag(SpeciationLikelihoodParser.CALIBRATION,
                            new Attribute[]{
                                    new Attribute.Default<String>(SpeciationLikelihoodParser.CORRECTION, CalibrationPoints.CorrectionType.EXACT.toString())
                            });
                    writer.writeOpenTag(SpeciationLikelihoodParser.POINT);

                    writer.writeIDref(TaxaParser.TAXA, t.getId());
                    writeDistribution(nodeCalib, true, writer);

                    writer.writeCloseTag(SpeciationLikelihoodParser.POINT);
                    writer.writeCloseTag(SpeciationLikelihoodParser.CALIBRATION);

                    if (!options.treeModelOptions.isNodeCalibrated(nodeCalib)) {

                    }
                } else {
                    throw new IllegalArgumentException("Calibrated Yule model is only applied to 1 calibrated node with monophyletic for species tree at moment !");
                }

            }

        } else {
            throw new IllegalArgumentException("Get wrong species tree prior using *BEAST : "
                    + nodeHeightPrior.toString());
        }

        // <sp> tree
        writer.writeOpenTag(SpeciesTreeModelParser.SPECIES_TREE);
        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);
        writer.writeCloseTag(SpeciesTreeModelParser.SPECIES_TREE);

        writer.writeCloseTag(SpeciationLikelihoodParser.SPECIATION_LIKELIHOOD);
    }

    private void writeSpeciesTreeRootHeight(XMLWriter writer) {
        writer.writeComment("Species Tree: tmrcaStatistic");

        writer.writeOpenTag(TMRCAStatisticParser.TMRCA_STATISTIC, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SpeciesTreeModelParser.SPECIES_TREE + "." + TreeModelParser.ROOT_HEIGHT),
                new Attribute.Default<String>(AttributeParser.NAME, SpeciesTreeModelParser.SPECIES_TREE + "." + TreeModelParser.ROOT_HEIGHT)});

        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);

        writer.writeOpenTag(TMRCAStatisticParser.MRCA);
        writer.writeIDref(TaxaParser.TAXA, ALL_SPECIES);
        writer.writeCloseTag(TMRCAStatisticParser.MRCA);
        writer.writeCloseTag(TMRCAStatisticParser.TMRCA_STATISTIC);

    }

    private void writeGeneUnderSpecies(XMLWriter writer) {

        writer.writeComment("Species Tree: Coalescent likelihood for gene trees under species tree");

        // speciesCoalescent id="coalescent"
        writer.writeOpenTag(MultiSpeciesCoalescentParser.SPECIES_COALESCENT, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, TraitData.TRAIT_SPECIES + "." + COALESCENT)});

        writer.writeIDref(TraitData.TRAIT_SPECIES, TraitData.TRAIT_SPECIES);
        writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);

        writer.writeCloseTag(MultiSpeciesCoalescentParser.SPECIES_COALESCENT);

        // exponentialDistributionModel id="pdist"
//        writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, new Attribute[]{
//                new Attribute.Default<String>(XMLParser.ID, PDIST)});
//
//        writer.writeOpenTag(DistributionModelParser.MEAN);
//
//        Parameter para = options.getParameter(TraitGuesser.Traits.TRAIT_SPECIES + "." + options.POP_MEAN);
//
//        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
//                new Attribute.Default<String>(XMLParser.ID, TraitGuesser.Traits.TRAIT_SPECIES + "." + options.POP_MEAN),
//                new Attribute.Default<String>(ParameterParser.VALUE, Double.toString(para.initial))}, true);
//
//        writer.writeCloseTag(DistributionModelParser.MEAN);
//
//        writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL);

//        if (options.speciesTreePrior == TreePriorType.SPECIES_YULE) {

        writer.writeComment("Species tree prior: gama2 + gamma4");
        writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, SPOPS)});

        // change exponential + gamma2 into gama2 + gamma4
        // <distribution0>
        writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);
//        writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, PDIST);
        writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
        writer.writeOpenTag(DistributionModelParser.SHAPE);
        writer.writeText("2");
        writer.writeCloseTag(DistributionModelParser.SHAPE);

        writer.writeOpenTag(DistributionModelParser.SCALE);

        Parameter para = options.starBEASTOptions.getParameter(TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
        writeParameter(para, 1, writer);
        writer.writeCloseTag(DistributionModelParser.SCALE);

        writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
        writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION0);

        // <distribution1>
        writer.writeOpenTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);
        writer.writeOpenTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);

        writer.writeOpenTag(DistributionModelParser.SHAPE);
        writer.writeText("4");
        writer.writeCloseTag(DistributionModelParser.SHAPE);

        writer.writeOpenTag(DistributionModelParser.SCALE);
        writer.writeIDref(ParameterParser.PARAMETER, TraitData.TRAIT_SPECIES + "." + options.starBEASTOptions.POP_MEAN);
        writer.writeCloseTag(DistributionModelParser.SCALE);

        writer.writeCloseTag(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL);
        writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION1);

        // <data>
        writer.writeOpenTag(MixedDistributionLikelihoodParser.DATA);

        writer.writeIDref(ParameterParser.PARAMETER, SpeciesTreeModelParser.SPECIES_TREE + "." + SPLIT_POPS);

        writer.writeCloseTag(MixedDistributionLikelihoodParser.DATA);

        // <indicators>
        writer.writeOpenTag(MixedDistributionLikelihoodParser.INDICATORS);
        // Needs special treatment - you have to generate "NS" ones and 2(N-1) zeros, where N is the number of species.
        // N "1", 2(N-1) "0"
        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{new Attribute.Default<String>(ParameterParser.VALUE, getIndicatorsParaValue())}, true);

        writer.writeCloseTag(MixedDistributionLikelihoodParser.INDICATORS);

        writer.writeCloseTag(MixedDistributionLikelihoodParser.DISTRIBUTION_LIKELIHOOD);

//        } else {
//            // STPopulationPrior id="stp" log_root="true"
//            writer.writeOpenTag(SpeciesTreeBMPrior.STPRIOR, new Attribute[]{
//                    new Attribute.Default<String>(XMLParser.ID, STP),
//                    new Attribute.Default<String>(SpeciesTreeBMPrior.LOG_ROOT, "true")});
//            writer.writeIDref(SpeciesTreeModelParser.SPECIES_TREE, SP_TREE);
//
//            writer.writeOpenTag(SpeciesTreeBMPrior.TIPS);
//
//            writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, PDIST);
//
//            writer.writeCloseTag(SpeciesTreeBMPrior.TIPS);
//
//            writer.writeOpenTag(SpeciesTreeBMPrior.STSIGMA);
//
//            writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
//                    // <parameter id="stsigma" value="1" />
//                    new Attribute.Default<String>(XMLParser.ID, SpeciesTreeBMPrior.STSIGMA.toLowerCase()),
//                    new Attribute.Default<String>(ParameterParser.VALUE, "1")}, true);
//
//            writer.writeCloseTag(SpeciesTreeBMPrior.STSIGMA);
//
//            writer.writeCloseTag(SpeciesTreeBMPrior.STPRIOR);
//        }
    }

    private String getIndicatorsParaValue() {
        String v = "";

        // CONTINUOUS_CONSTANT    N  1      2(N-1) 0
        // CONTINUOUS             N  1      2N-1   0
        // CONSTANT                         2N-1   0
        if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS_CONSTANT
                || options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS) {
            for (int i = 0; i < numOfSpecies; i++) {
                if (i == (numOfSpecies - 1)) {
                    v = v + "1"; // N 1
                } else {
                    v = v + "1 "; // N 1
                }
            }
        }

        if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS_CONSTANT) {
            for (int i = 0; i < (numOfSpecies - 1); i++) {
                v = v + " 0 0"; // 2(N-1) 0
            }
        } else if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONTINUOUS) {
            v = v + " 0"; // 1   0
            for (int i = 0; i < (numOfSpecies - 1); i++) {
                v = v + " 0 0"; // 2(N-1) 0
            }
        } else if (options.getPartitionTreePriors().get(0).getPopulationSizeModel() == PopulationSizeModelType.CONSTANT) {
            v = v + "0"; // 1   0
            for (int i = 0; i < (numOfSpecies - 1); i++) {
                v = v + " 0 0"; // 2(N-1) 0
            }
        }

        return v;
    }
}
