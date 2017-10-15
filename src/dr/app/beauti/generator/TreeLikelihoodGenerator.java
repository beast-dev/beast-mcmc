/*
 * TreeLikelihoodGenerator.java
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

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.evomodelxml.treelikelihood.MarkovJumpsTreeLikelihoodParser;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.components.ancestralstates.AncestralStatesComponentOptions;
import dr.app.beauti.options.*;
import dr.app.beauti.types.MicroSatModelType;
import dr.app.beauti.util.XMLWriter;
import dr.evolution.datatype.DataType;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.oldevomodel.substmodel.AsymmetricQuadraticModel;
import dr.oldevomodel.substmodel.LinearBiasModel;
import dr.oldevomodel.substmodel.TwoPhaseModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.tree.MicrosatelliteSamplerTreeModelParser;
import dr.oldevomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.oldevomodelxml.treelikelihood.MicrosatelliteSamplerTreeLikelihoodParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.SitePatternsParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class TreeLikelihoodGenerator extends Generator {

    private static final boolean GENERATE_TREE_LIKELIHOOD = false;

    public TreeLikelihoodGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    public void writeAllTreeLikelihoods(XMLWriter writer) throws GeneratorException {
        for (List<PartitionData> partitions : options.multiPartitionLists)  {
            writeTreeDataLikelihood(partitions, writer);
            writer.writeText("");
        }

        for (AbstractPartitionData partition : options.otherPartitions) {
            // generate tree likelihoods for the other data partitions
            if (partition.getTaxonList() != null) {
                if (partition instanceof PartitionData) {
                    writeTreeLikelihood((PartitionData) partition, writer);
                    writer.writeText("");
                } else if (partition instanceof PartitionPattern) { // microsat
                    writeTreeLikelihood((PartitionPattern) partition, writer);
                    writer.writeText("");
                } else {
                    throw new GeneratorException("Find unrecognized partition:\n" + partition.getName());
                }
            }
        }
    }

    public void writeTreeDataLikelihood(List<PartitionData> partitions, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partitions.get(0).getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partitions.get(0).getPartitionTreeModel();
        PartitionClockModel clockModel = partitions.get(0).getPartitionClockModel();

        String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
        String idString = prefix + "treeLikelihood";

        Attribute[] attributes = new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, idString),
                new Attribute.Default<Boolean>(TreeDataLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood())
        };

        writer.writeComment("Likelihood for tree given sequence data");
        writer.writeOpenTag(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, attributes);

        for (PartitionData partition : partitions) {
            substModel = partition.getPartitionSubstitutionModel();
            PartitionClockModel cm = partition.getPartitionClockModel();
            if (clockModel == null) {
                clockModel = cm;
            } else if (clockModel != cm) {
                throw new RuntimeException("All the partitions in a TreeDataLikelihood should share the same clock model.");
            }

            if (substModel.getCodonHeteroPattern() != null) {

                for (int num = 1; num <= substModel.getCodonPartitionCount(); num++) {
                    String prefix1 = partition.getPrefix() + substModel.getPrefixCodon(num);

                    writer.writeOpenTag(TreeDataLikelihoodParser.PARTITION);
                    writeCodonPatternsRef(prefix1, num, substModel.getCodonPartitionCount(), writer);

                    writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix(num) + SiteModel.SITE_MODEL);

                    writer.writeCloseTag(TreeDataLikelihoodParser.PARTITION);
                }


            } else {
                String prefix1 = partition.getPrefix();
                writer.writeOpenTag(TreeDataLikelihoodParser.PARTITION);
                writer.writeIDref(SitePatternsParser.PATTERNS, prefix1 + SitePatternsParser.PATTERNS);
                writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
                writer.writeCloseTag(TreeDataLikelihoodParser.PARTITION);
            }


        }

        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);


        writer.writeCloseTag(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD);

    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param partition the partition  to write likelihood block for
     * @param writer    the writer
     */
    public void writeTreeLikelihood(PartitionData partition, XMLWriter writer) {

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model.isDolloModel()) {
            return; // DolloComponent will add tree likelihood
        }

        String treeLikelihoodTag = TreeLikelihoodParser.TREE_LIKELIHOOD;
        if (ancestralStatesOptions.usingAncestralStates(partition)) {
            treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
            if (ancestralStatesOptions.isCountingStates(partition)) {
                // State change counting uses the MarkovJumpsTreeLikelihood but
                // dNdS robust counting doesn't as it has its own counting code...
                if (!ancestralStatesOptions.dNdSRobustCounting(partition)) {
                    treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
                }
            }
        }

        if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {

            for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                writeTreeLikelihood(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, i, partition, writer);
            }


        } else {
            writeTreeLikelihood(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, -1, partition, writer);
        }
    }

    /**
     * Write the tree likelihood XML block.
     *
     * @param id        the id of the tree likelihood
     * @param num       the likelihood number
     * @param partition the partition to write likelihood block for
     * @param writer    the writer
     */
    private void writeTreeLikelihood(String tag, String id, int num, PartitionData partition, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        writer.writeComment("Likelihood for tree given sequence data");

        String prefix;
        if (num > 0) {
            prefix = partition.getPrefix() + substModel.getPrefixCodon(num);
        } else {
            prefix = partition.getPrefix();
        }

        String idString = prefix + id;

        Attribute[] attributes;
        if (tag.equals(MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD)) {
            AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                    .getComponentOptions(AncestralStatesComponentOptions.class);
            boolean saveCompleteHistory = ancestralStatesOptions.isCompleteHistoryLogging(partition);
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood()),
                    new Attribute.Default<Boolean>(MarkovJumpsTreeLikelihoodParser.USE_UNIFORMIZATION, true),
                    new Attribute.Default<Integer>(MarkovJumpsTreeLikelihoodParser.NUMBER_OF_SIMULANTS, 1),
                    new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
                    new Attribute.Default<String>(MarkovJumpsTreeLikelihoodParser.SAVE_HISTORY, saveCompleteHistory ? "true" : "false"),
            };
        } else if (tag.equals(TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD)) {
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood()),
                    new Attribute.Default<String>(AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG_NAME, prefix + AncestralStateTreeLikelihoodParser.RECONSTRUCTION_TAG),
            };
        } else {
            attributes = new Attribute[]{
                    new Attribute.Default<String>(XMLParser.ID, idString),
                    new Attribute.Default<Boolean>(TreeLikelihoodParser.USE_AMBIGUITIES, substModel.isUseAmbiguitiesTreeLikelihood())
            };
        }

        writer.writeOpenTag(tag, attributes);

        if (!options.samplePriorOnly) {
            if (num > 0) {
                writeCodonPatternsRef(prefix, num, substModel.getCodonPartitionCount(), writer);
            } else {
                writer.writeIDref(SitePatternsParser.PATTERNS, prefix + SitePatternsParser.PATTERNS);
            }
        } else {
            // We just need to use the dummy alignment
            writer.writeIDref(AlignmentParser.ALIGNMENT, partition.getAlignment().getId());
        }

        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);

        if (num > 0) {
            writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix(num) + SiteModel.SITE_MODEL);
        } else {
            writer.writeIDref(GammaSiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
        }

        ClockModelGenerator.writeBranchRatesModelRef(clockModel, writer);

        generateInsertionPoint(ComponentGenerator.InsertionPoint.IN_TREE_LIKELIHOOD, partition, prefix, writer);

        writer.writeCloseTag(tag);
    }

    /**
     * Write a list of idrefs to all treelikelihoods
     * @param writer
     */
    public void writeTreeLikelihoodReferences(XMLWriter writer) {
        for (List<PartitionData> partitions : options.multiPartitionLists)  {
            // TreeDataLikelihoods are labelled by their tree and clock models
            PartitionTreeModel treeModel = partitions.get(0).getPartitionTreeModel();
            PartitionClockModel clockModel = partitions.get(0).getPartitionClockModel();

            String prefix = treeModel.getPrefix() + clockModel.getPrefix(); // use the treemodel prefix
            String idString = prefix + "treeLikelihood";

            writer.writeIDref(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD, idString);
        }

        for (AbstractPartitionData partition : options.otherPartitions) {
            // generate tree likelihoods for the other data partitions

            // TreeLikelihood are labelled by their partition
            AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                    .getComponentOptions(AncestralStatesComponentOptions.class);

            PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

            if (model.isDolloModel()) {
                return; // DolloComponent will add tree likelihood
            }

            String treeLikelihoodTag = TreeLikelihoodParser.TREE_LIKELIHOOD;
            if (ancestralStatesOptions.usingAncestralStates(partition)) {
                treeLikelihoodTag = TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD;
                if (ancestralStatesOptions.isCountingStates(partition)) {
                    // State change counting uses the MarkovJumpsTreeLikelihood but
                    // dNdS robust counting doesn't as it has its own counting code...
                    if (!ancestralStatesOptions.dNdSRobustCounting(partition)) {
                        treeLikelihoodTag = MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD;
                    }
                }
            }

            if (model.getDataType().getType() == DataType.NUCLEOTIDES && model.getCodonHeteroPattern() != null) {

                for (int i = 1; i <= model.getCodonPartitionCount(); i++) {
                    writeTreeLikelihoodRef(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, i, partition, writer);
                }


            } else {
                writeTreeLikelihoodRef(treeLikelihoodTag, TreeLikelihoodParser.TREE_LIKELIHOOD, -1, partition, writer);
            }
        }
    }

    private void writeTreeLikelihoodRef(String tag, String id, int num, AbstractPartitionData partition, XMLWriter writer) {

        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();

        String prefix;
        if (num > 0) {
            prefix = partition.getPrefix() + substModel.getPrefixCodon(num);
        } else {
            prefix = partition.getPrefix();
        }

        String idString = prefix + id;

        writer.writeIDref(tag, idString);
    }

    /**
     * Write Microsatellite Sampler tree likelihood XML block.
     *
     * @param partition the partition  to write likelihood block for
     * @param writer    the writer
     */
    public void writeTreeLikelihood(PartitionPattern partition, XMLWriter writer) {
        PartitionSubstitutionModel substModel = partition.getPartitionSubstitutionModel();
//        PartitionTreeModel treeModel = partition.getPartitionTreeModel();
        PartitionClockModel clockModel = partition.getPartitionClockModel();

        writer.writeComment("Microsatellite Sampler Tree Likelihood");

        writer.writeOpenTag(MicrosatelliteSamplerTreeLikelihoodParser.TREE_LIKELIHOOD,
                new Attribute[]{new Attribute.Default<String>(XMLParser.ID,
                        partition.getPrefix() + MicrosatelliteSamplerTreeLikelihoodParser.TREE_LIKELIHOOD)});

        writeMicrosatSubstModelRef(substModel, writer);

        writer.writeIDref(MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL,
                partition.getName() + "." + MicrosatelliteSamplerTreeModelParser.TREE_MICROSATELLITE_SAMPLER_MODEL);

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES, clockModel.getPrefix()
                        + BranchRateModel.BRANCH_RATES);
                break;
            case UNCORRELATED:
            case RANDOM_LOCAL_CLOCK:
            case AUTOCORRELATED:
                throw new UnsupportedOperationException("Microsatellite only supports strict clock model");

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        writer.writeCloseTag(MicrosatelliteSamplerTreeLikelihoodParser.TREE_LIKELIHOOD);
    }

    public void writeMicrosatSubstModelRef(PartitionSubstitutionModel model, XMLWriter writer) {
        if (model.getPhase() != MicroSatModelType.Phase.ONE_PHASE) {
            writer.writeIDref(TwoPhaseModel.TWO_PHASE_MODEL, model.getPrefix() + TwoPhaseModel.TWO_PHASE_MODEL);
        } else if (model.getMutationBias() != MicroSatModelType.MutationalBias.UNBIASED) {
            writer.writeIDref(LinearBiasModel.LINEAR_BIAS_MODEL, model.getPrefix() + LinearBiasModel.LINEAR_BIAS_MODEL);
        } else {
            writer.writeIDref(AsymmetricQuadraticModel.ASYMQUAD_MODEL, model.getPrefix() + AsymmetricQuadraticModel.ASYMQUAD_MODEL);
        }
    }

    /**
     * Is the data partition suitable for using the MultipartitionTreeDataLikelihoodDelegate.
     *
     * @param partition the partition  to write likelihood block for
     */
    public boolean canUseMultiPartition(AbstractPartitionData partition) {

        AncestralStatesComponentOptions ancestralStatesOptions = (AncestralStatesComponentOptions) options
                .getComponentOptions(AncestralStatesComponentOptions.class);

        PartitionSubstitutionModel model = partition.getPartitionSubstitutionModel();

        if (model.isDolloModel()) {
            return false;
        }

        if (ancestralStatesOptions.usingAncestralStates(partition)) {
            return false;
        }

        if (model.getDataType().getType() != DataType.NUCLEOTIDES && model.getDataType().getType() != DataType.AMINO_ACIDS) {
            return false;
        }

        return true;
    }

}
