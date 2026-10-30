/*
 * BeastSettingsImporter.java
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

package dr.app.beauti.util;

import dr.app.beauti.options.AbstractPartitionData;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionClockModel;
import dr.app.beauti.options.PartitionSubstitutionModel;
import dr.app.beauti.options.PartitionTreeModel;
import dr.app.beauti.options.PartitionTreePrior;
import dr.app.beauti.types.ClockDistributionType;
import dr.app.beauti.types.ClockType;
import dr.app.beauti.types.FrequencyPolicyType;
import dr.app.beauti.types.TreePriorParameterizationType;
import dr.app.beauti.types.TreePriorType;
import dr.evomodelxml.branchratemodel.ContinuousBranchRatesParser;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.evomodelxml.branchratemodel.MixtureModelBranchRatesParser;
import dr.evomodelxml.branchratemodel.RandomLocalClockModelParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.coalescent.GMRFSkyrideLikelihoodParser;
import dr.evomodelxml.coalescent.demographicmodel.ExponentialGrowthModelParser;
import dr.evomodelxml.speciation.Gernhard08BirthDeathModelParser;
import dr.inference.distribution.GammaDistributionModel;
import dr.inferencexml.distribution.LogNormalDistributionModelParser;
import dr.evomodel.substmodel.aminoacid.AminoAcidModelType;
import dr.evomodel.substmodel.nucleotide.NucModelType;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodelxml.substmodel.EmpiricalAminoAcidModelParser;
import dr.evomodelxml.substmodel.FrequencyModelParser;
import dr.evomodelxml.substmodel.GTRParser;
import dr.inference.model.ParameterParser;
import dr.inferencexml.operators.SimpleOperatorScheduleParser;
import dr.evomodelxml.siteratemodel.SiteModelParser;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.evomodelxml.treelikelihood.MarkovJumpsTreeLikelihoodParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.MergePatternsParser;
import dr.evoxml.SitePatternsParser;
import dr.inferencexml.MCMCParser;
import dr.oldevomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.inferencexml.loggers.CheckpointLoggerParser;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.XMLParser;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * Reads the settings of an analysis out of a BEAST XML file and applies them to a set of BEAUti options,
 * so that a file BEAUti generated can be opened again as the analysis it describes.
 * <p>
 * BEAUti writes its XML with fixed element names and a strict convention for ids
 * ({partition}.{codon position}.{name}), which is what makes this possible. A file written by hand, or by
 * another program, is read as far as it follows those conventions: anything not understood is noted in the
 * {@link ImportReport} and left at BEAUti's default rather than treated as an error.
 *
 * @author Andrew Rambaut
 */
public class BeastSettingsImporter {

    /**
     * The codon position prefix BEAUti gives the merged first and second positions.
     */
    private static final String CODON_1_AND_2 = "CP1+2.";

    /**
     * Every clock model's id ends with this, whatever kind of clock it is.
     */
    private static final String BRANCH_RATES = "branchRates";

    private final Element root;
    private final BeautiOptions options;
    private final ImportReport report;

    /**
     * Every element of the document that has an id, by id.
     */
    private final Map<String, Element> elementsById = new HashMap<String, Element>();

    public BeastSettingsImporter(Element root, BeautiOptions options, ImportReport report) {
        this.root = root;
        this.options = options;
        this.report = report;
        indexElements(root);
    }

    /**
     * Reads every setting this importer understands and applies it to the options.
     */
    public void importSettings() {
        importPartitions();
        importLinkage();
        importSubstitutionModels();
        importClockModels();
        importTreePriors();
        importTipDates();
        importMCMCSettings();
    }

    // ---- the clock models -------------------------------------------------------------------------

    /**
     * Reads the kind of molecular clock each clock model uses from the branch rates model named after it.
     */
    private void importClockModels() {
        for (PartitionClockModel model : options.getPartitionClockModels()) {
            Element branchRates = elementById(model.getPrefix() + BRANCH_RATES);
            if (branchRates == null) {
                report.notRecovered("the clock model of '" + model.getName() + "'");
                continue;
            }
            importClockModel(model, branchRates);
        }
    }

    private void importClockModel(PartitionClockModel model, Element branchRates) {
        String name = branchRates.getName();

        if (name.equalsIgnoreCase(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES)) {
            setClockType(model, ClockType.STRICT_CLOCK);
        } else if (name.equalsIgnoreCase(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES)
                || name.equalsIgnoreCase(ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES)) {
            model.setContinuousQuantile(name.equalsIgnoreCase(ContinuousBranchRatesParser.CONTINUOUS_BRANCH_RATES));
            setClockType(model, ClockType.UNCORRELATED);
            importClockDistribution(model, branchRates);
        } else if (name.equalsIgnoreCase(MixtureModelBranchRatesParser.MIXTURE_MODEL_BRANCH_RATES)) {
            model.setPerformModelAveraging(true);
            setClockType(model, ClockType.UNCORRELATED);
        } else if (name.equalsIgnoreCase(RandomLocalClockModelParser.LOCAL_BRANCH_RATES)) {
            setClockType(model, ClockType.RANDOM_LOCAL_CLOCK);
        } else if (name.equalsIgnoreCase(LocalClockModelParser.LOCAL_CLOCK_MODEL)) {
            setClockType(model, ClockType.FIXED_LOCAL_CLOCK);
        } else {
            report.notRecovered("the clock model of '" + model.getName() + "', which the file gives as " + name);
        }
    }

    /**
     * The distribution the branch rates are drawn from is the one child of the branch rates model that is
     * a distribution.
     */
    private void importClockDistribution(PartitionClockModel model, Element branchRates) {
        Element distribution = firstChild(branchRates, DiscretizedBranchRatesParser.DISTRIBUTION);
        if (distribution == null) {
            return;
        }
        for (Object child : distribution.getChildren()) {
            String name = ((Element) child).getName();
            if (name.equalsIgnoreCase(LogNormalDistributionModelParser.LOGNORMAL_DISTRIBUTION_MODEL)) {
                model.setClockDistributionType(ClockDistributionType.LOGNORMAL);
                return;
            } else if (name.equalsIgnoreCase(GammaDistributionModel.GAMMA_DISTRIBUTION_MODEL)) {
                model.setClockDistributionType(ClockDistributionType.GAMMA);
                return;
            } else if (name.toLowerCase().startsWith("exponentialdistribution")) {
                model.setClockDistributionType(ClockDistributionType.EXPONENTIAL);
                return;
            }
        }
    }

    private void setClockType(PartitionClockModel model, ClockType type) {
        model.setClockType(type);
        report.recovered("the clock model of '" + model.getName() + "', " + type);
    }

    // ---- the tree priors --------------------------------------------------------------------------

    /**
     * Reads the tree prior, which is named after the prior it belongs to. The demographic model written to
     * seed the starting tree is named differently, so it is not mistaken for the prior itself.
     */
    private void importTreePriors() {
        for (PartitionTreePrior prior : options.getPartitionTreePriors()) {
            String prefix = prior.getPrefix();
            if (setTreePriorIfPresent(prior, prefix + "constant", TreePriorType.CONSTANT)
                    || setTreePriorIfPresent(prior, prefix + "exponential", TreePriorType.EXPONENTIAL)
                    || setTreePriorIfPresent(prior, prefix + "logistic", TreePriorType.LOGISTIC)
                    || setTreePriorIfPresent(prior, prefix + "expansion", TreePriorType.EXPANSION)
                    || setTreePriorIfPresent(prior, prefix + "yule", TreePriorType.YULE)
                    || setTreePriorIfPresent(prior, prefix + "birthDeath", TreePriorType.BIRTH_DEATH)
                    || setTreePriorIfPresent(prior, prefix + "bdss", TreePriorType.BIRTH_DEATH_SERIAL_SAMPLING)
                    || setTreePriorIfPresent(prior, prefix + "skyride", TreePriorType.GMRF_SKYRIDE)
                    || setTreePriorIfPresent(prior, prefix + "skygrid", TreePriorType.SKYGRID)) {
                continue;
            }
            report.notRecovered("the tree prior of '" + prior.getName() + "'");
        }
    }

    /**
     * @return true if the file has this model, in which case the prior is set from it
     */
    private boolean setTreePriorIfPresent(PartitionTreePrior prior, String id, TreePriorType type) {
        Element element = elementById(id);
        if (element == null) {
            return false;
        }

        TreePriorType actualType = type;
        if (type == TreePriorType.BIRTH_DEATH
                && firstChild(element, Gernhard08BirthDeathModelParser.SAMPLE_PROB) != null) {
            actualType = TreePriorType.BIRTH_DEATH_INCOMPLETE_SAMPLING;
        }
        prior.setNodeHeightPrior(actualType);
        report.recovered("the tree prior of '" + prior.getName() + "', " + actualType);

        if (type == TreePriorType.EXPONENTIAL || type == TreePriorType.LOGISTIC || type == TreePriorType.EXPANSION) {
            boolean doublingTime = firstChild(element, ExponentialGrowthModelParser.DOUBLING_TIME) != null;
            prior.setParameterization(doublingTime
                    ? TreePriorParameterizationType.DOUBLING_TIME : TreePriorParameterizationType.GROWTH_RATE);
        }
        if (type == TreePriorType.SKYGRID) {
            importSkygrid(prior, element);
        }
        return true;
    }

    /**
     * The number of grid points is the dimension of the population sizes, and the cut-off is the time the
     * grid spans.
     */
    private void importSkygrid(PartitionTreePrior prior, Element skygrid) {
        Element populationSizes = firstChild(skygrid, GMRFSkyrideLikelihoodParser.POPULATION_PARAMETER);
        Integer count = populationSizes == null ? null
                : integerAttribute(firstChild(populationSizes, ParameterParser.PARAMETER), "dimension");
        if (count != null) {
            prior.setSkyGridCount(count);
        }

        Element cutOff = firstChild(skygrid, GMRFSkyrideLikelihoodParser.CUT_OFF);
        Double interval = cutOff == null ? null
                : doubleAttribute(firstChild(cutOff, ParameterParser.PARAMETER), "value");
        if (interval != null) {
            prior.setSkyGridInterval(interval);
        }
    }

    // ---- the dates on the tips --------------------------------------------------------------------

    /**
     * The dates themselves come in with the taxa. Whether the tips are calibrated with them is shown by
     * the statistic for the age of the root, which is only written for a tree whose tips are dated.
     */
    private void importTipDates() {
        if (!options.useTipDates) {
            return;
        }
        options.calculateTipHeights();

        boolean anyCalibrated = false;
        for (PartitionTreeModel treeModel : options.getPartitionTreeModels()) {
            boolean calibrated = elementById(treeModel.getPrefix() + "age(root)") != null;
            treeModel.setTipCalibrations(calibrated);
            anyCalibrated = anyCalibrated || calibrated;
        }
        if (anyCalibrated) {
            report.recovered("the dates of the tips");
        }
    }

    // ---- the substitution and site models ---------------------------------------------------------

    /**
     * Reads each substitution model and the site heterogeneity that goes with it. The site model of a
     * partition is named by the same convention BEAUti writes it with, so the model's own prefix finds it.
     */
    private void importSubstitutionModels() {
        for (PartitionSubstitutionModel model : options.getPartitionSubstitutionModels()) {
            int positions = model.getCodonPartitionCount();
            boolean found = false;
            for (int num = 1; num <= positions; num++) {
                Element siteModel = elementById(model.getPrefix(num) + SiteModelParser.SITE_MODEL);
                if (siteModel == null) {
                    continue;
                }
                found = true;
                if (num == 1) {
                    importSubstitutionModel(model, siteModel);
                    importSiteHeterogeneity(model, siteModel);
                }
            }
            if (!found) {
                report.notRecovered("the substitution model of '" + model.getName() + "'");
            }
        }
    }

    /**
     * Reads which model of sequence evolution this is, and how its base frequencies are set.
     */
    private void importSubstitutionModel(PartitionSubstitutionModel model, Element siteModel) {
        Element substitutionModel = firstChild(siteModel, SiteModelParser.SUBSTITUTION_MODEL);
        if (substitutionModel == null) {
            report.notRecovered("the substitution model of '" + model.getName() + "'");
            return;
        }

        Element referenced = null;
        String id = null;
        for (Object child : substitutionModel.getChildren()) {
            Element element = (Element) child;
            id = element.getAttributeValue(XMLParser.IDREF);
            referenced = id == null ? element : elementById(id);
            break;
        }
        if (referenced == null) {
            report.notRecovered("the substitution model of '" + model.getName() + "'");
            return;
        }

        String name = referenced.getName();
        if (name.equalsIgnoreCase(NucModelType.HKY.getXMLName())) {
            // JC is written as an HKY model whose kappa is fixed at one; only the id says which it is
            boolean isJC = id != null && id.endsWith("jc");
            setNucleotideModel(model, isJC ? NucModelType.JC : NucModelType.HKY);
        } else if (name.equalsIgnoreCase(NucModelType.TN93.getXMLName())) {
            setNucleotideModel(model, NucModelType.TN93);
        } else if (name.equalsIgnoreCase(GTRParser.GTR_MODEL)) {
            setNucleotideModel(model, NucModelType.GTR);
        } else if (name.equalsIgnoreCase(EmpiricalAminoAcidModelParser.EMPIRICAL_AMINO_ACID_MODEL)) {
            setAminoAcidModel(model, referenced.getAttributeValue("type"));
        } else {
            report.notRecovered("the substitution model of '" + model.getName() + "', which the file gives as " + name);
            return;
        }

        importFrequencyPolicy(model, referenced);
    }

    private void setNucleotideModel(PartitionSubstitutionModel model, NucModelType type) {
        model.setNucSubstitutionModel(type);
        report.recovered("the substitution model of '" + model.getName() + "', " + type);
    }

    private void setAminoAcidModel(PartitionSubstitutionModel model, String type) {
        if (type == null) {
            report.notRecovered("which amino acid model '" + model.getName() + "' uses");
            return;
        }
        for (AminoAcidModelType candidate : AminoAcidModelType.values()) {
            if (candidate.getXMLName().equalsIgnoreCase(type)) {
                model.setAaSubstitutionModel(candidate);
                report.recovered("the substitution model of '" + model.getName() + "', " + candidate);
                return;
            }
        }
        report.notRecovered("the amino acid model '" + type + "' of '" + model.getName() + "', which BEAUti does not offer");
    }

    /**
     * Empirical frequencies are shown by the alignment being referred to inside the frequency model.
     * Estimated and all-equal frequencies are written identically; the difference is that estimated ones
     * are sampled, so they have an operator.
     */
    private void importFrequencyPolicy(PartitionSubstitutionModel model, Element substitutionModel) {
        Element frequencies = firstChild(substitutionModel, FrequencyModelParser.FREQUENCIES);
        Element frequencyModel = frequencies == null ? null : firstChild(frequencies, FrequencyModelParser.FREQUENCY_MODEL);
        if (frequencyModel == null) {
            return;
        }

        if (firstChild(frequencyModel, AlignmentParser.ALIGNMENT) != null
                || firstChild(frequencyModel, SitePatternsParser.PATTERNS) != null) {
            model.setFrequencyPolicy(FrequencyPolicyType.EMPIRICAL);
            report.recovered("empirical base frequencies for '" + model.getName() + "'");
            return;
        }

        String frequenciesId = parameterIdIn(frequencyModel);
        if (frequenciesId != null && !hasOperatorOn(frequenciesId)) {
            model.setFrequencyPolicy(FrequencyPolicyType.ALLEQUAL);
            report.recovered("fixed, equal base frequencies for '" + model.getName() + "'");
        } else {
            model.setFrequencyPolicy(FrequencyPolicyType.ESTIMATED);
        }
    }

    /**
     * Reads the gamma categories and the proportion of invariant sites.
     */
    private void importSiteHeterogeneity(PartitionSubstitutionModel model, Element siteModel) {
        Element gammaShape = firstChild(siteModel, SiteModelParser.GAMMA_SHAPE);
        model.setGammaHetero(gammaShape != null);
        if (gammaShape != null) {
            Integer categories = integerAttribute(gammaShape, SiteModelParser.GAMMA_CATEGORIES);
            if (categories != null) {
                model.setRateCategories(categories);
            }
            model.setGammaHeteroEqualWeights("equal".equalsIgnoreCase(
                    gammaShape.getAttributeValue(SiteModelParser.DISCRETIZATION)));
            report.recovered("gamma rate heterogeneity for '" + model.getName() + "', " +
                    model.getRateCategories() + " categories");
        }

        boolean invariant = firstChild(siteModel, SiteModelParser.PROPORTION_INVARIANT) != null;
        model.setInvarHetero(invariant);
        if (invariant) {
            report.recovered("invariant sites for '" + model.getName() + "'");
        }
    }

    /**
     * @return the id of the first parameter defined anywhere under this element
     */
    private String parameterIdIn(Element element) {
        if (element.getName().equalsIgnoreCase(ParameterParser.PARAMETER)) {
            String id = element.getAttributeValue(XMLParser.ID);
            if (id != null) {
                return id;
            }
        }
        for (Object child : element.getChildren()) {
            String id = parameterIdIn((Element) child);
            if (id != null) {
                return id;
            }
        }
        return null;
    }

    /**
     * @return true if any operator in the file acts on the parameter of this id
     */
    private boolean hasOperatorOn(String parameterId) {
        Element operators = firstElement(SimpleOperatorScheduleParser.OPERATOR_SCHEDULE);
        if (operators == null) {
            return false;
        }
        for (Object child : operators.getChildren()) {
            if (refersTo((Element) child, parameterId)) {
                return true;
            }
        }
        return false;
    }

    private boolean refersTo(Element element, String id) {
        if (id.equals(element.getAttributeValue(XMLParser.IDREF))) {
            return true;
        }
        for (Object child : element.getChildren()) {
            if (refersTo((Element) child, id)) {
                return true;
            }
        }
        return false;
    }

    // ---- the partitions: their names and how the sites are divided among them ---------------------

    /**
     * Each partition was created from one alignment and is named after it; its real name, and whether it
     * is split into codon positions, are carried by the ids of the pattern lists that refer to that
     * alignment.
     */
    private void importPartitions() {
        for (AbstractPartitionData partition : new ArrayList<AbstractPartitionData>(options.dataPartitions)) {
            List<Element> patterns = patternListsFor(partition.getName());
            if (patterns.isEmpty()) {
                report.notRecovered("the name of partition '" + partition.getName() + "'");
                continue;
            }

            String codonHeteroPattern = codonHeteroPattern(patterns);
            if (codonHeteroPattern != null) {
                partition.getPartitionSubstitutionModel().setCodonHeteroPattern(codonHeteroPattern);
                report.recovered("codon positions " + codonHeteroPattern + " for partition '" +
                        partitionName(patterns.get(0)) + "'");
            }

            String name = partitionName(patterns.get(0));
            if (name != null && !name.equals(partition.getName()) && !options.hasPartitionData(name)) {
                options.renamePartition(partition, name);
            }
        }
    }

    /**
     * @return the pattern lists that take their data from the alignment of this name
     */
    private List<Element> patternListsFor(String alignmentId) {
        List<Element> found = new ArrayList<Element>();
        for (Element patterns : allElements(SitePatternsParser.PATTERNS, MergePatternsParser.MERGE_PATTERNS)) {
            if (patterns.getAttributeValue(XMLParser.ID) == null) {
                continue; // the halves of a merged pair, which are reached through their parent
            }
            if (alignmentId.equals(alignmentIdOf(patterns))) {
                found.add(patterns);
            }
        }
        return found;
    }

    /**
     * @return the id of the alignment a pattern list takes its data from, looking inside a merged pair
     */
    private String alignmentIdOf(Element patterns) {
        Element alignment = firstChild(patterns, AlignmentParser.ALIGNMENT);
        if (alignment != null) {
            return alignment.getAttributeValue(XMLParser.IDREF);
        }
        for (Element child : children(patterns, SitePatternsParser.PATTERNS)) {
            String id = alignmentIdOf(child);
            if (id != null) {
                return id;
            }
        }
        return null;
    }

    /**
     * @return "123" or "112" if these pattern lists divide the partition into codon positions, else null
     */
    private String codonHeteroPattern(List<Element> patterns) {
        boolean merged = false;
        int positions = 0;
        for (Element element : patterns) {
            String prefix = codonPrefix(element.getAttributeValue(XMLParser.ID));
            if (prefix == null) {
                continue;
            }
            positions++;
            if (prefix.equals(CODON_1_AND_2)) {
                merged = true;
            }
        }
        if (positions == 0) {
            return null;
        }
        return merged ? "112" : "123";
    }

    /**
     * @return the name of the partition a pattern list belongs to, or null when the file does not say
     * (BEAUti leaves the prefix off when there is only one partition)
     */
    private String partitionName(Element patterns) {
        String id = patterns.getAttributeValue(XMLParser.ID);
        if (id == null || !id.endsWith(SitePatternsParser.PATTERNS)) {
            return null;
        }
        String prefix = id.substring(0, id.length() - SitePatternsParser.PATTERNS.length());
        String codonPrefix = codonPrefix(id);
        if (codonPrefix != null) {
            prefix = prefix.substring(0, prefix.length() - codonPrefix.length());
        }
        return prefix.isEmpty() ? null : prefix.substring(0, prefix.length() - 1);
    }

    /**
     * @return the codon position part of an id ("CP1.", "CP1+2." ...), or null if it has none
     */
    private static String codonPrefix(String id) {
        if (id == null) {
            return null;
        }
        String withoutName = id.substring(0, id.length() - SitePatternsParser.PATTERNS.length());
        for (String prefix : new String[]{CODON_1_AND_2, "CP1.", "CP2.", "CP3."}) {
            if (withoutName.endsWith(prefix)) {
                return prefix;
            }
        }
        return null;
    }

    // ---- which partitions share a substitution, clock or tree model -------------------------------

    /**
     * Reads from the tree likelihoods which partitions share a model, and gives each group its own model,
     * named as the file names it. Everything starts out sharing one model of each kind.
     */
    private void importLinkage() {
        Map<AbstractPartitionData, String> substitutionModelIds = new LinkedHashMap<AbstractPartitionData, String>();
        Map<AbstractPartitionData, String> treeModelIds = new LinkedHashMap<AbstractPartitionData, String>();
        Map<AbstractPartitionData, String> clockModelIds = new LinkedHashMap<AbstractPartitionData, String>();

        for (Element likelihood : allElements(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD,
                TreeLikelihoodParser.TREE_LIKELIHOOD, TreeLikelihoodParser.ANCESTRAL_TREE_LIKELIHOOD,
                MarkovJumpsTreeLikelihoodParser.MARKOV_JUMP_TREE_LIKELIHOOD)) {
            String treeModelId = idrefOf(likelihood, DefaultTreeModel.TREE_MODEL);
            String clockModelId = branchRatesIdOf(likelihood);

            List<Element> partitionElements = children(likelihood, TreeDataLikelihoodParser.PARTITION);
            if (partitionElements.isEmpty()) {
                // the older form: the likelihood is itself for one partition
                partitionElements = java.util.Collections.singletonList(likelihood);
            }
            for (Element partitionElement : partitionElements) {
                AbstractPartitionData partition = partitionOf(partitionElement);
                if (partition == null) {
                    continue;
                }
                String siteModelId = idrefOf(partitionElement, SiteModelParser.SITE_MODEL);
                if (siteModelId != null) {
                    substitutionModelIds.put(partition, siteModelId);
                }
                if (treeModelId != null) {
                    treeModelIds.put(partition, treeModelId);
                }
                if (clockModelId != null) {
                    clockModelIds.put(partition, clockModelId);
                }
            }
        }

        splitSubstitutionModels(substitutionModelIds);
        splitTreeModels(treeModelIds);
        splitClockModels(clockModelIds);
    }

    private void splitSubstitutionModels(Map<AbstractPartitionData, String> modelIds) {
        for (Map.Entry<String, List<AbstractPartitionData>> group : groupByModel(modelIds).entrySet()) {
            String name = modelName(group.getKey(), SiteModelParser.SITE_MODEL);
            List<AbstractPartitionData> partitions = group.getValue();
            PartitionSubstitutionModel model = partitions.get(0).getPartitionSubstitutionModel();
            if (name != null && !name.equals(model.getName())) {
                model = new PartitionSubstitutionModel(options, name, model);
            }
            for (AbstractPartitionData partition : partitions) {
                partition.setPartitionSubstitutionModel(model);
            }
        }
        options.clearDataPartitionCaches();
    }

    private void splitTreeModels(Map<AbstractPartitionData, String> modelIds) {
        for (Map.Entry<String, List<AbstractPartitionData>> group : groupByModel(modelIds).entrySet()) {
            String name = modelName(group.getKey(), DefaultTreeModel.TREE_MODEL);
            List<AbstractPartitionData> partitions = group.getValue();
            PartitionTreeModel treeModel = partitions.get(0).getPartitionTreeModel();
            if (name != null && !name.equals(treeModel.getName())) {
                treeModel = new PartitionTreeModel(options, name, treeModel);
            }
            for (AbstractPartitionData partition : partitions) {
                partition.setPartitionTreeModel(treeModel);
                PartitionClockModel clockModel = partition.getPartitionClockModel();
                if (clockModel != null) {
                    clockModel.setPartitionTreeModel(treeModel);
                }
            }
        }
        options.clearDataPartitionCaches();
    }

    private void splitClockModels(Map<AbstractPartitionData, String> modelIds) {
        for (Map.Entry<String, List<AbstractPartitionData>> group : groupByModel(modelIds).entrySet()) {
            String name = modelName(group.getKey(), BRANCH_RATES);
            List<AbstractPartitionData> partitions = group.getValue();
            PartitionClockModel clockModel = partitions.get(0).getPartitionClockModel();
            if (clockModel != null && name != null && !name.equals(clockModel.getName())) {
                clockModel = new PartitionClockModel(options, name, clockModel);
            }
            for (AbstractPartitionData partition : partitions) {
                if (clockModel != null) {
                    partition.setPartitionClockModel(clockModel);
                }
            }
        }
        options.clearDataPartitionCaches();
    }

    private Map<String, List<AbstractPartitionData>> groupByModel(Map<AbstractPartitionData, String> modelIds) {
        Map<String, List<AbstractPartitionData>> groups = new LinkedHashMap<String, List<AbstractPartitionData>>();
        for (Map.Entry<AbstractPartitionData, String> entry : modelIds.entrySet()) {
            List<AbstractPartitionData> group = groups.get(entry.getValue());
            if (group == null) {
                group = new ArrayList<AbstractPartitionData>();
                groups.put(entry.getValue(), group);
            }
            group.add(entry.getKey());
        }
        return groups;
    }

    /**
     * @return the name of the model an id refers to, or null when the file leaves the prefix off because
     * there is only one of that kind
     */
    private static String modelName(String id, String baseName) {
        if (id == null || !id.endsWith(baseName)) {
            return null;
        }
        String prefix = id.substring(0, id.length() - baseName.length());
        // a site model of a codon-partitioned model carries the position as well, which is not its name
        for (String codonPrefix : new String[]{CODON_1_AND_2, "CP1.", "CP2.", "CP3."}) {
            if (prefix.endsWith(codonPrefix)) {
                prefix = prefix.substring(0, prefix.length() - codonPrefix.length());
            }
        }
        return prefix.isEmpty() ? null : prefix.substring(0, prefix.length() - 1);
    }

    /**
     * @return the partition whose pattern list this element refers to
     */
    private AbstractPartitionData partitionOf(Element element) {
        String patternsId = idrefOf(element, SitePatternsParser.PATTERNS);
        if (patternsId == null) {
            patternsId = idrefOf(element, MergePatternsParser.MERGE_PATTERNS);
        }
        if (patternsId == null) {
            return null;
        }
        Element patterns = elementById(patternsId);
        String alignmentId = patterns == null ? null : alignmentIdOf(patterns);
        for (AbstractPartitionData partition : options.dataPartitions) {
            if (partition.getName().equals(alignmentId) || partition.getName().equals(partitionName(patterns))) {
                return partition;
            }
        }
        return null;
    }

    /**
     * @return the id of the branch rates model a tree likelihood uses, whatever kind of clock it is
     */
    private String branchRatesIdOf(Element likelihood) {
        for (Object child : likelihood.getChildren()) {
            Element element = (Element) child;
            String idref = element.getAttributeValue(XMLParser.IDREF);
            if (idref != null && idref.endsWith(BRANCH_RATES)) {
                return idref;
            }
        }
        return null;
    }

    private static String idrefOf(Element element, String childName) {
        Element child = firstChild(element, childName);
        return child == null ? null : child.getAttributeValue(XMLParser.IDREF);
    }

    // ---- the MCMC block, the logs and the file names ----------------------------------------------

    private void importMCMCSettings() {
        Element mcmc = firstElement(MCMCParser.MCMC);
        if (mcmc == null) {
            report.notRecovered("the MCMC settings (the file has no " + MCMCParser.MCMC + " element)");
            return;
        }

        Integer chainLength = integerAttribute(mcmc, MCMCParser.CHAIN_LENGTH);
        if (chainLength != null) {
            options.chainLength = chainLength;
            report.recovered("the chain length, " + chainLength);
        }

        Boolean autoOptimize = booleanAttribute(mcmc, MCMCParser.AUTO_OPTIMIZE);
        if (autoOptimize != null) {
            options.autoOptimize = autoOptimize;
        }

        // the operator analysis is written as the file it goes to, and only when it is on
        String operatorAnalysisFileName = mcmc.getAttributeValue(MCMCParser.OPERATOR_ANALYSIS);
        options.operatorAnalysis = operatorAnalysisFileName != null;
        if (operatorAnalysisFileName != null) {
            options.operatorAnalysisFileName = operatorAnalysisFileName;
        }

        importLogSettings(mcmc);
    }

    private void importLogSettings(Element mcmc) {
        String stem = null;

        for (Element log : children(mcmc, LoggerParser.LOG)) {
            String id = log.getAttributeValue(XMLParser.ID);
            Integer logEvery = integerAttribute(log, LoggerParser.LOG_EVERY);
            if (logEvery == null) {
                continue;
            }
            if ("screenLog".equals(id)) {
                options.echoEvery = logEvery;
            } else {
                // the log to a file: its name gives the stem every other file name is built from
                options.logEvery = logEvery;
                String fileName = log.getAttributeValue(LoggerParser.FILE_NAME);
                if (fileName != null) {
                    options.logFileName = fileName;
                    stem = trimExtension(fileName, ".log");
                }
                Boolean overwrite = booleanAttribute(log, LoggerParser.ALLOW_OVERWRITE_LOG);
                if (overwrite != null) {
                    options.allowOverwriteLog = overwrite;
                }
                report.recovered("the log interval, " + logEvery);
            }
        }

        // a second tree log per tree, with branch lengths in substitutions
        options.substTreeLog = false;
        for (Element logTree : children(mcmc, "logTree")) {
            if ("substitutions".equals(logTree.getAttributeValue("branchLengths"))) {
                options.substTreeLog = true;
            }
        }

        Element checkpoint = firstChild(mcmc, CheckpointLoggerParser.LOG_CHECKPOINT);
        if (checkpoint != null) {
            Integer checkpointEvery = integerAttribute(checkpoint, CheckpointLoggerParser.CHECKPOINT_EVERY);
            if (checkpointEvery != null) {
                options.checkpointEvery = checkpointEvery;
            }
        }

        if (stem != null) {
            options.fileNameStem = stem;
            options.updateFileNames(stem.endsWith(".txt"));
            report.recovered("the file name stem, " + stem);
        }
    }

    // ---- reading the document ---------------------------------------------------------------------

    private void indexElements(Element element) {
        String id = element.getAttributeValue(XMLParser.ID);
        if (id != null) {
            elementsById.put(id, element);
        }
        for (Object child : element.getChildren()) {
            indexElements((Element) child);
        }
    }

    /**
     * @return the element with this id, or null if the file has none
     */
    protected Element elementById(String id) {
        return elementsById.get(id);
    }

    /**
     * @return the element an idref attribute points at, or null
     */
    protected Element referenced(Element element, String attributeName) {
        String idref = element.getAttributeValue(attributeName);
        return idref == null ? null : elementsById.get(idref);
    }

    /**
     * @return the first element of this name anywhere in the document, or null
     */
    protected Element firstElement(String name) {
        return firstElement(root, name);
    }

    /**
     * @return every element of the document with one of these names, in document order
     */
    protected List<Element> allElements(String... names) {
        List<Element> found = new ArrayList<Element>();
        collectElements(root, names, found);
        return found;
    }

    private static void collectElements(Element element, String[] names, List<Element> found) {
        for (String name : names) {
            if (element.getName().equalsIgnoreCase(name)) {
                found.add(element);
                break;
            }
        }
        for (Object child : element.getChildren()) {
            collectElements((Element) child, names, found);
        }
    }

    private static Element firstElement(Element element, String name) {
        if (element.getName().equalsIgnoreCase(name)) {
            return element;
        }
        for (Object child : element.getChildren()) {
            Element found = firstElement((Element) child, name);
            if (found != null) {
                return found;
            }
        }
        return null;
    }

    protected static Element firstChild(Element element, String name) {
        List<Element> found = children(element, name);
        return found.isEmpty() ? null : found.get(0);
    }

    @SuppressWarnings("unchecked")
    protected static List<Element> children(Element element, String name) {
        List<Element> found = new ArrayList<Element>();
        for (Object child : element.getChildren()) {
            Element childElement = (Element) child;
            if (childElement.getName().equalsIgnoreCase(name)) {
                found.add(childElement);
            }
        }
        return found;
    }

    protected static Integer integerAttribute(Element element, String name) {
        String value = element.getAttributeValue(name);
        try {
            return value == null ? null : Integer.valueOf(value.trim());
        } catch (NumberFormatException e) {
            return null;
        }
    }

    protected static Double doubleAttribute(Element element, String name) {
        String value = element.getAttributeValue(name);
        try {
            return value == null ? null : Double.valueOf(value.trim());
        } catch (NumberFormatException e) {
            return null;
        }
    }

    protected static Boolean booleanAttribute(Element element, String name) {
        String value = element.getAttributeValue(name);
        if (value == null) {
            return null;
        }
        return Boolean.valueOf(value.trim().equalsIgnoreCase("true"));
    }

    private static String trimExtension(String fileName, String extension) {
        String name = fileName;
        if (name.endsWith(".txt")) {
            name = name.substring(0, name.length() - ".txt".length());
        }
        if (name.endsWith(extension)) {
            name = name.substring(0, name.length() - extension.length());
        }
        return name;
    }
}
