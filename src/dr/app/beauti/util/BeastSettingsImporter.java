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
import dr.evomodel.tree.DefaultTreeModel;
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
        importMCMCSettings();
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
