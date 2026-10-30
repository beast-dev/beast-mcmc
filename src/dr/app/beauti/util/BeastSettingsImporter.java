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

import dr.app.beauti.options.BeautiOptions;
import dr.inferencexml.MCMCParser;
import dr.inferencexml.loggers.CheckpointLoggerParser;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.XMLParser;
import org.jdom.Element;

import java.util.ArrayList;
import java.util.HashMap;
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
        importMCMCSettings();
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
