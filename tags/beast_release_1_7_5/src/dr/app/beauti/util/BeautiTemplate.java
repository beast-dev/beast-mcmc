/*
 * BeautiTemplate.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
package dr.app.beauti.util;

import dr.app.beast.BeastVersion;
import dr.app.beauti.types.PriorType;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.ModelOptions;
import dr.app.beauti.options.Operator;
import dr.app.beauti.options.Parameter;
import dr.evolution.util.Date;
import dr.evolution.util.Taxa;
import dr.evolution.util.Units;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.xml.XMLParser;
import org.jdom.Document;
import org.jdom.Element;


/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @author Walter Xie
 * @version $Id: BeautiTemplate.java, rambaut Exp $
 */
public class BeautiTemplate extends ModelOptions {

	private final BeautiOptions options;

    public BeautiTemplate(BeautiOptions options) {
    	this.options = options;

    }


    /**
     * Write options from a file
     *
     * @param guessDates guess dates?
     * @return the Document
     */
    public Document create(boolean guessDates) {

        final BeastVersion version = new BeastVersion();
        Element root = new Element("beauti");
        root.setAttribute("version", version.getVersion());

        Element dataElement = new Element("data");

//        dataElement.addContent(createChild("fileNameStem", fileNameStem));

//        dataElement.addContent(createChild("datesUnits", options.datesUnits));
//        dataElement.addContent(createChild("datesDirection", options.datesDirection));
        dataElement.addContent(createChild("translation", options.translation));
        //TODO:
//        dataElement.addContent(createChild("startingTreeType", startingTreeType.name()));

        dataElement.addContent(createChild("guessDates", guessDates));
        dataElement.addContent(createChild("guessType", options.dateGuesser.guessType.name()));
        dataElement.addContent(createChild("fromLast", options.dateGuesser.fromLast));
        dataElement.addContent(createChild("order", options.dateGuesser.order));
        dataElement.addContent(createChild("prefix", options.dateGuesser.prefix));
        dataElement.addContent(createChild("offset", options.dateGuesser.offset));
        dataElement.addContent(createChild("unlessLessThan", options.dateGuesser.unlessLessThan));
        dataElement.addContent(createChild("offset2", options.dateGuesser.offset2));

        root.addContent(dataElement);

        Element taxaElement = new Element(TaxaParser.TAXA);

        for (Taxa taxonSet : options.taxonSets) {
            Element taxonSetElement = new Element("taxonSet");
            taxonSetElement.addContent(createChild(XMLParser.ID, taxonSet.getId()));
            taxonSetElement.addContent(createChild("enforceMonophyly",
            		options.taxonSetsMono.get(taxonSet) ? "true" : "false"));
            for (int j = 0; j < taxonSet.getTaxonCount(); j++) {
                Element taxonElement = new Element(TaxonParser.TAXON);
                taxonElement.addContent(createChild(XMLParser.ID, taxonSet.getTaxon(j).getId()));
                taxonSetElement.addContent(taxonElement);
            }
            taxaElement.addContent(taxonSetElement);
        }

        root.addContent(taxaElement);

//        for (PartitionSubstitutionModel model : partitionModels) {
//
//            Element modelElement = new Element("model");
//
//            /*modelElement.addContent(createChild("nucSubstitutionModel", nucSubstitutionModel));
//                           modelElement.addContent(createChild("aaSubstitutionModel", aaSubstitutionModel));
//                           modelElement.addContent(createChild("binarySubstitutionModel", binarySubstitutionModel));
//                           modelElement.addContent(createChild("frequencyPolicy", frequencyPolicy));
//                           modelElement.addContent(createChild("gammaHetero", gammaHetero));
//                           modelElement.addContent(createChild("gammaCategories", gammaCategories));
//                           modelElement.addContent(createChild("invarHetero", invarHetero));
//                           modelElement.addContent(createChild("codonHeteroPattern", codonHeteroPattern));
//                           modelElement.addContent(createChild("maximumTipHeight", maximumTipHeight));
//                           modelElement.addContent(createChild("hasSetFixedSubstitutionRate", hasSetFixedSubstitutionRate));
//                           modelElement.addContent(createChild("meanSubstitutionRate", meanSubstitutionRate));
//                           modelElement.addContent(createChild("fixedSubstitutionRate", fixedSubstitutionRate));
//                           modelElement.addContent(createChild("unlinkedSubstitutionModel", unlinkedSubstitutionModel));
//                           modelElement.addContent(createChild("unlinkedHeterogeneityModel", unlinkedHeterogeneityModel));
//                           modelElement.addContent(createChild("unlinkedFrequencyModel", unlinkedFrequencyModel));
//                           modelElement.addContent(createChild("clockModel", clockModel));
//                           modelElement.addContent(createChild("nodeHeightPrior", nodeHeightPrior));
//                           modelElement.addContent(createChild("parameterization", parameterization));
//                           modelElement.addContent(createChild("skylineGroupCount", skylineGroupCount));
//                           modelElement.addContent(createChild("skylineModel", skylineModel));
//                           modelElement.addContent(createChild("fixedTree", fixedTree)); */
//
//            root.addContent(modelElement);
//        }

        Element priorsElement = new Element("priors");

        for (String name : getParameters().keySet()) {
            Parameter parameter = getParameters().get(name);
            Element e = new Element(name);
            e.addContent(createChild("initial", parameter.initial));
            e.addContent(createChild("priorType", parameter.priorType));
            e.addContent(createChild("priorEdited", parameter.isPriorEdited()));
//            e.addContent(createChild("uniformLower", parameter.uniformLower));
//            e.addContent(createChild("uniformUpper", parameter.uniformUpper));
//            e.addContent(createChild("exponentialMean", parameter.exponentialMean));
//            e.addContent(createChild("exponentialOffset", parameter.exponentialOffset));
//            e.addContent(createChild("normalMean", parameter.normalMean));
//            e.addContent(createChild("normalStdev", parameter.normalStdev));
//            e.addContent(createChild("logNormalMean", parameter.logNormalMean));
//            e.addContent(createChild("logNormalStdev", parameter.logNormalStdev));
//            e.addContent(createChild("logNormalOffset", parameter.logNormalOffset));
//            e.addContent(createChild("gammaAlpha", parameter.gammaAlpha));
//            e.addContent(createChild("gammaBeta", parameter.gammaBeta));
//            e.addContent(createChild("gammaOffset", parameter.gammaOffset));
            priorsElement.addContent(e);
        }

        for (Taxa taxonSet : options.taxonSets) {
            Parameter statistic = getStatistics().get(taxonSet);
            Element e = new Element(statistic.getXMLName());
            e.addContent(createChild("initial", statistic.initial));
            e.addContent(createChild("priorType", statistic.priorType));
            e.addContent(createChild("priorEdited", statistic.isPriorEdited()));
//            e.addContent(createChild("uniformLower", statistic.uniformLower));
//            e.addContent(createChild("uniformUpper", statistic.uniformUpper));
//            e.addContent(createChild("exponentialMean", statistic.exponentialMean));
//            e.addContent(createChild("exponentialOffset", statistic.exponentialOffset));
//            e.addContent(createChild("normalMean", statistic.normalMean));
//            e.addContent(createChild("normalStdev", statistic.normalStdev));
//            e.addContent(createChild("logNormalMean", statistic.logNormalMean));
//            e.addContent(createChild("logNormalStdev", statistic.logNormalStdev));
//            e.addContent(createChild("logNormalOffset", statistic.logNormalOffset));
//            e.addContent(createChild("gammaAlpha", statistic.gammaAlpha));
//            e.addContent(createChild("gammaBeta", statistic.gammaBeta));
//            e.addContent(createChild("gammaOffset", statistic.gammaOffset));
            priorsElement.addContent(e);
        }

        root.addContent(priorsElement);

        Element operatorsElement = new Element("operators");

        operatorsElement.addContent(createChild("autoOptimize", options.autoOptimize));
        for (String name : getOperators().keySet()) {
            Operator operator = getOperators().get(name);
            Element e = new Element(name);
            e.addContent(createChild("tuning", operator.tuning));
            e.addContent(createChild("tuningEdited", operator.tuningEdited));
            e.addContent(createChild("weight", operator.weight));
            e.addContent(createChild("inUse", operator.inUse));
            operatorsElement.addContent(e);
        }

        root.addContent(operatorsElement);

        Element mcmcElement = new Element("mcmc");

        mcmcElement.addContent(createChild("chainLength", options.chainLength));
        mcmcElement.addContent(createChild("logEvery", options.logEvery));
        mcmcElement.addContent(createChild("echoEvery", options.echoEvery));
        //if (logFileName != null) mcmcElement.addContent(createChild("logFileName", logFileName));
        //if (treeFileName != null) mcmcElement.addContent(createChild("treeFileName", treeFileName));
        //mcmcElement.addContent(createChild("mapTreeLog", mapTreeLog));
        //if (mapTreeFileName != null) mcmcElement.addContent(createChild("mapTreeFileName", mapTreeFileName));
        mcmcElement.addContent(createChild("substTreeLog", options.substTreeLog));
        //if (substTreeFileName != null) mcmcElement.addContent(createChild("substTreeFileName", substTreeFileName));

        root.addContent(mcmcElement);

        return new Document(root);
    }

    private Element createChild(String name, String value) {
        Element e = new Element(name);
        if (value != null) {
            e.setText(value);
        }
        return e;
    }

    private Element createChild(String name, int value) {
        Element e = new Element(name);
        e.setText(Integer.toString(value));
        return e;
    }

    private Element createChild(String name, PriorType value) {
        Element e = new Element(name);
        e.setText(value.name());
        return e;
    }

    private Element createChild(String name, double value) {
        Element e = new Element(name);
        e.setText(Double.toString(value));
        return e;
    }

    private Element createChild(String name, boolean value) {
        Element e = new Element(name);
        e.setText(value ? "true" : "false");
        return e;
    }

    /**
     * Read options from a file
     *
     * @param document the Document
     * @throws dr.xml.XMLParseException if there is a problem with XML parsing
     */
    public void parse(Document document) throws dr.xml.XMLParseException {

        Element root = document.getRootElement();
        if (!root.getName().equals("beauti")) {
            throw new dr.xml.XMLParseException("This document does not appear to be a BEAUti file");
        }

        Element taxaElement = root.getChild(TaxaParser.TAXA);
        Element modelElement = root.getChild("model");
        Element priorsElement = root.getChild("priors");
        Element operatorsElement = root.getChild("operators");
        Element mcmcElement = root.getChild("mcmc");
        /*
                  if (taxaElement != null) {
                      for (Object ts : taxaElement.getChildren("taxonSet")) {
                          Element taxonSetElement = (Element) ts;

                          String id = getStringChild(taxonSetElement, XMLParser.ID, "");
                          final Taxa taxonSet = new Taxa(id);

                          Boolean enforceMonophyly = Boolean.valueOf(getStringChild(taxonSetElement, "enforceMonophyly", "false"));
                          for (Object o : taxonSetElement.getChildren("taxon")) {
                              Element taxonElement = (Element) o;
                              String taxonId = getStringChild(taxonElement, XMLParser.ID, "");
                              int index = taxonList.getTaxonIndex(taxonId);
                              if (index != -1) {
                                  taxonSet.addTaxon(taxonList.getTaxon(index));
                              }
                          }
                          taxonSets.add(taxonSet);
                          taxonSetsMono.put(taxonSet, enforceMonophyly);
                      }
                  }

                  if (modelElement != null) {
                      nucSubstitutionModel = getIntegerChild(modelElement, "nucSubstitutionModel", HKY);
                      aaSubstitutionModel = getIntegerChild(modelElement, "aaSubstitutionModel", BLOSUM_62);
                      binarySubstitutionModel = getIntegerChild(modelElement, "binarySubstitutionModel", BIN_SIMPLE);
                      frequencyPolicy = getIntegerChild(modelElement, "frequencyPolicy", ESTIMATED);
                      gammaHetero = getBooleanChild(modelElement, "gammaHetero", false);
                      gammaCategories = getIntegerChild(modelElement, "gammaCategories", 5);
                      invarHetero = getBooleanChild(modelElement, "invarHetero", false);
                      codonHeteroPattern = (getBooleanChild(modelElement, "codonHetero", false) ? "123" : null);
                      codonHeteroPattern = getStringChild(modelElement, "codonHeteroPattern", null);
                      maximumTipHeight = getDoubleChild(modelElement, "maximumTipHeight", 0.0);
                      fixedSubstitutionRate = getBooleanChild(modelElement, "fixedSubstitutionRate", false);
                      hasSetFixedSubstitutionRate = getBooleanChild(modelElement, "hasSetFixedSubstitutionRate", false);
                      meanSubstitutionRate = getDoubleChild(modelElement, "meanSubstitutionRate", 1.0);
                      unlinkedSubstitutionModel = getBooleanChild(modelElement, "unlinkedSubstitutionModel", false);
                      unlinkedHeterogeneityModel = getBooleanChild(modelElement, "unlinkedHeterogeneityModel", false);
                      unlinkedFrequencyModel = getBooleanChild(modelElement, "unlinkedFrequencyModel", false);

                      clockModel = getIntegerChild(modelElement, "clockModel", clockModel);

                      // the old name was "coalescentModel" so try to read this first
                      nodeHeightPrior = getIntegerChild(modelElement, "coalescentModel", CONSTANT);
                      nodeHeightPrior = getIntegerChild(modelElement, "nodeHeightPrior", nodeHeightPrior);
                      // we don't allow no nodeHeightPrior in BEAUti so switch it to Yule:
                      if (nodeHeightPrior == NONE_TREE_PRIOR) nodeHeightPrior = YULE;

                      parameterization = getIntegerChild(modelElement, "parameterization", GROWTH_RATE);
                      skylineGroupCount = getIntegerChild(modelElement, "skylineGroupCount", 10);
                      skylineModel = getIntegerChild(modelElement, "skylineModel", CONSTANT_SKYLINE);
                      fixedTree = getBooleanChild(modelElement, "fixedTree", false);
                  }

                  if (operatorsElement != null) {
                      autoOptimize = getBooleanChild(operatorsElement, "autoOptimize", true);
                      for (String name : operators.keySet()) {
                          Operator operator = operators.get(name);
                          Element e = operatorsElement.getChild(name);
                          if (e == null) {
                              throw new XMLParseException("Operators element, " + name + " missing");
                          }

                          operator.tuning = getDoubleChild(e, "tuning", 1.0);
                          operator.tuningEdited = getBooleanChild(e, "tuningEdited", false);
                          operator.weight = getDoubleChild(e, "weight", 1);
                          operator.inUse = getBooleanChild(e, "inUse", true);
                      }
                  }

                  if (priorsElement != null) {
                      for (String name : parameters.keySet()) {
                          Parameter parameter = parameters.get(name);
                          Element e = priorsElement.getChild(name);
                          if (e == null) {
                              throw new XMLParseException("Priors element, " + name + " missing");
                          }

                          parameter.initial = getDoubleChild(e, "initial", 1.0);
                          parameter.priorType = PriorType.valueOf(getStringChild(e, "priorType", PriorType.UNIFORM_PRIOR.name()));
                          parameter.priorEdited = getBooleanChild(e, "priorEdited", false);
                          parameter.uniformLower = Math.max(getDoubleChild(e, "uniformLower", parameter.uniformLower), parameter.lower);
                          parameter.uniformUpper = Math.min(getDoubleChild(e, "uniformUpper", parameter.uniformUpper), parameter.upper);
                          parameter.exponentialMean = getDoubleChild(e, "exponentialMean", parameter.exponentialMean);
                          parameter.exponentialOffset = getDoubleChild(e, "exponentialOffset", parameter.exponentialOffset);
                          parameter.normalMean = getDoubleChild(e, "normalMean", parameter.normalMean);
                          parameter.normalStdev = getDoubleChild(e, "normalStdev", parameter.normalStdev);
                          parameter.logNormalMean = getDoubleChild(e, "logNormalMean", parameter.logNormalMean);
                          parameter.logNormalStdev = getDoubleChild(e, "logNormalStdev", parameter.logNormalStdev);
                          parameter.logNormalOffset = getDoubleChild(e, "logNormalOffset", parameter.logNormalOffset);
                          parameter.gammaAlpha = getDoubleChild(e, "gammaAlpha", parameter.gammaAlpha);
                          parameter.gammaBeta = getDoubleChild(e, "gammaBeta", parameter.gammaBeta);
                          parameter.gammaOffset = getDoubleChild(e, "gammaOffset", parameter.gammaOffset);
                      }

                      for (Taxa taxonSet : taxonSets) {
                          Parameter statistic = statistics.get(taxonSet);
                          if (statistic == null) {
                              statistic = new Parameter(this, taxonSet, "tMRCA for taxon set ");
                              statistics.put(taxonSet, statistic);
                          }
                          Element e = priorsElement.getChild(statistic.getXMLName());
                          statistic.initial = getDoubleChild(e, "initial", 1.0);
                          statistic.priorType = PriorType.valueOf(getStringChild(e, "priorType", PriorType.UNIFORM_PRIOR.name()));
                          statistic.priorEdited = getBooleanChild(e, "priorEdited", false);
                          statistic.uniformLower = getDoubleChild(e, "uniformLower", statistic.uniformLower);
                          statistic.uniformUpper = getDoubleChild(e, "uniformUpper", statistic.uniformUpper);
                          statistic.exponentialMean = getDoubleChild(e, "exponentialMean", statistic.exponentialMean);
                          statistic.exponentialOffset = getDoubleChild(e, "exponentialOffset", statistic.exponentialOffset);
                          statistic.normalMean = getDoubleChild(e, "normalMean", statistic.normalMean);
                          statistic.normalStdev = getDoubleChild(e, "normalStdev", statistic.normalStdev);
                          statistic.logNormalMean = getDoubleChild(e, "logNormalMean", statistic.logNormalMean);
                          statistic.logNormalStdev = getDoubleChild(e, "logNormalStdev", statistic.logNormalStdev);
                          statistic.logNormalOffset = getDoubleChild(e, "logNormalOffset", statistic.logNormalOffset);
                          statistic.gammaAlpha = getDoubleChild(e, "gammaAlpha", statistic.gammaAlpha);
                          statistic.gammaBeta = getDoubleChild(e, "gammaBeta", statistic.gammaBeta);
                          statistic.gammaOffset = getDoubleChild(e, "gammaOffset", statistic.gammaOffset);
                      }

                  }


                  if (mcmcElement != null) {
                      upgmaStartingTree = getBooleanChild(mcmcElement, "upgmaStartingTree", true);
                      chainLength = getIntegerChild(mcmcElement, "chainLength", 100000000);
                      logEvery = getIntegerChild(mcmcElement, "logEvery", 1000);
                      echoEvery = getIntegerChild(mcmcElement, "echoEvery", 1000);
                      logFileName = getStringChild(mcmcElement, "logFileName", null);
                      treeFileName = getStringChild(mcmcElement, "treeFileName", null);
                      mapTreeLog = getBooleanChild(mcmcElement, "mapTreeLog", false);
                      mapTreeFileName = getStringChild(mcmcElement, "mapTreeFileName", null);
                      substTreeLog = getBooleanChild(mcmcElement, "substTreeLog", false);
                      substTreeFileName = getStringChild(mcmcElement, "substTreeFileName", null);
                  }      */
    }

    private String getStringChild(Element element, String childName, String defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null || value.length() == 0) return defaultValue;
        return value;
    }

    private int getIntegerChild(Element element, String childName, int defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return Integer.parseInt(value);
    }

    private double getDoubleChild(Element element, String childName, double defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return Double.parseDouble(value);
    }

    private boolean getBooleanChild(Element element, String childName, boolean defaultValue) {
        String value = element.getChildTextTrim(childName);
        if (value == null) return defaultValue;
        return value.equals("true");
    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

}
