/*
 * TimeTreeExportParser.java
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

package dr.evomodelxml.tree;

import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TimeTreeExporter;
import dr.evomodel.tree.TreeLogger;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.loggers.LoggerParser;
import dr.util.FileHelpers;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Logs the trees once, as MCMC logs them, and writes the tables and the run record, when parsed (after the
 * maximizer has run):
 * <pre>
 * &lt;timeTreeExport fileNameStem="run"&gt;
 *     &lt;logTree idref="timeTreeLog"/&gt;
 *     &lt;maximizeWrtParameter idref="maximizer"/&gt;
 *     &lt;treeLikelihoods&gt;&lt;likelihood idref="likelihood"/&gt;&lt;/treeLikelihoods&gt;
 *     &lt;priors&gt;&lt;prior idref="prior"/&gt;&lt;/priors&gt;
 *     &lt;timeTree&gt;
 *         &lt;treeModel idref="treeModel"/&gt;
 *         &lt;strictClockBranchRates idref="branchRates"/&gt;
 *         &lt;taxa idref="ingroup"/&gt;
 *     &lt;/timeTree&gt;
 *     &lt;notOptimized reason="fixed in the Priors tab"&gt;&lt;parameter idref="kappa"/&gt;&lt;/notOptimized&gt;
 * &lt;/timeTreeExport&gt;
 * </pre>
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeExportParser extends AbstractXMLObjectParser {

    public static final String TIME_TREE_EXPORT = "timeTreeExport";
    public static final String FILE_NAME_STEM = "fileNameStem";
    public static final String TREE_LIKELIHOODS = "treeLikelihoods";
    public static final String PRIORS = "priors";
    public static final String TIME_TREE = "timeTree";
    public static final String NOT_OPTIMIZED = "notOptimized";
    public static final String REASON = "reason";
    public static final String CAN_OPTIMIZE = "canOptimize";
    public static final String NOT_IDENTIFIED = "notIdentified";

    public String getParserName() {
        return TIME_TREE_EXPORT;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        File stem = FileHelpers.getFile(xo.getStringAttribute(FILE_NAME_STEM), System.getProperty("file.name.prefix"));

        MaximizerWrtParameter maximizer = (MaximizerWrtParameter) xo.getChild(MaximizerWrtParameter.class);
        Likelihood likelihood = xo.hasChildNamed(TREE_LIKELIHOODS) ?
                (Likelihood) xo.getChild(TREE_LIKELIHOODS).getChild(Likelihood.class) : null;
        Likelihood prior = xo.hasChildNamed(PRIORS) ? (Likelihood) xo.getChild(PRIORS).getChild(Likelihood.class) : null;

        List<TimeTreeExporter.ExportedTree> trees = new ArrayList<TimeTreeExporter.ExportedTree>();
        List<TimeTreeExporter.FixedParameters> fixed = new ArrayList<TimeTreeExporter.FixedParameters>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (!(xo.getChild(i) instanceof XMLObject)) {
                continue;
            }
            XMLObject child = (XMLObject) xo.getChild(i);
            if (child.getName().equals(TIME_TREE)) {
                Tree tree = (Tree) child.getChild(Tree.class);
                List<TaxonList> taxonSets = new ArrayList<TaxonList>();
                for (int j = 0; j < child.getChildCount(); ++j) {
                    // a tree is a taxon list too
                    if (child.getChild(j) instanceof TaxonList && !(child.getChild(j) instanceof Tree)) {
                        taxonSets.add((TaxonList) child.getChild(j));
                    }
                }
                trees.add(new TimeTreeExporter.ExportedTree(tree,
                        (BranchRateModel) child.getChild(BranchRateModel.class), taxonSets));
            } else if (child.getName().equals(NOT_OPTIMIZED)) {
                List<Parameter> parameters = new ArrayList<Parameter>();
                for (int j = 0; j < child.getChildCount(); ++j) {
                    parameters.add((Parameter) child.getChild(j));
                }
                fixed.add(new TimeTreeExporter.FixedParameters(child.getStringAttribute(REASON), parameters,
                        child.getAttribute(CAN_OPTIMIZE, true)));
            }
        }

        List<TreeLogger> treeLoggers = new ArrayList<TreeLogger>();
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (xo.getChild(i) instanceof TreeLogger) {
                treeLoggers.add((TreeLogger) xo.getChild(i));
            }
        }

        TimeTreeExporter exporter = new TimeTreeExporter(stem, trees, treeLoggers, maximizer, likelihood, prior, fixed,
                xo.hasAttribute(NOT_IDENTIFIED) ? xo.getStringAttribute(NOT_IDENTIFIED) : null);

        boolean overwrite = xo.getAttribute(LoggerParser.ALLOW_OVERWRITE_LOG, false);
        if (System.getProperty("log.allow.overwrite") != null) {
            overwrite = Boolean.parseBoolean(System.getProperty("log.allow.overwrite", "false"));
        }
        for (File file : exporter.getFiles()) {
            if (file.exists() && !overwrite) {
                throw new XMLParseException("\nThe file " + file.getName() + " already exists in the working " +
                        "directory.\nTo allow it to be overwritten, use the '-overwrite' command line option when " +
                        "running BEAST or select the option in the Run Options dialog box as appropriate.");
            }
        }

        try {
            exporter.export();
        } catch (FileNotFoundException e) {
            throw new XMLParseException("timeTreeExport cannot write a file: " + e.getMessage());
        }
        Logger.getLogger("dr.evomodel").info("timeTreeExport: logged " + treeLoggers.size() + " tree(s); wrote " +
                exporter.getFiles().size() + " files to " + stem.getAbsolutePath() + ".*");
        return exporter;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Logs trees once (logTree elements) and writes tables of dates, rates and parameters, and a run record, " +
                "after a maximization.";
    }

    public Class getReturnType() {
        return TimeTreeExporter.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILE_NAME_STEM),
            AttributeRule.newBooleanRule(LoggerParser.ALLOW_OVERWRITE_LOG, true),
            AttributeRule.newStringRule(NOT_IDENTIFIED, true, "Why the maximum-likelihood fit is not identified: " +
                    "the information criteria are not reported"),
            new ElementRule(TreeLogger.class, 0, Integer.MAX_VALUE),
            new ElementRule(MaximizerWrtParameter.class, true),
            new ElementRule(TREE_LIKELIHOODS, new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}, true),
            new ElementRule(PRIORS, new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}, true),
            new ElementRule(TIME_TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class),
                    new ElementRule(BranchRateModel.class, true),
                    new ElementRule(TaxonList.class, 0, Integer.MAX_VALUE),
            }, 1, Integer.MAX_VALUE),
            new ElementRule(NOT_OPTIMIZED, new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(REASON),
                    AttributeRule.newBooleanRule(CAN_OPTIMIZE, true, "false when the parameters are held only " +
                            "because they cannot be optimized yet: the estimates are then not a full MLE"),
                    new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
            }, 0, Integer.MAX_VALUE),
    };
}
