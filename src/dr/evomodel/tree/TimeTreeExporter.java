/*
 * TimeTreeExporter.java
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

package dr.evomodel.tree;

import dr.app.beast.BeastVersion;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Date;
import dr.evolution.util.TaxonList;
import dr.evolution.util.TimeScale;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.stats.InformationCriteria;
import dr.util.Identifiable;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.TimeZone;

/**
 * Writes the result of a maximization as files a user can read (plan 8 O1): the trees, through BEAST's own tree
 * loggers (&lt;logTree&gt;, logged once as MCMC logs them: the time tree with the clock's rates, the same tree in
 * substitutions), tables of node dates, branch rates and parameters, and a run record. Intervals are not computed here: the se, lower, upper and ci_method columns are NA until a method that owns
 * them (Hessian, bootstrap) fills them.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeExporter {

    public static final String NA = "NA";
    public static final String STATEMENT = "All results are conditional on point estimates and on the topology.";

    /**
     * A tree to export, with the branch rates of its clock (null for none) and the taxon sets whose tMRCAs the dates
     * table lists.
     */
    public static final class ExportedTree {
        final Tree tree;
        final BranchRateModel rates;
        final List<TaxonList> taxonSets;

        public ExportedTree(Tree tree, BranchRateModel rates, List<TaxonList> taxonSets) {
            this.tree = tree;
            this.rates = rates;
            this.taxonSets = taxonSets;
        }
    }

    /**
     * Parameters the maximizer did not move, and why: the user's choice, or a limitation (canOptimize false: they could
     * not be optimized yet, so the estimates are not a full MLE or MAP).
     */
    public static final class FixedParameters {
        final String reason;
        final List<Parameter> parameters;
        final boolean canOptimize;

        public FixedParameters(String reason, List<Parameter> parameters) {
            this(reason, parameters, true);
        }

        public FixedParameters(String reason, List<Parameter> parameters, boolean canOptimize) {
            this.reason = reason;
            this.parameters = parameters;
            this.canOptimize = canOptimize;
        }
    }

    private final File stem;
    private final List<ExportedTree> trees;
    private final MaximizerWrtParameter maximizer;
    private final Likelihood likelihood;
    private final Likelihood prior;
    private final List<FixedParameters> fixed;
    private final List<TreeLogger> treeLoggers;
    private final String notIdentified;

    /**
     * @param stem        the tables and the run record are this path plus a suffix (.dates.tsv, .summary.txt, ...)
     * @param treeLoggers the tree loggers to log once
     * @param maximizer   null when nothing was maximized
     * @param likelihood  the tree likelihoods, or null
     * @param prior       the priors, or null
     * @param notIdentified why a maximum-likelihood fit is not identified (no information criteria), or null
     */
    public TimeTreeExporter(File stem, List<ExportedTree> trees, List<TreeLogger> treeLoggers,
                            MaximizerWrtParameter maximizer, Likelihood likelihood, Likelihood prior,
                            List<FixedParameters> fixed, String notIdentified) {
        this.notIdentified = notIdentified;
        this.stem = stem;
        this.trees = trees;
        this.treeLoggers = treeLoggers;
        this.maximizer = maximizer;
        this.likelihood = likelihood;
        this.prior = prior;
        this.fixed = fixed;
    }

    /**
     * @return the tables and the run record, so they can be checked before anything is written; the tree loggers
     * open their own files
     */
    public List<File> getFiles() {
        List<File> files = new ArrayList<File>();
        for (ExportedTree tree : trees) {
            files.add(file(tree, ".dates.tsv"));
            if (tree.rates != null) {
                files.add(file(tree, ".rates.tsv"));
            }
        }
        if (maximizer != null) {
            files.add(file(null, ".params.tsv"));
        }
        files.add(file(null, ".summary.txt"));
        return files;
    }

    public void export() throws FileNotFoundException {
        // as MCMC logs them, for a single state
        for (TreeLogger logger : treeLoggers) {
            logger.startLogging();
            logger.log(0L); // MCLogger.log(int) does nothing
            logger.stopLogging();
        }
        for (ExportedTree tree : trees) {
            NodeTable table = new NodeTable(tree);
            writeDates(tree, table);
            if (tree.rates != null) {
                writeRates(tree, table);
            }
        }
        if (maximizer != null) {
            writeParameters();
        }
        writeSummary();
    }

    // <stem>.<suffix> for a single tree, <stem>.<tree id>.<suffix> for several
    private File file(ExportedTree tree, String suffix) {
        String name = stem.getName();
        if (tree != null && trees.size() > 1) {
            name += "." + treeName(tree.tree);
        }
        return new File(stem.getParentFile(), name + suffix);
    }

    private static String treeName(Tree tree) {
        return tree instanceof Identifiable && ((Identifiable) tree).getId() != null ?
                ((Identifiable) tree).getId() : "tree";
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Node heights, dates, names and clades

    /**
     * Dates come from this tree's own tips: date(v) = D - h(v), where D is the time of height 0. The global
     * Taxon.getMostRecentDate() also counts the taxa of other trees.
     */
    final class NodeTable {
        final Tree tree;
        final String[] names;
        final String[] clades;
        final double zero; // decimal year of height 0, NaN without tip dates
        final boolean backwards; // ages before the present rather than dates
        final double toYears; // years per unit of height: the tip dates' units, from which the tip heights come

        NodeTable(ExportedTree exported) {
            this.tree = exported.tree;
            final int count = tree.getNodeCount();
            names = new String[count];
            clades = new String[count];

            double zero = Double.NaN;
            double toYears = Double.NaN;
            boolean backwards = false;
            double youngest = Double.POSITIVE_INFINITY;
            for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
                NodeRef tip = tree.getExternalNode(i);
                Date date = tree.getNodeTaxon(tip).getDate();
                if (date == null || Double.isNaN(yearsPer(date.getUnits()))) {
                    zero = Double.NaN;
                    break;
                }
                double height = tree.getNodeHeight(tip);
                if (height < youngest) {
                    youngest = height;
                    toYears = yearsPer(date.getUnits());
                    // getAbsoluteTimeValue is minus the age for dates given as ages
                    zero = (date.getAbsoluteTimeValue() + height) * toYears;
                    backwards = date.isBackwards();
                }
            }
            this.zero = zero;
            this.toYears = toYears;
            this.backwards = backwards;

            List<List<String>> tips = new ArrayList<List<String>>();
            for (int i = 0; i < count; ++i) {
                tips.add(null);
            }
            collect(tree.getRoot(), tips);
        }

        // names and clades (sorted tip names) of the subtree below node
        private void collect(NodeRef node, List<List<String>> tips) {
            final int number = node.getNumber();
            List<String> below = new ArrayList<String>();
            if (tree.isExternal(node)) {
                names[number] = tree.getNodeTaxon(node).getId();
                below.add(names[number]);
            } else {
                names[number] = String.format(Locale.ROOT, "NODE_%07d", number - tree.getExternalNodeCount() + 1);
                for (int i = 0; i < tree.getChildCount(node); ++i) {
                    NodeRef child = tree.getChild(node, i);
                    collect(child, tips);
                    below.addAll(tips.get(child.getNumber()));
                }
                Collections.sort(below);
            }
            tips.set(number, below);
            StringBuilder sb = new StringBuilder();
            for (String tip : below) {
                sb.append(sb.length() > 0 ? "," : "").append(tip);
            }
            clades[number] = sb.toString();
        }

        boolean isDated() {
            return !Double.isNaN(zero);
        }

        double height(NodeRef node) {
            return tree.getNodeHeight(node);
        }

        /**
         * @return the decimal year of the node, or for ages its age in years
         */
        double date(NodeRef node) {
            final double date = zero - height(node) * toYears;
            return backwards ? -date : date;
        }
    }

    private static double yearsPer(Units.Type units) {
        switch (units) {
            case YEARS:
                return 1.0;
            case MONTHS:
                return 1.0 / 12.0;
            case WEEKS:
                return 7.0 / TimeScale.DAYS_PER_YEAR;
            case DAYS:
                return 1.0 / TimeScale.DAYS_PER_YEAR;
            default:
                return Double.NaN;
        }
    }

    /**
     * The calendar date of a decimal year, in BEAST's convention (dr.evolution.util.Date: 2013.0 is 1 January 2013,
     * the fraction is of that year's length in GMT).
     */
    public static String isoDate(double decimalYear) {
        if (Double.isNaN(decimalYear) || decimalYear < 1.0 || decimalYear >= 10000.0) {
            return NA;
        }
        final int year = (int) Math.floor(decimalYear);
        Calendar calendar = Calendar.getInstance(TimeZone.getTimeZone("GMT"));
        calendar.clear();
        calendar.set(year, Calendar.JANUARY, 1, 0, 0, 0);
        final long start = calendar.getTimeInMillis();
        calendar.set(year + 1, Calendar.JANUARY, 1, 0, 0, 0);
        final long end = calendar.getTimeInMillis();
        calendar.setTimeInMillis(start + Math.round((decimalYear - year) * (end - start)));
        SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd", Locale.ROOT);
        format.setTimeZone(TimeZone.getTimeZone("GMT"));
        return format.format(calendar.getTime());
    }

    private static double substitutionLength(Tree tree, BranchRateModel rates, NodeRef node) {
        return rates.getBranchRate(tree, node) * (tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node));
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Tables

    /**
     * One row per node, then one per taxon set at its MRCA. With dates: date (ISO) and numeric_date (decimal year);
     * with ages: age in years and date NA; without either: NA.
     */
    private void writeDates(ExportedTree exported, NodeTable table) throws FileNotFoundException {
        final Tree tree = exported.tree;
        PrintWriter writer = new PrintWriter(file(exported, ".dates.tsv"));
        writer.println("#node\tclade\tdate\t" + (table.backwards ? "age" : "numeric_date") +
                "\theight\tlower\tupper\tci_method");
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            writeDate(writer, table.names[i], table, tree.getNode(i));
        }
        for (TaxonList taxa : exported.taxonSets) {
            NodeRef mrca;
            try {
                mrca = TreeUtils.getCommonAncestorNode(tree, TreeUtils.getLeavesForTaxa(tree, taxa));
            } catch (TreeUtils.MissingTaxonException e) {
                throw new IllegalArgumentException("Taxon set " + taxa.getId() + " has a taxon that is not in tree " +
                        treeName(tree));
            }
            writeDate(writer, "tmrca(" + taxa.getId() + ")", table, mrca);
        }
        writer.close();
    }

    private static void writeDate(PrintWriter writer, String name, NodeTable table, NodeRef node) {
        String date = NA;
        String numeric = NA;
        if (table.isDated()) {
            numeric = Double.toString(table.date(node));
            if (!table.backwards) {
                date = isoDate(table.date(node));
            }
        }
        writer.println(name + "\t" + table.clades[node.getNumber()] + "\t" + date + "\t" + numeric + "\t" +
                table.height(node) + "\t" + NA + "\t" + NA + "\t" + NA);
    }

    private void writeRates(ExportedTree exported, NodeTable table) throws FileNotFoundException {
        final Tree tree = exported.tree;
        PrintWriter writer = new PrintWriter(file(exported, ".rates.tsv"));
        writer.println("branch\tparent\ttime_length\trate\tsubst_length\trate_lower\trate_upper\tci_method");
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) {
                continue;
            }
            NodeRef parent = tree.getParent(node);
            writer.println(table.names[i] + "\t" + table.names[parent.getNumber()] + "\t" +
                    (tree.getNodeHeight(parent) - tree.getNodeHeight(node)) + "\t" +
                    exported.rates.getBranchRate(tree, node) + "\t" + substitutionLength(tree, exported.rates, node) +
                    "\t" + NA + "\t" + NA + "\t" + NA);
        }
        writer.close();
    }

    private void writeParameters() throws FileNotFoundException {
        PrintWriter writer = new PrintWriter(file(null, ".params.tsv"));
        writer.println("name\tvalue\ttransformed\ttransform\tstatus\tgradient\tse\tlower\tupper\tci_method");
        for (MaximizerWrtParameter.Coordinate coordinate : maximizer.getOptimum()) {
            writer.println(coordinate.name + "\t" + coordinate.value + "\t" + coordinate.transformed + "\t" +
                    coordinate.transform + "\t" + coordinate.getStatus() + "\t" + coordinate.gradient + "\t" + NA +
                    "\t" + NA + "\t" + NA + "\t" + NA);
        }
        for (FixedParameters group : fixed) {
            for (Parameter parameter : group.parameters) {
                for (int i = 0; i < parameter.getDimension(); ++i) {
                    writer.println(parameter.getDimensionName(i) + "\t" + parameter.getParameterValue(i) + "\t" + NA +
                            "\t" + NA + "\tfixed:" + group.reason + "\t" + NA + "\t" + NA + "\t" + NA + "\t" + NA +
                            "\t" + NA);
                }
            }
        }
        writer.close();
    }

    // ---------------------------------------------------------------------------------------------------------------
    // Run record

    /**
     * @return MAP, ML, prior, or other when the maximized density is neither
     */
    public String objectiveName() {
        if (maximizer == null) {
            return NA;
        }
        Likelihood objective = maximizer.getLikelihood();
        if (objective == likelihood) {
            return "ML";
        }
        if (objective == prior) {
            return "prior";
        }
        if (objective instanceof CompoundLikelihood) {
            List<Likelihood> terms = ((CompoundLikelihood) objective).getLikelihoods();
            if (likelihood != null && prior != null && terms.size() == 2 && terms.contains(likelihood) &&
                    terms.contains(prior)) {
                return "MAP";
            }
        }
        return "other";
    }

    private static double logDensity(Likelihood density) {
        return density != null ? density.getLogLikelihood() : Double.NaN;
    }

    private static String format(double value) {
        return Double.isNaN(value) ? NA : Double.toString(value);
    }

    /**
     * A partition's log likelihood and the number of sites behind it (the sum of the pattern weights).
     */
    public static final class PartitionTerm {
        public final String name;
        public final double logLikelihood;
        public final double sites;

        PartitionTerm(String name, double logLikelihood, double sites) {
            this.name = name;
            this.logLikelihood = logLikelihood;
            this.sites = sites;
        }
    }

    /**
     * @return the log likelihood of each partition, named by its patterns: each tree likelihood, split by partition
     * when one likelihood computes several (the BEAGLE 3 multi-partition delegate)
     */
    public List<PartitionTerm> getPartitionTerms() {
        List<PartitionTerm> terms = new ArrayList<PartitionTerm>();
        if (likelihood != null) {
            likelihood.getLogLikelihood();
            addPartitionTerms(likelihood, terms);
        }
        // a partition split over several BEAGLE instances is one term per piece, all with the partition's name
        Map<String, PartitionTerm> byName = new LinkedHashMap<String, PartitionTerm>();
        for (PartitionTerm term : terms) {
            PartitionTerm sum = byName.get(term.name);
            byName.put(term.name, sum == null ? term : new PartitionTerm(term.name,
                    sum.logLikelihood + term.logLikelihood, sum.sites + term.sites));
        }
        return new ArrayList<PartitionTerm>(byName.values());
    }

    private static void addPartitionTerms(Likelihood term, List<PartitionTerm> terms) {
        if (term instanceof CompoundLikelihood) {
            for (Likelihood child : ((CompoundLikelihood) term).getLikelihoods()) {
                addPartitionTerms(child, terms);
            }
            return;
        }
        if (term instanceof TreeDataLikelihood) {
            DataLikelihoodDelegate delegate = ((TreeDataLikelihood) term).getDataLikelihoodDelegate();
            if (delegate instanceof MultiPartitionDataLikelihoodDelegate) {
                MultiPartitionDataLikelihoodDelegate multi = (MultiPartitionDataLikelihoodDelegate) delegate;
                double[] byPartition = multi.getLogLikelihoodsByPartition();
                List<PatternList> patterns = multi.getPatternLists();
                for (int i = 0; i < patterns.size(); ++i) {
                    terms.add(new PartitionTerm(name(patterns.get(i), term.getId() + "[" + i + "]"), byPartition[i],
                            sites(patterns.get(i))));
                }
                return;
            }
            if (delegate instanceof BeagleDataLikelihoodDelegate) {
                PatternList patterns = ((BeagleDataLikelihoodDelegate) delegate).getPatternList();
                terms.add(new PartitionTerm(name(patterns, term.getId()), term.getLogLikelihood(), sites(patterns)));
                return;
            }
        }
        terms.add(new PartitionTerm(term.getId() != null ? term.getId() : term.prettyName(), term.getLogLikelihood(),
                Double.NaN));
    }

    private static String name(PatternList patterns, String fallback) {
        return patterns instanceof Identifiable && ((Identifiable) patterns).getId() != null ?
                ((Identifiable) patterns).getId() : fallback;
    }

    private static double sites(PatternList patterns) {
        double sum = 0.0;
        for (double weight : patterns.getPatternWeights()) {
            sum += weight;
        }
        return sum;
    }

    private Map<String, Double> priorTerms() {
        Map<String, Double> terms = new LinkedHashMap<String, Double>();
        if (prior instanceof CompoundLikelihood) {
            int i = 0;
            for (Likelihood term : ((CompoundLikelihood) prior).getLikelihoods()) {
                String name = term.getId() != null ? term.getId() : term.prettyName();
                terms.put(terms.containsKey(name) ? name + "#" + i : name, term.getLogLikelihood());
                ++i;
            }
        }
        return terms;
    }

    private void writeSummary() throws FileNotFoundException {
        PrintWriter writer = new PrintWriter(file(null, ".summary.txt"));
        BeastVersion version = new BeastVersion();
        writer.println("BEAST " + version.getVersionString() + " (" + BeastVersion.getRevision() + ")");
        String commandLine = System.getProperty("command_line");
        writer.println("Command line: " + (commandLine != null ? commandLine.trim() : NA));
        writer.println("Seed: " + MathUtils.getSeed());
        writer.println("Objective: " + objectiveName());

        double logL = logDensity(likelihood);
        double logPrior = logDensity(prior);
        int held = 0;
        if (maximizer != null) {
            writer.println("Maximizer: " + maximizer.describeSettings());
            writer.println("Stop reason: " + maximizer.getReason() + " (" + maximizer.getStatus() + ")");
            writer.println("Iterations: " + maximizer.getIterations() + " (restarts: " + maximizer.getRestarts() + ")");
            writer.println("Wall time: " + maximizer.getSeconds() + " s");
            writer.println("F (objective at the optimum, minimized): " + maximizer.getMinimumValue());
            for (MaximizerWrtParameter.Coordinate coordinate : maximizer.getOptimum()) {
                if (coordinate.held) {
                    ++held;
                }
            }
        }
        writer.println();

        writer.println("Log likelihood: " + format(logL));
        double sites = 0.0;
        for (PartitionTerm term : getPartitionTerms()) {
            writer.println("  " + term.name + ": " + term.logLikelihood + (Double.isNaN(term.sites) ? "" :
                    " (sites: " + formatCount(term.sites) + ")"));
            sites += term.sites;
        }
        writer.println("Log prior: " + format(logPrior));
        for (Map.Entry<String, Double> term : priorTerms().entrySet()) {
            writer.println("  " + term.getKey() + ": " + term.getValue());
        }
        writer.println();

        if (maximizer != null) {
            writer.println("k (free coordinates, node heights included): " + maximizer.getFreeCoordinateCount() +
                    (held > 0 ? " (" + held + " held at a bound)" : ""));
        }
        writer.println("n (weighted sites in the likelihood): " + (Double.isNaN(sites) ? NA : formatCount(sites)));
        writer.println();

        Fit fit = getFit();
        if (fit.reason != null) {
            writer.println("Information criteria: not reported: " + fit.reason + ".");
        } else {
            writer.println("Information criteria (maximum likelihood; lower is better):");
            writer.println("  AIC  = " + fit.aic);
            writer.println("  AICc = " + (Double.isNaN(fit.aicc) ? "undefined (n <= k + 1)" : Double.toString(fit.aicc)));
            writer.println("  BIC  = " + fit.bic);
            if (held > 0) {
                writer.println("  k includes " + held + " coordinate(s) held at a bound: the penalty assumes an " +
                        "interior optimum");
            }
            writer.println("  k counts the internal node heights and the clock, not 2n - 3 branch lengths: compare " +
                    "these criteria between TreeLime fits of the same data.");
        }
        List<String> limited = new ArrayList<String>();
        for (FixedParameters group : fixed) {
            if (!group.canOptimize) {
                for (Parameter parameter : group.parameters) {
                    limited.add(parameter.getId());
                }
            }
        }
        if (!limited.isEmpty()) {
            writer.println("Not a full " + ("ML".equals(objectiveName()) ? "MLE" : "optimum") + ": " +
                    String.join(", ", limited) + " " + (limited.size() > 1 ? "are" : "is") +
                    " held only because TreeLime cannot optimize " + (limited.size() > 1 ? "them" : "it") +
                    " yet; k counts " + (limited.size() > 1 ? "them" : "it") + " as fixed.");
        }
        writer.println();

        List<String> warnings = new ArrayList<String>();
        if (maximizer != null) {
            for (MaximizerWrtParameter.Held hold : maximizer.getHeld()) {
                warnings.add("held at the lower bound: " + hold);
            }
            if (held > 0) {
                warnings.add("coordinates held at a bound are counted in k; information criteria and chi-squared " +
                        "references assume an interior optimum");
            }
            if (!Double.isNaN(maximizer.getWorstGradientCheckError())) {
                warnings.add("worst gradient check error: " + maximizer.getWorstGradientCheckError());
            }
            if (maximizer.includesJacobian()) {
                warnings.add("includeJacobian is true: F includes the transform's log-Jacobian, so -F is not a log " +
                        "density of the model's parameters");
            }
        }
        writer.println("Warnings:" + (warnings.isEmpty() ? " none" : ""));
        for (String warning : warnings) {
            writer.println("  " + warning);
        }
        if (!fixed.isEmpty()) {
            writer.println("Not optimized:");
            for (FixedParameters group : fixed) {
                StringBuilder sb = new StringBuilder();
                for (Parameter parameter : group.parameters) {
                    sb.append(sb.length() > 0 ? ", " : "").append(parameter.getId());
                }
                writer.println("  " + group.reason + (group.canOptimize ? "" : " (cannot be optimized yet)") + ": " + sb);
            }
        }
        writer.println();
        writer.println(STATEMENT);
        writer.close();
    }

    /**
     * The information criteria of a maximum-likelihood fit, or why there are none.
     */
    public static final class Fit {
        public final int k;
        public final double n;
        public final double logLikelihood;
        public final double aic;
        public final double aicc; // NaN when n <= k + 1
        public final double bic;
        public final String reason; // why they are not reported, or null

        Fit(int k, double n, double logLikelihood, String reason) {
            this.k = k;
            this.n = n;
            this.logLikelihood = logLikelihood;
            this.reason = reason;
            boolean reported = reason == null;
            this.aic = reported ? InformationCriteria.aic(logLikelihood, k) : Double.NaN;
            this.aicc = reported ? InformationCriteria.aicc(logLikelihood, k, n) : Double.NaN;
            this.bic = reported ? InformationCriteria.bic(logLikelihood, k, n) : Double.NaN;
        }
    }

    /**
     * k is the number of free coordinates the maximizer moved (node heights included; held ones counted, the
     * simplex's dependent one not), n the number of sites behind the likelihood. Only a maximum-likelihood fit that is
     * identified has information criteria.
     */
    public Fit getFit() {
        if (maximizer == null || likelihood == null) {
            return new Fit(0, Double.NaN, Double.NaN, "nothing was maximized");
        }
        int k = maximizer.getFreeCoordinateCount();
        double n = 0.0;
        for (PartitionTerm term : getPartitionTerms()) {
            n += term.sites;
        }
        double logL = likelihood.getLogLikelihood();
        String reason = null;
        if (!"ML".equals(objectiveName())) {
            reason = "they are defined for a maximum-likelihood fit, and this objective includes the priors; refit " +
                    "with the maximum-likelihood objective";
        } else if (notIdentified != null) {
            reason = "the fit is not identified: " + notIdentified;
        } else if (Double.isNaN(n)) {
            reason = "the number of sites is unknown";
        }
        return new Fit(k, n, logL, reason);
    }

    private static String formatCount(double count) {
        return count == Math.rint(count) ? Long.toString((long) count) : Double.toString(count);
    }
}
