/*
 * TimeTreeSearchParser.java
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

package dr.evomodelxml.tree.search;

import dr.evolution.util.Taxa;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.search.BranchRateWarmStart;
import dr.evomodel.tree.search.SearchSettings;
import dr.evomodel.tree.search.TimeTreeSearch;
import dr.evomodel.tree.search.TopologyConstraints;
import dr.evomodel.tree.search.TrialEvaluator;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.MaximizerWrtParameter;
import dr.evomodel.treedatalikelihood.discrete.NodeHeightProxyParameter;
import dr.evomodelxml.continuous.hmc.MaximizeWrtParameterParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * {@code <timeTreeSearch>}: a deterministic topology search for a TreeLime time tree (plan 7 §4.3). It needs the
 * tree, the objective, a round maximizer over node heights (through a node-height proxy) and branch rates (or a
 * strict clock's rate), and the final maximizer over everything; both maximizers must maximize the objective, must
 * not have run on parse, and must not include the transform's Jacobian (which would favour some tree shapes).
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class TimeTreeSearchParser extends AbstractXMLObjectParser {

    public static final String NAME = "timeTreeSearch";
    private static final String OBJECTIVE = "objective";
    private static final String ROUND_MAXIMIZER = "roundMaximizer";
    private static final String FINAL_MAXIMIZER = "finalMaximizer";
    private static final String BRANCH_RATES = "branchRates";
    private static final String CONSTRAINTS = "constraints";
    private static final String SEARCH_ON_PARSE = "searchOnParse";
    private static final String MODE = "mode";
    private static final String SEED = "seed";
    private static final String EPSILON = "epsilon";
    private static final String MOVE_EPSILON = "moveEpsilon";
    private static final String DELTA_T1 = "deltaT1";
    private static final String DELTA_KEEP = "deltaKeep";
    private static final String T1_SWEEPS = "t1Sweeps";
    private static final String T2_CANDIDATES = "t2Candidates";
    private static final String T2_ITERATIONS = "t2Iterations";
    private static final String SPEED_NNI = "speedNNI";
    private static final String BATCH = "batch";
    private static final String TOPLIST_SIZE = "toplistSize";
    private static final String SECOND_CHANCE_WINDOW = "secondChanceWindow";
    private static final String SECOND_CHANCE_ITERATIONS = "secondChanceIterations";
    private static final String ROOT_SCAN = "rootScan";
    private static final String PRESET = "preset";
    private static final String MOVES = "moves";
    private static final String SPR_MAX_RADIUS = "sprMaxRadius";
    private static final String SPR_RADIUS = "sprRadius";
    private static final String SUBTREE_CUTOFF = "subtreeCutoff";
    private static final String SWAP_MAX_DISTANCE = "swapMaxDistance";
    private static final String ROOT_MAX_DEPTH = "rootMaxDepth";
    private static final String STARTING_TREE = "startingTree";
    private static final String N_STOP = "nStop";
    private static final String PERTURB_STRENGTH = "perturbStrength";
    private static final String PERTURB_OPERATORS = "perturbOperators";
    private static final String LEAP_SIZE = "leapSize";
    private static final String CANDIDATE_SET_SIZE = "candidateSetSize";
    private static final String PARENT_POOL = "parentPool";
    private static final String NNI_DELTA = "nniDelta";
    private static final String CANDIDATES_FILE_NAME = "candidatesFileName";
    private static final String MAX_ROUNDS = "maxRounds";
    private static final String MAX_RUN_TIME = "maxRunTime";
    private static final String MAX_EVALUATIONS = "maxEvaluations";
    private static final String RATE_WARM_START = "rateWarmStart";
    private static final String VERIFY_EVERY = "verifyEvery";
    private static final String PRINT_SCREEN = "printToScreen";
    private static final String TREE_FILE_NAME = "treeFileName";
    private static final String LOG_FILE_NAME = "logFileName";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
        final Likelihood objective = (Likelihood) xo.getChild(OBJECTIVE).getChild(Likelihood.class);
        final MaximizerWrtParameter round =
                (MaximizerWrtParameter) xo.getChild(ROUND_MAXIMIZER).getChild(MaximizerWrtParameter.class);
        final MaximizerWrtParameter fin =
                (MaximizerWrtParameter) xo.getChild(FINAL_MAXIMIZER).getChild(MaximizerWrtParameter.class);

        if (!(tree instanceof DefaultTreeModel)) {
            throw new XMLParseException("The topology search needs a plain tree model, not a fixed or constrained " +
                    "one (" + tree.getClass().getSimpleName() + ")");
        }
        checkMaximizer(ROUND_MAXIMIZER, round, objective);
        checkMaximizer(FINAL_MAXIMIZER, fin, objective);

        final Set<Model> models = TrialEvaluator.models(objective.getModel());
        for (Model model : models) {
            final String name = model.getClass().getSimpleName();
            if (name.equals("LocalClockModel") || name.equals("RandomLocalClockModel")) {
                throw new XMLParseException("The topology search does not support local clocks (" + name + "): " +
                        "their rates are attached to clades");
            }
        }

        ArbitraryBranchRates rates = null;
        if (xo.hasChildNamed(BRANCH_RATES)) {
            rates = (ArbitraryBranchRates) xo.getChild(BRANCH_RATES).getChild(ArbitraryBranchRates.class);
        } else {
            for (Model model : models) {
                if (model instanceof ArbitraryBranchRates && ((ArbitraryBranchRates) model).getTree() == tree) {
                    rates = (ArbitraryBranchRates) model;
                }
            }
        }
        if (rates != null && likelihoodOnly(objective)) {
            throw new XMLParseException("Maximum likelihood with free per-branch rates does not identify the root or " +
                    "the node times: every rooting and dating reproduces the same branch lengths. Use the joint " +
                    "density (MAP) or a strict clock");
        }

        // the round maximizer: node heights through a proxy of this tree, and parameters TreeModel's store window
        // brings back
        boolean heights = false;
        List<Parameter> leaves = new ArrayList<Parameter>();
        leaves(round.getParameter(), leaves);
        for (Parameter parameter : leaves) {
            if (parameter instanceof NodeHeightProxyParameter) {
                if (((NodeHeightProxyParameter) parameter).getTree() != tree) {
                    throw new XMLParseException("The round maximizer's node heights belong to another tree");
                }
                heights = true;
            } else if (!isStored(parameter, models)) {
                throw new XMLParseException("The round maximizer's parameter " + parameter.getId() + " is not part " +
                        "of a model of the objective, so a rejected trial would not restore it; give the round " +
                        "maximizer node heights and branch rates (or a clock rate) only");
            }
        }
        if (!heights) {
            throw new XMLParseException("The round maximizer does not optimize the node heights (a <nodeHeightProxyParameter> " +
                    "in its gradient): the search needs them, and they have analytic gradients only");
        }

        List<Set<String>> extra = new ArrayList<Set<String>>();
        if (xo.hasChildNamed(CONSTRAINTS)) {
            for (Object child : xo.getChild(CONSTRAINTS).getAllChildren(Taxa.class)) {
                Taxa taxa = (Taxa) child;
                Set<String> ids = new HashSet<String>();
                for (int i = 0; i < taxa.getTaxonCount(); ++i) {
                    ids.add(taxa.getTaxonId(i));
                }
                extra.add(ids);
            }
        }
        TopologyConstraints constraints;
        try {
            constraints = new TopologyConstraints(tree, objective, extra);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }

        // a preset's values are the defaults of the attributes below
        SearchSettings settings = xo.hasAttribute(PRESET) ?
                SearchSettings.preset(xo.getStringAttribute(PRESET)) : new SearchSettings();
        if (settings == null) {
            throw new XMLParseException(PRESET + " is fast, standard or thorough, not " +
                    xo.getStringAttribute(PRESET));
        }
        settings.epsilon = xo.getAttribute(EPSILON, settings.epsilon);
        settings.moveEpsilon = xo.getAttribute(MOVE_EPSILON, settings.moveEpsilon);
        settings.deltaT1 = xo.getAttribute(DELTA_T1, settings.deltaT1);
        settings.deltaKeep = xo.getAttribute(DELTA_KEEP, settings.deltaKeep);
        settings.t1Sweeps = xo.getAttribute(T1_SWEEPS, settings.t1Sweeps);
        settings.t2Candidates = xo.getAttribute(T2_CANDIDATES, settings.t2Candidates);
        settings.t2Iterations = xo.getAttribute(T2_ITERATIONS, settings.t2Iterations);
        settings.speedNni = xo.getAttribute(SPEED_NNI, settings.speedNni);
        settings.batch = xo.getAttribute(BATCH, settings.batch);
        settings.toplistSize = xo.getAttribute(TOPLIST_SIZE, settings.toplistSize);
        settings.secondChanceWindow = xo.getAttribute(SECOND_CHANCE_WINDOW, settings.secondChanceWindow);
        settings.secondChanceIterations = xo.getAttribute(SECOND_CHANCE_ITERATIONS, settings.secondChanceIterations);
        settings.rootScan = xo.getAttribute(ROOT_SCAN, settings.rootScan);
        final String moves = xo.getAttribute(MOVES, settings.moves()).replace(" ", "");
        if (!moves.matches("nni(,spr)?(,swap)?(,root)?")) {
            throw new XMLParseException(MOVES + " is nni, optionally followed by ,spr, ,swap and ,root, not " + moves);
        }
        settings.spr = moves.contains("spr");
        settings.swap = moves.contains("swap");
        settings.rootMoves = moves.contains("root");
        settings.rootMaxDepth = xo.getAttribute(ROOT_MAX_DEPTH, settings.rootMaxDepth);
        final String start = xo.getAttribute(STARTING_TREE, settings.startingTree());
        if (!start.equals("current") && !start.equals("nj") && !start.equals("both")) {
            throw new XMLParseException(STARTING_TREE + " is current, nj or both, not " + start);
        }
        settings.currentStart = !start.equals("nj");
        settings.njStart = !start.equals("current");
        settings.swapMaxDistance = xo.getAttribute(SWAP_MAX_DISTANCE, settings.swapMaxDistance);
        settings.nStop = xo.getAttribute(N_STOP, settings.nStop);
        settings.perturbStrength = xo.getAttribute(PERTURB_STRENGTH, settings.perturbStrength);
        settings.perturbOperators = xo.getAttribute(PERTURB_OPERATORS, settings.perturbOperators);
        settings.leapSize = xo.getAttribute(LEAP_SIZE, settings.leapSize);
        settings.candidateSetSize = xo.getAttribute(CANDIDATE_SET_SIZE, settings.candidateSetSize);
        settings.parentPool = xo.getAttribute(PARENT_POOL, settings.parentPool);
        settings.nniDelta = xo.getAttribute(NNI_DELTA, settings.nniDelta);
        settings.candidatesFileName = xo.hasAttribute(CANDIDATES_FILE_NAME) ?
                xo.getStringAttribute(CANDIDATES_FILE_NAME) : null;
        settings.sprMaxRadius = xo.getAttribute(SPR_MAX_RADIUS, settings.sprMaxRadius);
        settings.sprRadius = xo.getAttribute(SPR_RADIUS, settings.sprRadius);
        settings.subtreeCutoff = xo.getAttribute(SUBTREE_CUTOFF, settings.subtreeCutoff);
        settings.maxRounds = xo.getAttribute(MAX_ROUNDS, 0);
        settings.maxRunTime = xo.hasAttribute(MAX_RUN_TIME) ?
                MaximizeWrtParameterParser.parseRunTime(xo.getStringAttribute(MAX_RUN_TIME)) : 0.0;
        settings.maxEvaluations = xo.getLongIntegerAttribute(MAX_EVALUATIONS, 0);
        settings.verifyEvery = xo.getAttribute(VERIFY_EVERY, 0);
        settings.seed = xo.hasAttribute(SEED) ? xo.getLongIntegerAttribute(SEED) : -1;
        settings.printToScreen = xo.getAttribute(PRINT_SCREEN, true);
        settings.treeFileName = xo.hasAttribute(TREE_FILE_NAME) ? xo.getStringAttribute(TREE_FILE_NAME) : null;
        settings.logFileName = xo.hasAttribute(LOG_FILE_NAME) ? xo.getStringAttribute(LOG_FILE_NAME) : null;
        try {
            settings.rateWarmStart = BranchRateWarmStart.Mode.parse(
                    xo.getAttribute(RATE_WARM_START, BranchRateWarmStart.Mode.KEEP.getName()));
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
        final String mode = xo.getAttribute(MODE, "search");
        if (!mode.equals("search") && !mode.equals("polish")) {
            throw new XMLParseException(MODE + " is search or polish, not " + mode);
        }
        if (settings.epsilon < 0 || settings.moveEpsilon < 0 || settings.deltaT1 < 0 || settings.deltaKeep < 0 ||
                settings.t1Sweeps < 0 || settings.t2Candidates < 0 || settings.t2Iterations < 0 ||
                settings.maxRounds < 0 || settings.maxEvaluations < 0 || settings.verifyEvery < 0 ||
                settings.toplistSize < 0 || settings.sprMaxRadius < 0 || settings.sprRadius < 0 ||
                settings.subtreeCutoff < 0 || settings.swapMaxDistance < 0 || settings.rootMaxDepth < 0 || settings.nStop < 0 ||
                settings.perturbStrength < 0 || settings.leapSize <= 0 || settings.candidateSetSize < 1 ||
                settings.parentPool < 1 || settings.secondChanceWindow < 0 || settings.secondChanceIterations < 0) {
            throw new XMLParseException("The search's thresholds, counts and budgets must not be negative");
        }

        TimeTreeSearch search;
        try {
            search = new TimeTreeSearch(tree, objective, round, fin, rates, constraints, settings);
        } catch (IllegalArgumentException e) {
            throw new XMLParseException(e.getMessage());
        }
        if (xo.getAttribute(SEARCH_ON_PARSE, true)) {
            search.run();
        }
        return search;
    }

    private static void checkMaximizer(String role, MaximizerWrtParameter maximizer, Likelihood objective)
            throws XMLParseException {
        if (maximizer.getLikelihood() != objective) {
            throw new XMLParseException("The " + role + " maximizes another density than the search's objective");
        }
        if (maximizer.wasExecuted()) {
            throw new XMLParseException("The " + role + " has already run: give it maximizeOnParse=\"false\"");
        }
        if (maximizer.includesJacobian()) {
            throw new XMLParseException("The " + role + " includes the transform's Jacobian, which depends on the " +
                    "topology: use includeJacobian=\"false\"");
        }
    }

    private static boolean isStored(Parameter parameter, Set<Model> models) {
        for (Model model : models) {
            if (model instanceof AbstractModel) {
                AbstractModel abstractModel = (AbstractModel) model;
                for (int i = 0; i < abstractModel.getVariableCount(); ++i) {
                    if (abstractModel.getVariable(i) == parameter) {
                        return true;
                    }
                }
            }
        }
        return false;
    }

    private static void leaves(Parameter parameter, List<Parameter> leaves) {
        if (parameter instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) parameter;
            for (int i = 0; i < compound.getParameterCount(); ++i) {
                leaves(compound.getParameter(i), leaves);
            }
        } else {
            leaves.add(parameter);
        }
    }

    // the objective is tree likelihoods only: maximum likelihood
    private static boolean likelihoodOnly(Likelihood likelihood) {
        if (likelihood instanceof CompoundLikelihood) {
            for (Likelihood term : ((CompoundLikelihood) likelihood).getLikelihoods()) {
                if (!likelihoodOnly(term)) {
                    return false;
                }
            }
            return true;
        }
        return likelihood instanceof TreeDataLikelihood;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "A deterministic topology search for a time tree: rooted NNI and clamped-height SPR moves, with " +
                "height repair and a rate warm start, scored in tiers and accepted only when the objective rises.";
    }

    @Override
    public Class getReturnType() {
        return TimeTreeSearch.class;
    }

    @Override
    public String getParserName() {
        return NAME;
    }

    private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TreeModel.class),
            new ElementRule(OBJECTIVE, new XMLSyntaxRule[]{new ElementRule(Likelihood.class)}),
            new ElementRule(ROUND_MAXIMIZER, new XMLSyntaxRule[]{new ElementRule(MaximizerWrtParameter.class)},
                    "maximizes node heights and branch rates (or a clock rate) after topology moves"),
            new ElementRule(FINAL_MAXIMIZER, new XMLSyntaxRule[]{new ElementRule(MaximizerWrtParameter.class)},
                    "maximizes every parameter, at the start, on a new best and at the end"),
            new ElementRule(BRANCH_RATES, new XMLSyntaxRule[]{new ElementRule(ArbitraryBranchRates.class)},
                    "per-branch rates that moves warm-start (default: found in the objective)", true),
            new ElementRule(CONSTRAINTS, new XMLSyntaxRule[]{new ElementRule(Taxa.class, 0, Integer.MAX_VALUE)},
                    "taxon sets that must stay clades, besides the objective's monophyly statistics", true),
            AttributeRule.newBooleanRule(SEARCH_ON_PARSE, true, "run while parsing (default true)"),
            new StringAttributeRule(MODE, "search (default) or polish (both hill-climb from the current tree)", true),
            AttributeRule.newLongIntegerRule(SEED, true, "seed for the random number generator"),
            AttributeRule.newDoubleRule(EPSILON, true, "a round must gain more than this to continue (default 0.1)"),
            AttributeRule.newDoubleRule(MOVE_EPSILON, true, "a move must gain more than this (default 0.01)"),
            AttributeRule.newDoubleRule(DELTA_T1, true, "local optimization for moves whose first score is within this of the incumbent (default 50)"),
            AttributeRule.newDoubleRule(DELTA_KEEP, true, "a re-fit (T2) for moves within this of the incumbent when none gains (default 50)"),
            AttributeRule.newIntegerRule(T1_SWEEPS, true, "sweeps of the local optimization (default 2; 0: off)"),
            AttributeRule.newIntegerRule(T2_CANDIDATES, true, "moves re-fitted when none gains (default 3)"),
            AttributeRule.newIntegerRule(T2_ITERATIONS, true, "iterations of a T2 re-fit (default 20)"),
            AttributeRule.newBooleanRule(SPEED_NNI, true, "after a change, score the NNIs near it first (default true)"),
            AttributeRule.newBooleanRule(BATCH, true, "apply all non-conflicting gaining NNIs at once (default true)"),
            AttributeRule.newIntegerRule(TOPLIST_SIZE, true, "the closest moves kept for the re-fit when none gains (default 20)"),
            AttributeRule.newDoubleRule(SECOND_CHANCE_WINDOW, true, "a re-fitted move this close below the incumbent gets a longer re-fit (default 2)"),
            AttributeRule.newIntegerRule(SECOND_CHANCE_ITERATIONS, true, "iterations of that longer re-fit (default 80; 0: off)"),
            new StringAttributeRule(PRESET, "fast: moves nni from the current tree; standard: moves nni,spr,root from the nj tree; thorough: moves nni,spr,swap,root from the nj tree with nStop 20. The other attributes override the preset's values", true),
            new StringAttributeRule(MOVES, "nni (default, or the preset's), optionally followed by ,spr, ,swap and ,root (in that order): clamped-height SPR rounds, subtree swaps (wide exchanges), then root-changing SPR, when no NNI gains", true),
            AttributeRule.newIntegerRule(N_STOP, true, "perturb and hill-climb again until this many climbs in a row bring no new best (default 0: no perturbation)"),
            AttributeRule.newDoubleRule(PERTURB_STRENGTH, true, "operator moves per perturbation, times (taxa - 2) (default 0.3)"),
            new StringAttributeRule(PERTURB_OPERATORS, "the BEAST operators perturbation draws from: nni, narrow, wide, leap, wb (default nni,wide,leap,wb)", true),
            AttributeRule.newDoubleRule(LEAP_SIZE, true, "subtree leap's distance scale, relative to the root height (default 0.1)"),
            AttributeRule.newIntegerRule(CANDIDATE_SET_SIZE, true, "the best distinct trees kept (default 20)"),
            AttributeRule.newIntegerRule(PARENT_POOL, true, "perturb one of this many best trees (default 5)"),
            AttributeRule.newBooleanRule(NNI_DELTA, true, "at the end, F minus the best NNI alternative at each node (default true)"),
            new StringAttributeRule(CANDIDATES_FILE_NAME, "the NEXUS file for the candidate trees", true),
            new StringAttributeRule(STARTING_TREE, "current (default, or the preset's): the tree as it is; nj: the neighbour-joining tree of the alignment, rooted by root-to-tip regression on the dates and dated; both: a climb from each, the NJ one second", true),
            AttributeRule.newIntegerRule(ROOT_MAX_DEPTH, true, "root-changing SPR: the deepest new place of the old root below the new one, in branches (default 3)"),
            AttributeRule.newIntegerRule(SWAP_MAX_DISTANCE, true, "the farthest apart the parents of swapped subtrees may be, in nodes (default 4)"),
            AttributeRule.newIntegerRule(SPR_MAX_RADIUS, true, "the farthest SPR regraft, in branches (default 10)"),
            AttributeRule.newIntegerRule(SPR_RADIUS, true, "the fast SPR rounds' radius (default 0: found by trying 5, 10, ...)"),
            AttributeRule.newDoubleRule(SUBTREE_CUTOFF, true, "SPR targets beyond one scoring this times the mean drop below the incumbent are skipped (default 1)"),
            AttributeRule.newBooleanRule(ROOT_SCAN, true, "after a new best, try the root positions one and two edges away (default true)"),
            AttributeRule.newIntegerRule(MAX_ROUNDS, true, "at most this many rounds (default: the number of internal nodes)"),
            new StringAttributeRule(MAX_RUN_TIME, "stop after this much wall-clock time, e.g. 30min (default: no limit)", true),
            AttributeRule.newLongIntegerRule(MAX_EVALUATIONS, true, "stop after this many objective evaluations (default: no limit)"),
            new StringAttributeRule(RATE_WARM_START, "keep (default: rates stay with their nodes) or preserveSubstitutions", true),
            AttributeRule.newIntegerRule(VERIFY_EVERY, true, "recompute everything after every this many trials and compare (default 0: never)"),
            AttributeRule.newBooleanRule(PRINT_SCREEN, true, "a line per round (default true)"),
            new StringAttributeRule(TREE_FILE_NAME, "the NEXUS file for the final tree", true),
            new StringAttributeRule(LOG_FILE_NAME, "the tab-separated log, one row per phase or round", true),
    };
}
