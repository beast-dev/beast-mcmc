/*
 * BeautiOptions.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.app.oldbeauti;

import dr.evolution.alignment.Alignment;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.datatype.AminoAcids;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.TwoStates;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.*;
import dr.evomodel.coalescent.VariableDemographicModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodelxml.speciation.BirthDeathModelParser;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.evoxml.AlignmentParser;
import dr.evoxml.TaxaParser;
import dr.evoxml.TaxonParser;
import dr.util.NumberFormatter;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.jdom.Document;
import org.jdom.Element;

import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class BeautiOptions {
    final public static String LOCAL_CLOCK = "localClock";
    final public static String UCLD_MEAN = "ucld.mean";
    final public static String UCLD_STDEV = "ucld.stdev";
    final public static String UCED_MEAN = "uced.mean";


    public BeautiOptions() {
        double demoWeights = 3.0;
        double substWeights = 1.0;
        double rateWeights = 3.0;
        double branchWeights = 30.0;
        double treeWeights = 15.0;

        createParameter("tree", "The tree");
        createParameter("treeModel.internalNodeHeights", "internal node heights of the tree (except the root)");
        createParameter("treeModel.allInternalNodeHeights", "internal node heights of the tree");
        createParameter("treeModel.rootHeight", "root height of the tree", true, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleParameter("constant.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

        createScaleParameter("exponential.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("exponential.growthRate", "coalescent growth rate parameter", GROWTH_RATE_SCALE, 0.0, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);
        createParameter("exponential.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("logistic.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("logistic.t50", "logistic shape parameter", T50_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createScaleParameter("expansion.popSize", "coalescent population size parameter", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.growthRate", "coalescent logistic growth rate parameter", GROWTH_RATE_SCALE, 0.001, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.doublingTime", "coalescent doubling time parameter", TIME_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("expansion.ancestralProportion", "ancestral population proportion", NONE, 0.1, 0.0, 1.0);
        createParameter("skyline.popSize", "Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("skyline.groupSize", "Bayesian Skyline group sizes");

        createParameter("demographic.popSize", "Extended Bayesian Skyline population sizes", TIME_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("demographic.indicators", "Extended Bayesian Skyline population switch");
        createScaleParameter("demographic.populationMean", "Extended Bayesian Skyline population prior mean", TIME_SCALE, 1, 0, Double.POSITIVE_INFINITY);
        {
            final Parameter p = createStatistic("demographic.populationSizeChanges", "Average number of population change points", true);
            p.priorType = PriorType.POISSON_PRIOR;
            p.poissonMean = Math.log(2);
        }
        createParameter("yule.birthRate", "Yule speciation process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);

//        createParameter("birthDeath.birthRate", "Birth-Death speciation process birth rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
//        createParameter("birthDeath.deathRate", "Birth-Death speciation process death rate", BIRTH_RATE_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, "Birth-Death speciation process rate", BIRTH_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, "Death/Birth speciation process relative death rate", BIRTH_RATE_SCALE, 0.5, 0.0, 1.0);
        //createParameter("birthDeath.samplingProportion", "Birth-Death speciation process sampling proportion", NONE, 1.0, 0.0, 1.0);

        createParameter("clock.rate", "substitution rate", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(UCED_MEAN, "uncorrelated exponential relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(UCLD_MEAN, "uncorrelated lognormal relaxed clock mean", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(UCLD_STDEV, "uncorrelated lognormal relaxed clock stdev", LOG_STDEV_SCALE, 0.1, 0.0, Double.POSITIVE_INFINITY);
        createParameter("branchRates.categories", "relaxed clock branch rate categories");
        createParameter(LOCAL_CLOCK + "." + "rates", "random local clock rates", SUBSTITUTION_RATE_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter(LOCAL_CLOCK + "." + "changes", "random local clock rate change indicator");

        //Substitution model parameters
        createParameter("hky.frequencies", "HKY base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky1.frequencies", "HKY base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky2.frequencies", "HKY base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("hky3.frequencies", "HKY base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("hky.kappa", "HKY transition-transversion parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky1.kappa", "HKY transition-transversion parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky2.kappa", "HKY transition-transversion parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("hky3.kappa", "HKY transition-transversion parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("gtr.frequencies", "GTR base frequencies", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr1.frequencies", "GTR base frequencies for codon position 1", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr2.frequencies", "GTR base frequencies for codon position 2", UNITY_SCALE, 0.25, 0.0, 1.0);
        createParameter("gtr3.frequencies", "GTR base frequencies for codon position 3", UNITY_SCALE, 0.25, 0.0, 1.0);

        createScaleParameter("gtr.ac", "GTR A-C substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.ag", "GTR A-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.at", "GTR A-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.cg", "GTR C-G substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr.gt", "GTR G-T substitution parameter", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr1.ac", "GTR A-C substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.ag", "GTR A-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.at", "GTR A-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.cg", "GTR C-G substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr1.gt", "GTR G-T substitution parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr2.ac", "GTR A-C substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.ag", "GTR A-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.at", "GTR A-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.cg", "GTR C-G substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr2.gt", "GTR G-T substitution parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createScaleParameter("gtr3.ac", "GTR A-C substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.ag", "GTR A-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.at", "GTR A-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.cg", "GTR C-G substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);
        createScaleParameter("gtr3.gt", "GTR G-T substitution parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 1.0E-8, Double.POSITIVE_INFINITY);

        createParameter("bsimple.frequencies", "Binary Simple frequencies", UNITY_SCALE, 0.5, 0.0, 1.0);

        createParameter("bcov.frequencies", "Binary Covarion frequencies of the visible states", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.hfrequencies", "Binary Covarion frequencies of the hidden rates", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.alpha", "Binary Covarion rate of evolution in slow mode", UNITY_SCALE, 0.5, 0.0, 1.0);
        createParameter("bcov.s", "Binary Covarion rate of flipping between slow and fast modes", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, 100.0);

        createParameter(SiteModel.SITE_MODEL + "." + "alpha", "gamma shape parameter", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel1.alpha", "gamma shape parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel2.alpha", "gamma shape parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel3.alpha", "gamma shape parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 0.5, 0.0, Double.POSITIVE_INFINITY);

        createParameter(SiteModel.SITE_MODEL + "." + "pInv", "proportion of invariant sites parameter", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel1.pInv", "proportion of invariant sites parameter for codon position 1", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel2.pInv", "proportion of invariant sites parameter for codon position 2", NONE, 0.5, 0.0, 1.0);
        createParameter("siteModel3.pInv", "proportion of invariant sites parameter for codon position 3", NONE, 0.5, 0.0, 1.0);

        createParameter("siteModel1.mu", "relative rate parameter for codon position 1", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel2.mu", "relative rate parameter for codon position 2", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("siteModel3.mu", "relative rate parameter for codon position 3", SUBSTITUTION_PARAMETER_SCALE, 1.0, 0.0, Double.POSITIVE_INFINITY);
        createParameter("allMus", "All the relative rates");

        // These are statistics which could have priors on...
        createStatistic("meanRate", "The mean rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
        createStatistic(RateStatisticParser.COEFFICIENT_OF_VARIATION, "The variation in rate of evolution over the whole tree", 0.0, Double.POSITIVE_INFINITY);
        createStatistic("covariance", "The covariance in rates of evolution on each lineage with their ancestral lineages", Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY);

        createOperator("constant.popSize", SCALE, 0.75, demoWeights);
        createOperator("exponential.popSize", SCALE, 0.75, demoWeights);
        createOperator("exponential.growthRate", RANDOM_WALK, 1.0, demoWeights);
        createOperator("exponential.doublingTime", SCALE, 0.75, demoWeights);
        createOperator("logistic.popSize", SCALE, 0.75, demoWeights);
        createOperator("logistic.growthRate", SCALE, 0.75, demoWeights);
        createOperator("logistic.doublingTime", SCALE, 0.75, demoWeights);
        createOperator("logistic.t50", SCALE, 0.75, demoWeights);
        createOperator("expansion.popSize", SCALE, 0.75, demoWeights);
        createOperator("expansion.growthRate", SCALE, 0.75, demoWeights);
        createOperator("expansion.doublingTime", SCALE, 0.75, demoWeights);
        createOperator("expansion.ancestralProportion", SCALE, 0.75, demoWeights);
        createOperator("skyline.popSize", SCALE, 0.75, demoWeights * 5);
        createOperator("skyline.groupSize", INTEGER_DELTA_EXCHANGE, 1.0, demoWeights * 2);

        createOperator("demographic.populationMean", SCALE, 0.9, demoWeights);
        createOperator("demographic.indicators", BITFLIP, 1, 2 * treeWeights);
        // hack pass distribution in name
        createOperator("demographic.popSize", "demographic.populationMeanDist", "", "demographic.popSize", "demographic.indicators", SAMPLE_NONACTIVE, 1, 5 * demoWeights);
        createOperator("demographic.scaleActive", "demographic.scaleActive", "", "demographic.popSize", "demographic.indicators", SCALE_WITH_INDICATORS, 0.5, 2 * demoWeights);

        createOperator("yule.birthRate", SCALE, 0.75, demoWeights);
//        createOperator("birthDeath.birthRate", SCALE, 0.75, demoWeights);
//        createOperator("birthDeath.deathRate", SCALE, 0.75, demoWeights);

        createOperator(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME, SCALE, 0.75, demoWeights);
        createOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME, SCALE, 0.75, demoWeights);
        //createOperator("birthDeath.samplingProportion", RANDOM_WALK, 0.75, demoWeights);

        createOperator("clock.rate", SCALE, 0.75, rateWeights);
        createOperator(UCED_MEAN, SCALE, 0.75, rateWeights);
        createOperator(UCLD_MEAN, SCALE, 0.75, rateWeights);
        createOperator(UCLD_STDEV, SCALE, 0.75, rateWeights);
//        createOperator("swapBranchRateCategories", "branchRates.categories", "Performs a swap of branch rate categories", "branchRates.categories", SWAP, 1, branchWeights);
        createOperator("randomWalkBranchRateCategories", "branchRates.categories", "Performs an integer random walk of branch rate categories", "branchRates.categories", INTEGER_RANDOM_WALK, 1, branchWeights);
        createOperator("unformBranchRateCategories", "branchRates.categories", "Performs an integer uniform draw of branch rate categories", "branchRates.categories", INTEGER_UNIFORM, 1, branchWeights);

        createOperator(LOCAL_CLOCK + "." + "rates", SCALE, 0.75, treeWeights);
        createOperator(LOCAL_CLOCK + "." + "changes", BITFLIP, 1, treeWeights);
        createOperator("treeBitMove", "Tree", "Swaps the rates and change locations of local clocks", "tree", TREE_BIT_MOVE, -1.0, treeWeights);

        createOperator("hky.kappa", SCALE, 0.75, substWeights);
        createOperator("hky1.kappa", SCALE, 0.75, substWeights);
        createOperator("hky2.kappa", SCALE, 0.75, substWeights);
        createOperator("hky3.kappa", SCALE, 0.75, substWeights);
        createOperator("hky.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky1.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky2.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("hky3.frequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator("gtr.ac", SCALE, 0.75, substWeights);
        createOperator("gtr.ag", SCALE, 0.75, substWeights);
        createOperator("gtr.at", SCALE, 0.75, substWeights);
        createOperator("gtr.cg", SCALE, 0.75, substWeights);
        createOperator("gtr.gt", SCALE, 0.75, substWeights);

        createOperator("gtr1.ac", SCALE, 0.75, substWeights);
        createOperator("gtr1.ag", SCALE, 0.75, substWeights);
        createOperator("gtr1.at", SCALE, 0.75, substWeights);
        createOperator("gtr1.cg", SCALE, 0.75, substWeights);
        createOperator("gtr1.gt", SCALE, 0.75, substWeights);

        createOperator("gtr2.ac", SCALE, 0.75, substWeights);
        createOperator("gtr2.ag", SCALE, 0.75, substWeights);
        createOperator("gtr2.at", SCALE, 0.75, substWeights);
        createOperator("gtr2.cg", SCALE, 0.75, substWeights);
        createOperator("gtr2.gt", SCALE, 0.75, substWeights);

        createOperator("gtr3.ac", SCALE, 0.75, substWeights);
        createOperator("gtr3.ag", SCALE, 0.75, substWeights);
        createOperator("gtr3.at", SCALE, 0.75, substWeights);
        createOperator("gtr3.cg", SCALE, 0.75, substWeights);
        createOperator("gtr3.gt", SCALE, 0.75, substWeights);

        createOperator("gtr.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr1.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr2.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("gtr3.frequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator("bcov.alpha", SCALE, 0.75, substWeights);
        createOperator("bcov.s", SCALE, 0.75, substWeights);
        createOperator("bcov.frequencies", DELTA_EXCHANGE, 0.01, substWeights);
        createOperator("bcov.hfrequencies", DELTA_EXCHANGE, 0.01, substWeights);

        createOperator(SiteModel.SITE_MODEL + "." + "alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel1.alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel2.alpha", SCALE, 0.75, substWeights);
        createOperator("siteModel3.alpha", SCALE, 0.75, substWeights);

        createOperator(SiteModel.SITE_MODEL + "." + "pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel1.pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel2.pInv", SCALE, 0.75, substWeights);
        createOperator("siteModel3.pInv", SCALE, 0.75, substWeights);

        createOperator("upDownRateHeights", "Substitution rate and heights", "Scales substitution rates inversely to node heights of the tree", "clock.rate", "treeModel.allInternalNodeHeights", UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCEDMeanHeights", "UCED mean and heights", "Scales UCED mean inversely to node heights of the tree", UCED_MEAN, "treeModel.allInternalNodeHeights", UP_DOWN, 0.75, rateWeights);
        createOperator("upDownUCLDMeanHeights", "UCLD mean and heights", "Scales UCLD mean inversely to node heights of the tree", UCLD_MEAN, "treeModel.allInternalNodeHeights", UP_DOWN, 0.75, rateWeights);
        createOperator("centeredMu", "Relative rates", "Scales codon position rates relative to each other maintaining mean", "allMus", CENTERED_SCALE, 0.75, substWeights);
        createOperator("deltaMu", "Relative rates", "Changes codon position rates relative to each other maintaining mean", "allMus", DELTA_EXCHANGE, 0.75, substWeights);

        createOperator("treeModel.rootHeight", SCALE, 0.75, demoWeights);
        createOperator("uniformHeights", "Internal node heights", "Draws new internal node heights uniformally", "treeModel.internalNodeHeights", UNIFORM, -1, branchWeights);

        createOperator("subtreeSlide", "Tree", "Performs the subtree-slide rearrangement of the tree", "tree", SUBTREE_SLIDE, 1.0, treeWeights);
        createOperator("narrowExchange", "Tree", "Performs local rearrangements of the tree", "tree", NARROW_EXCHANGE, -1, treeWeights);
        createOperator("wideExchange", "Tree", "Performs global rearrangements of the tree", "tree", WIDE_EXCHANGE, -1, demoWeights);
        createOperator("wilsonBalding", "Tree", "Performs the Wilson-Balding rearrangement of the tree", "tree", WILSON_BALDING, -1, demoWeights);
    }

    protected void createScaleParameter(String name, String description, int scale, double value, double lower, double upper) {
        Parameter p = createParameter(name, description, scale, value, lower, upper);
        p.priorType = PriorType.JEFFREYS_PRIOR;
    }

    protected Parameter createParameter(String name, String description, int scale, double value, double lower, double upper) {
        final Parameter parameter = new Parameter(name, description, scale, value, lower, upper);
        parameters.put(name, parameter);
        return parameter;
    }

    protected Parameter createParameter(String name, String description) {
        final Parameter parameter = new Parameter(name, description);
        parameters.put(name, parameter);
        return parameter;
    }

    protected Parameter createStatistic(String name, String description, boolean isDiscrete) {
        final Parameter parameter = new Parameter(name, description, isDiscrete);
        parameters.put(name, parameter);
        return parameter;
    }

    protected void createStatistic(String name, String description, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, lower, upper));
    }

    protected void createParameter(String name, String description, boolean isNodeHeight, double value, double lower, double upper) {
        parameters.put(name, new Parameter(name, description, isNodeHeight, value, lower, upper));
    }

    protected void createOperator(String parameterName, String type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(parameter.name, new Operator(parameterName, "", parameter, type, tuning, weight));
    }

    protected void createOperator(String key, String name, String description, String parameterName, String type, double tuning, double weight) {
        Parameter parameter = getParameter(parameterName);
        operators.put(key, new Operator(name, description, parameter, type, tuning, weight));
    }

    protected void createOperator(String key, String name, String description, String parameterName1, String parameterName2, String type, double tuning, double weight) {
        Parameter parameter1 = getParameter(parameterName1);
        Parameter parameter2 = getParameter(parameterName2);
        operators.put(key, new Operator(name, description, parameter1, parameter2, type, tuning, weight));
    }

    private double round(double value, int sf) {
        NumberFormatter formatter = new NumberFormatter(sf);
        try {
            return NumberFormat.getInstance().parse(formatter.format(value)).doubleValue();
        } catch (ParseException e) {
            return value;
        }
    }

    /**
     * return an list of operators that are required
     *
     * @return the parameter list
     */
    public ArrayList<Parameter> selectParameters() {

        ArrayList<Parameter> ops = new ArrayList<Parameter>();

        selectParameters(ops);
        selectStatistics(ops);

        double growthRateMaximum = 1E6;
        double birthRateMaximum = 1E6;
        double substitutionRateMaximum = 100;
        double logStdevMaximum = 10;
        double substitutionParameterMaximum = 100;
        double initialRootHeight = 1;
        double initialRate = 1;


        if (fixedSubstitutionRate) {
            double rate = meanSubstitutionRate;

            growthRateMaximum = 1E6 * rate;
            birthRateMaximum = 1E6 * rate;

            if (alignment != null) {
                initialRootHeight = meanDistance / rate;

                initialRootHeight = round(initialRootHeight, 2);
            }

        } else {
            if (maximumTipHeight > 0) {
                initialRootHeight = maximumTipHeight * 10.0;
            }

            initialRate = round((meanDistance * 0.2) / initialRootHeight, 2);
        }

        double timeScaleMaximum = round(initialRootHeight * 1000.0, 2);

        for (Parameter param : ops) {
            if (alignmentReset) param.priorEdited = false;

            if (!param.priorEdited) {
                switch (param.scale) {
                    case TIME_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = initialRootHeight;
                        break;
                    case T50_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(timeScaleMaximum, param.upper);
                        param.initial = initialRootHeight / 5.0;
                        break;
                    case GROWTH_RATE_SCALE:
                        param.uniformLower = Math.max(-growthRateMaximum, param.lower);
                        param.uniformUpper = Math.min(growthRateMaximum, param.upper);
                        break;
                    case BIRTH_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(birthRateMaximum, param.upper);
                        break;
                    case SUBSTITUTION_RATE_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionRateMaximum, param.upper);
                        param.initial = initialRate;
                        break;
                    case LOG_STDEV_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(logStdevMaximum, param.upper);
                        break;
                    case SUBSTITUTION_PARAMETER_SCALE:
                        param.uniformLower = Math.max(0.0, param.lower);
                        param.uniformUpper = Math.min(substitutionParameterMaximum, param.upper);
                        break;

                    case UNITY_SCALE:
                        param.uniformLower = 0.0;
                        param.uniformUpper = 1.0;
                        break;

                }
                if (param.isNodeHeight) {
                    param.lower = maximumTipHeight;
                    param.uniformLower = maximumTipHeight;
                    param.uniformUpper = timeScaleMaximum;
                    param.initial = initialRootHeight;
                }
            }
        }

        alignmentReset = false;

        return ops;
    }

    /**
     * return an list of operators that are required
     *
     * @return the operator list
     */
    public ArrayList<Operator> selectOperators() {

        ArrayList<Operator> ops = new ArrayList<Operator>();

        selectOperators(ops);

        double initialRootHeight = 1;

        if (fixedSubstitutionRate) {
            double rate = meanSubstitutionRate;

            if (alignment != null) {
                initialRootHeight = meanDistance / rate;
                initialRootHeight = round(initialRootHeight, 2);
            }

        } else {
            if (maximumTipHeight > 0) {
                initialRootHeight = maximumTipHeight * 10.0;
            }
        }

        Operator op = getOperator("subtreeSlide");
        if (!op.tuningEdited) {
            op.tuning = initialRootHeight / 10.0;
        }

        return ops;
    }

    /**
     * return a list of parameters that are required
     *
     * @param params the parameter list
     */
    private void selectParameters(ArrayList<Parameter> params) {

        if (alignment != null) {

            if (partitionCount > 1) {
                for (int i = 1; i <= partitionCount; i++) {
                    params.add(getParameter(SiteModel.SITE_MODEL + i + ".mu"));
                }
            }
            switch (dataType) {
                case DataType.NUCLEOTIDES:
                    switch (nucSubstitutionModel) {
                        case HKY:
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    params.add(getParameter("hky" + i + ".kappa"));
                                }
                            } else {
                                params.add(getParameter("hky.kappa"));
                            }
                            break;
                        case GTR:
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    params.add(getParameter("gtr" + i + ".ac"));
                                    params.add(getParameter("gtr" + i + ".ag"));
                                    params.add(getParameter("gtr" + i + ".at"));
                                    params.add(getParameter("gtr" + i + ".cg"));
                                    params.add(getParameter("gtr" + i + ".gt"));
                                }
                            } else {
                                params.add(getParameter("gtr.ac"));
                                params.add(getParameter("gtr.ag"));
                                params.add(getParameter("gtr.at"));
                                params.add(getParameter("gtr.cg"));
                                params.add(getParameter("gtr.gt"));
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }
                    break;

                case DataType.AMINO_ACIDS:
                    break;

                case DataType.TWO_STATES:
                case DataType.COVARION:
                    switch (binarySubstitutionModel) {
                        case BIN_SIMPLE:
                            break;

                        case BIN_COVARION:
                            params.add(getParameter("bcov.alpha"));
                            params.add(getParameter("bcov.s"));
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown binary substitution model");
                    }
                    break;

                default:
                    throw new IllegalArgumentException("Unknown data type");
            }

            // if gamma do shape move
            if (gammaHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        params.add(getParameter(SiteModel.SITE_MODEL + i + ".alpha"));
                    }
                } else {
                    params.add(getParameter(SiteModel.SITE_MODEL + "." + "alpha"));
                }
            }
            // if pinv do pinv move
            if (invarHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        params.add(getParameter(SiteModel.SITE_MODEL + i + ".pInv"));
                    }
                } else {
                    params.add(getParameter(SiteModel.SITE_MODEL + "." + "pInv"));
                }
            }

            // if not fixed then do mutation rate move and up/down move
            if (!fixedSubstitutionRate) {
                Parameter rateParam;

                if (clockModel == STRICT_CLOCK || clockModel == RANDOM_LOCAL_CLOCK) {
                    rateParam = getParameter("clock.rate");
                    params.add(rateParam);
                } else {
                    if (clockModel == UNCORRELATED_EXPONENTIAL) {
                        rateParam = getParameter(UCED_MEAN);
                        params.add(rateParam);
                    } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                        rateParam = getParameter(UCLD_MEAN);
                        params.add(rateParam);
                        params.add(getParameter(UCLD_STDEV));
                    } else {
                        throw new IllegalArgumentException("Unknown clock model");
                    }
                }

                rateParam.isFixed = false;
            } else {
                Parameter rateParam;
                if (clockModel == STRICT_CLOCK || clockModel == RANDOM_LOCAL_CLOCK) {
                    rateParam = getParameter("clock.rate");
                } else {
                    if (clockModel == UNCORRELATED_EXPONENTIAL) {
                        rateParam = getParameter(UCED_MEAN);
                    } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                        rateParam = getParameter(UCLD_MEAN);
                        params.add(getParameter(UCLD_STDEV));
                    } else {
                        throw new IllegalArgumentException("Unknown clock model");
                    }
                }
                rateParam.isFixed = true;
            }
        }

        if (nodeHeightPrior == CONSTANT) {
            params.add(getParameter("constant.popSize"));
        } else if (nodeHeightPrior == EXPONENTIAL) {
            params.add(getParameter("exponential.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("exponential.growthRate"));
            } else {
                params.add(getParameter("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == LOGISTIC) {
            params.add(getParameter("logistic.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("logistic.growthRate"));
            } else {
                params.add(getParameter("logistic.doublingTime"));
            }
            params.add(getParameter("logistic.t50"));
        } else if (nodeHeightPrior == EXPANSION) {
            params.add(getParameter("expansion.popSize"));
            if (parameterization == GROWTH_RATE) {
                params.add(getParameter("expansion.growthRate"));
            } else {
                params.add(getParameter("expansion.doublingTime"));
            }
            params.add(getParameter("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == SKYLINE) {
            params.add(getParameter("skyline.popSize"));
        } else if (nodeHeightPrior == EXTENDED_SKYLINE) {
            params.add(getParameter("demographic.populationSizeChanges"));
            params.add(getParameter("demographic.populationMean"));
        } else if (nodeHeightPrior == YULE) {
            params.add(getParameter("yule.birthRate"));
        } else if (nodeHeightPrior == BIRTH_DEATH) {
//            params.add(getParameter("birthDeath.birthRate"));
//            params.add(getParameter("birthDeath.deathRate"));
            params.add(getParameter(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            params.add(getParameter(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            // at present we are not allowing the sampling of samplingProportion
        }

        params.add(getParameter("treeModel.rootHeight"));
    }

    private void selectStatistics(ArrayList<Parameter> params) {

        if (taxonSets != null) {
            for (Taxa taxonSet : taxonSets) {
                Parameter statistic = statistics.get(taxonSet);
                if (statistic == null) {
                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
                    statistics.put(taxonSet, statistic);
                }
                params.add(statistic);
            }
        }

        if (clockModel == RANDOM_LOCAL_CLOCK) {
            if (localClockRateChangesStatistic == null) {
                localClockRateChangesStatistic = new Parameter("rateChanges", "number of random local clocks", true);
                localClockRateChangesStatistic.priorType = PriorType.POISSON_PRIOR;
                localClockRateChangesStatistic.poissonMean = 1.0;
                localClockRateChangesStatistic.poissonOffset = 0.0;
            }
            if (localClockRatesStatistic == null) {
                localClockRatesStatistic = new Parameter(LOCAL_CLOCK + "." + "rates", "random local clock rates", false);

                localClockRatesStatistic.priorType = PriorType.GAMMA_PRIOR;
                localClockRatesStatistic.gammaAlpha = 0.5;
                localClockRatesStatistic.gammaBeta = 2.0;
            }
            params.add(localClockRatesStatistic);
            params.add(localClockRateChangesStatistic);
        }

        if (clockModel != STRICT_CLOCK) {
            params.add(getParameter("meanRate"));
            params.add(getParameter(RateStatisticParser.COEFFICIENT_OF_VARIATION));
            params.add(getParameter("covariance"));
        }
    }

    protected Parameter getParameter(String name) {
        Parameter parameter = parameters.get(name);
        if (parameter == null) throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        return parameter;
    }

    /**
     * return a list of operators that are required
     *
     * @param ops the operator list
     */
    private void selectOperators(ArrayList<Operator> ops) {

        if (alignment != null) {
            switch (dataType) {
                case DataType.NUCLEOTIDES:

                    switch (nucSubstitutionModel) {
                        case HKY:
                            // if (frequencyPolicy == BeautiOptions.ESTIMATED || frequencyPolicy == BeautiOptions.EMPIRICAL){
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    ops.add(getOperator("hky" + i + ".kappa"));
                                }
                            } else {
                                ops.add(getOperator("hky.kappa"));
                            }
                            //}
                            if (frequencyPolicy == BeautiOptions.ESTIMATED) {
                                if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                    for (int i = 1; i <= partitionCount; i++) {
                                        ops.add(getOperator("hky" + i + ".frequencies"));
                                    }
                                } else {
                                    ops.add(getOperator("hky.frequencies"));
                                }
                            }
                            break;

                        case GTR:
                            //if (frequencyPolicy == BeautiOptions.ESTIMATED || frequencyPolicy == BeautiOptions.EMPIRICAL){
                            if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                for (int i = 1; i <= partitionCount; i++) {
                                    ops.add(getOperator("gtr" + i + ".ac"));
                                    ops.add(getOperator("gtr" + i + ".ag"));
                                    ops.add(getOperator("gtr" + i + ".at"));
                                    ops.add(getOperator("gtr" + i + ".cg"));
                                    ops.add(getOperator("gtr" + i + ".gt"));
                                }
                            } else {
                                ops.add(getOperator("gtr.ac"));
                                ops.add(getOperator("gtr.ag"));
                                ops.add(getOperator("gtr.at"));
                                ops.add(getOperator("gtr.cg"));
                                ops.add(getOperator("gtr.gt"));
                            }
                            //}

                            if (frequencyPolicy == BeautiOptions.ESTIMATED) {
                                if (partitionCount > 1 && unlinkedSubstitutionModel) {
                                    for (int i = 1; i <= partitionCount; i++) {
                                        ops.add(getOperator("gtr" + i + ".frequencies"));
                                    }
                                } else {
                                    ops.add(getOperator("gtr.frequencies"));
                                }
                            }
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown nucleotides substitution model");
                    }

                    break;

                case DataType.AMINO_ACIDS:
                    break;

                case DataType.TWO_STATES:
                case DataType.COVARION:
                    switch (binarySubstitutionModel) {
                        case BIN_SIMPLE:
                            break;

                        case BIN_COVARION:
                            ops.add(getOperator("bcov.alpha"));
                            ops.add(getOperator("bcov.s"));
                            ops.add(getOperator("bcov.frequencies"));
                            ops.add(getOperator("bcov.hfrequencies"));
                            break;

                        default:
                            throw new IllegalArgumentException("Unknown binary substitution model");
                    }
                    break;

                default:
                    throw new IllegalArgumentException("Unknown data type");
            }

            // if gamma do shape move
            if (gammaHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        ops.add(getOperator(SiteModel.SITE_MODEL + i + ".alpha"));
                    }
                } else {
                    ops.add(getOperator(SiteModel.SITE_MODEL + "." + "alpha"));
                }
            }
            // if pinv do pinv move
            if (invarHetero) {
                if (partitionCount > 1 && unlinkedHeterogeneityModel) {
                    for (int i = 1; i <= partitionCount; i++) {
                        ops.add(getOperator(SiteModel.SITE_MODEL + i + ".pInv"));
                    }
                } else {
                    ops.add(getOperator(SiteModel.SITE_MODEL + "." + "pInv"));
                }
            }

            if (partitionCount > 1) {
                if (!codonHeteroPattern.equals("112")) {
                    ops.add(getOperator("centeredMu"));
                }
                ops.add(getOperator("deltaMu"));
            }

            // if not fixed then do mutation rate move and up/down move
            if (!fixedSubstitutionRate) {
                if (clockModel == STRICT_CLOCK) {
                    ops.add(getOperator("clock.rate"));
                    ops.add(getOperator("upDownRateHeights"));
                } else if (clockModel == RANDOM_LOCAL_CLOCK) {
                    ops.add(getOperator("clock.rate"));
                    ops.add(getOperator("upDownRateHeights"));
                    ops.add(getOperator(LOCAL_CLOCK + "." + "rates"));
                    ops.add(getOperator(LOCAL_CLOCK + "." + "changes"));
                    ops.add(getOperator("treeBitMove"));
                } else {
                    if (clockModel == UNCORRELATED_EXPONENTIAL) {
                        ops.add(getOperator(UCED_MEAN));
                        ops.add(getOperator("upDownUCEDMeanHeights"));
                    } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                        ops.add(getOperator(UCLD_MEAN));
                        ops.add(getOperator(UCLD_STDEV));
                        ops.add(getOperator("upDownUCLDMeanHeights"));
                    } else {
                        throw new IllegalArgumentException("Unknown clock model");
                    }
//                    ops.add(getOperator("swapBranchRateCategories"));
                    ops.add(getOperator("randomWalkBranchRateCategories"));
                    ops.add(getOperator("unformBranchRateCategories"));
                }
            } else {
                if (clockModel == STRICT_CLOCK) {
                    // no parameter to operator on
                } else if (clockModel == RANDOM_LOCAL_CLOCK) {
                    ops.add(getOperator(LOCAL_CLOCK + "." + "rates"));
                    ops.add(getOperator(LOCAL_CLOCK + "." + "changes"));
                    ops.add(getOperator("treeBitMove"));
                } else {
                    if (clockModel == UNCORRELATED_EXPONENTIAL) {
                        // no parameter to operator on
                    } else if (clockModel == UNCORRELATED_LOGNORMAL) {
                        ops.add(getOperator(UCLD_STDEV));
                    } else {
                        throw new IllegalArgumentException("Unknown clock model");
                    }
 //                   ops.add(getOperator("swapBranchRateCategories"));
                    ops.add(getOperator("randomWalkBranchRateCategories"));
                    ops.add(getOperator("unformBranchRateCategories"));
                }
            }
        }

        if (nodeHeightPrior == CONSTANT) {
            ops.add(getOperator("constant.popSize"));
        } else if (nodeHeightPrior == EXPONENTIAL) {
            ops.add(getOperator("exponential.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("exponential.growthRate"));
            } else {
                ops.add(getOperator("exponential.doublingTime"));
            }
        } else if (nodeHeightPrior == LOGISTIC) {
            ops.add(getOperator("logistic.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("logistic.growthRate"));
            } else {
                ops.add(getOperator("logistic.doublingTime"));
            }
            ops.add(getOperator("logistic.t50"));
        } else if (nodeHeightPrior == EXPANSION) {
            ops.add(getOperator("expansion.popSize"));
            if (parameterization == GROWTH_RATE) {
                ops.add(getOperator("expansion.growthRate"));
            } else {
                ops.add(getOperator("expansion.doublingTime"));
            }
            ops.add(getOperator("expansion.ancestralProportion"));
        } else if (nodeHeightPrior == SKYLINE) {
            ops.add(getOperator("skyline.popSize"));
            ops.add(getOperator("skyline.groupSize"));
        } else if (nodeHeightPrior == EXTENDED_SKYLINE) {
            ops.add(getOperator("demographic.populationMean"));
            ops.add(getOperator("demographic.popSize"));
            ops.add(getOperator("demographic.indicators"));
            ops.add(getOperator("demographic.scaleActive"));
        } else if (nodeHeightPrior == YULE) {
            ops.add(getOperator("yule.birthRate"));
        } else if (nodeHeightPrior == BIRTH_DEATH) {
//            ops.add(getOperator("birthDeath.birthRate"));
//            ops.add(getOperator("birthDeath.deathRate"));
            ops.add(getOperator(BirthDeathModelParser.MEAN_GROWTH_RATE_PARAM_NAME));
            ops.add(getOperator(BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME));
            // at present we are not allowing the sampling of samplingProportion
        }

        ops.add(getOperator("treeModel.rootHeight"));
        ops.add(getOperator("uniformHeights"));

        // if not a fixed tree then sample tree space
        if (!fixedTree) {
            ops.add(getOperator("subtreeSlide"));
            ops.add(getOperator("narrowExchange"));
            ops.add(getOperator("wideExchange"));
            ops.add(getOperator("wilsonBalding"));
        }
    }

    private Operator getOperator(String name) {
        Operator operator = operators.get(name);
        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");
        return operator;
    }

    /**
     * Read options from a file
     *
     * @param includeData include a data block?
     * @param guessDates  guess dates?
     * @return the Document
     */
    public Document create(boolean includeData, boolean guessDates) {

        Element root = new Element("beauti");
        root.setAttribute("version", version);

        Element dataElement = new Element("data");

        //dataElement.addContent(createChild("fileNameStem", fileNameStem));

        dataElement.addContent(createChild("datesUnits", datesUnits));
        dataElement.addContent(createChild("datesDirection", datesDirection));
        dataElement.addContent(createChild("translation", translation));
        dataElement.addContent(createChild("userTree", userTree));

        if (includeData && originalAlignment != null) {
            Element alignmentElement = new Element(AlignmentParser.ALIGNMENT);
            alignmentElement.addContent(createChild("dataType", originalAlignment.getDataType().getType()));
            for (int i = 0; i < originalAlignment.getTaxonCount(); i++) {
                Element taxonElement = new Element(TaxonParser.TAXON);
                taxonElement.addContent(createChild(XMLParser.ID, originalAlignment.getTaxonId(i)));
                dr.evolution.util.Date date = originalAlignment.getTaxon(i).getDate();
                if (date != null) {
                    taxonElement.addContent(createChild("date", date.getTimeValue()));
                }
                Sequence sequence = originalAlignment.getSequence(i);
                taxonElement.addContent(createChild("sequence", sequence.getSequenceString()));
                alignmentElement.addContent(taxonElement);
            }
            dataElement.addContent(alignmentElement);
        }

        dataElement.addContent(createChild("guessDates", guessDates));
        dataElement.addContent(createChild("guessDateFromOrder", guessDateFromOrder));
        dataElement.addContent(createChild("fromLast", fromLast));
        dataElement.addContent(createChild("order", order));
        dataElement.addContent(createChild("prefix", prefix));
        dataElement.addContent(createChild("offset", offset));
        dataElement.addContent(createChild("unlessLessThan", unlessLessThan));
        dataElement.addContent(createChild("offset2", offset2));

        root.addContent(dataElement);

        Element taxaElement = new Element(TaxaParser.TAXA);

        for (Taxa taxonSet : taxonSets) {
            Element taxonSetElement = new Element("taxonSet");
            taxonSetElement.addContent(createChild(XMLParser.ID, taxonSet.getId()));
            taxonSetElement.addContent(createChild("enforceMonophyly",
                    taxonSetsMono.get(taxonSet) ? "true" : "false"));
            for (int j = 0; j < taxonSet.getTaxonCount(); j++) {
                Element taxonElement = new Element(TaxonParser.TAXON);
                taxonElement.addContent(createChild(XMLParser.ID, taxonSet.getTaxon(j).getId()));
                taxonSetElement.addContent(taxonElement);
            }
            taxaElement.addContent(taxonSetElement);
        }

        root.addContent(taxaElement);

        Element modelElement = new Element("model");

        modelElement.addContent(createChild("nucSubstitutionModel", nucSubstitutionModel));
        modelElement.addContent(createChild("aaSubstitutionModel", aaSubstitutionModel));
        modelElement.addContent(createChild("binarySubstitutionModel", binarySubstitutionModel));
        modelElement.addContent(createChild("frequencyPolicy", frequencyPolicy));
        modelElement.addContent(createChild("gammaHetero", gammaHetero));
        modelElement.addContent(createChild("gammaCategories", gammaCategories));
        modelElement.addContent(createChild("invarHetero", invarHetero));
        modelElement.addContent(createChild("codonHeteroPattern", codonHeteroPattern));
        modelElement.addContent(createChild("maximumTipHeight", maximumTipHeight));
        modelElement.addContent(createChild("hasSetFixedSubstitutionRate", hasSetFixedSubstitutionRate));
        modelElement.addContent(createChild("meanSubstitutionRate", meanSubstitutionRate));
        modelElement.addContent(createChild("fixedSubstitutionRate", fixedSubstitutionRate));
        modelElement.addContent(createChild("unlinkedSubstitutionModel", unlinkedSubstitutionModel));
        modelElement.addContent(createChild("unlinkedHeterogeneityModel", unlinkedHeterogeneityModel));
        modelElement.addContent(createChild("unlinkedFrequencyModel", unlinkedFrequencyModel));
        modelElement.addContent(createChild("clockModel", clockModel));
        modelElement.addContent(createChild("nodeHeightPrior", nodeHeightPrior));
        modelElement.addContent(createChild("parameterization", parameterization));
        modelElement.addContent(createChild("skylineGroupCount", skylineGroupCount));
        modelElement.addContent(createChild("skylineModel", skylineModel));
        modelElement.addContent(createChild("fixedTree", fixedTree));

        root.addContent(modelElement);

        Element priorsElement = new Element("priors");

        for (String name : parameters.keySet()) {
            Parameter parameter = parameters.get(name);
            Element e = new Element(name);
            e.addContent(createChild("initial", parameter.initial));
            e.addContent(createChild("priorType", parameter.priorType));
            e.addContent(createChild("priorEdited", parameter.priorEdited));
            e.addContent(createChild("uniformLower", parameter.uniformLower));
            e.addContent(createChild("uniformUpper", parameter.uniformUpper));
            e.addContent(createChild("exponentialMean", parameter.exponentialMean));
            e.addContent(createChild("exponentialOffset", parameter.exponentialOffset));
            e.addContent(createChild("normalMean", parameter.normalMean));
            e.addContent(createChild("normalStdev", parameter.normalStdev));
            e.addContent(createChild("logNormalMean", parameter.logNormalMean));
            e.addContent(createChild("logNormalStdev", parameter.logNormalStdev));
            e.addContent(createChild("logNormalOffset", parameter.logNormalOffset));
            e.addContent(createChild("gammaAlpha", parameter.gammaAlpha));
            e.addContent(createChild("gammaBeta", parameter.gammaBeta));
            e.addContent(createChild("gammaOffset", parameter.gammaOffset));
            priorsElement.addContent(e);
        }

        for (Taxa taxonSet : taxonSets) {
            Parameter statistic = statistics.get(taxonSet);
            Element e = new Element(statistic.getXMLName());
            e.addContent(createChild("initial", statistic.initial));
            e.addContent(createChild("priorType", statistic.priorType));
            e.addContent(createChild("priorEdited", statistic.priorEdited));
            e.addContent(createChild("uniformLower", statistic.uniformLower));
            e.addContent(createChild("uniformUpper", statistic.uniformUpper));
            e.addContent(createChild("exponentialMean", statistic.exponentialMean));
            e.addContent(createChild("exponentialOffset", statistic.exponentialOffset));
            e.addContent(createChild("normalMean", statistic.normalMean));
            e.addContent(createChild("normalStdev", statistic.normalStdev));
            e.addContent(createChild("logNormalMean", statistic.logNormalMean));
            e.addContent(createChild("logNormalStdev", statistic.logNormalStdev));
            e.addContent(createChild("logNormalOffset", statistic.logNormalOffset));
            e.addContent(createChild("gammaAlpha", statistic.gammaAlpha));
            e.addContent(createChild("gammaBeta", statistic.gammaBeta));
            e.addContent(createChild("gammaOffset", statistic.gammaOffset));
            priorsElement.addContent(e);
        }

        root.addContent(priorsElement);

        Element operatorsElement = new Element("operators");

        operatorsElement.addContent(createChild("autoOptimize", autoOptimize));
        for (String name : operators.keySet()) {
            Operator operator = operators.get(name);
            Element e = new Element(name);
            e.addContent(createChild("tuning", operator.tuning));
            e.addContent(createChild("tuningEdited", operator.tuningEdited));
            e.addContent(createChild("weight", operator.weight));
            e.addContent(createChild("inUse", operator.inUse));
            operatorsElement.addContent(e);
        }

        root.addContent(operatorsElement);

        Element mcmcElement = new Element("mcmc");

        mcmcElement.addContent(createChild("upgmaStartingTree", upgmaStartingTree));
        mcmcElement.addContent(createChild("chainLength", chainLength));
        mcmcElement.addContent(createChild("logEvery", logEvery));
        mcmcElement.addContent(createChild("echoEvery", echoEvery));
        //if (logFileName != null) mcmcElement.addContent(createChild("logFileName", logFileName));
        //if (treeFileName != null) mcmcElement.addContent(createChild("treeFileName", treeFileName));
        //mcmcElement.addContent(createChild("mapTreeLog", mapTreeLog));
        //if (mapTreeFileName != null) mcmcElement.addContent(createChild("mapTreeFileName", mapTreeFileName));
        mcmcElement.addContent(createChild("substTreeLog", substTreeLog));
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

        Element dataElement = root.getChild("data");
        Element taxaElement = root.getChild(TaxaParser.TAXA);
        Element modelElement = root.getChild("model");
        Element priorsElement = root.getChild("priors");
        Element operatorsElement = root.getChild("operators");
        Element mcmcElement = root.getChild("mcmc");

        if (dataElement != null) {
            //fileNameStem = getStringChild(dataElement, "fileNameStem", "untitled");

            datesUnits = getIntegerChild(dataElement, "datesUnits", YEARS);
            datesDirection = getIntegerChild(dataElement, "datesDirection", FORWARDS);
            translation = getIntegerChild(dataElement, "translation", NONE);
            userTree = getBooleanChild(dataElement, "userTree", false);

            Units.Type theUnits = Units.Type.SUBSTITUTIONS;
            if (datesUnits == YEARS) theUnits = Units.Type.YEARS;
            if (datesUnits == MONTHS) theUnits = Units.Type.MONTHS;
            if (datesUnits == DAYS) theUnits = Units.Type.DAYS;

            Element alignmentElement = dataElement.getChild(AlignmentParser.ALIGNMENT);
            if (alignmentElement != null) {
                originalAlignment = new SimpleAlignment();

                int dataType = getIntegerChild(alignmentElement, "dataType", DataType.NUCLEOTIDES);
                switch (dataType) {
                    case DataType.NUCLEOTIDES:
                        originalAlignment.setDataType(Nucleotides.INSTANCE);
                        break;
                    case DataType.AMINO_ACIDS:
                        originalAlignment.setDataType(AminoAcids.INSTANCE);
                        break;
                    case DataType.TWO_STATES:
                        originalAlignment.setDataType(TwoStates.INSTANCE);
                        break;
                    default:
                        originalAlignment.setDataType(Nucleotides.INSTANCE);
                }

                for (Object o : alignmentElement.getChildren(TaxonParser.TAXON)) {
                    Element taxonElement = (Element) o;

                    String id = getStringChild(taxonElement, XMLParser.ID, "");
                    Taxon taxon = new Taxon(id);

                    if (taxonElement.getChild("date") != null) {
                        double dateValue = getDoubleChild(taxonElement, "date", 0.0);

                        if (datesDirection == FORWARDS) {
                            taxon.setDate(Date.createTimeSinceOrigin(dateValue, theUnits, 0.0));
                        } else {
                            taxon.setDate(Date.createTimeAgoFromOrigin(dateValue, theUnits, 0.0));
                        }
                    }
                    String seqString = getStringChild(taxonElement, "sequence", "");
                    Sequence sequence = new Sequence(taxon, seqString);

                    originalAlignment.addSequence(sequence);
                }
                taxonList = originalAlignment;
                alignment = originalAlignment;
            }

            guessDates = getBooleanChild(dataElement, "guessDates", false);
            guessDateFromOrder = getBooleanChild(dataElement, "guessDateFromOrder", false);
            fromLast = getBooleanChild(dataElement, "fromLast", false);
            order = getIntegerChild(dataElement, "order", 0);
            prefix = getStringChild(dataElement, "prefix", "");
            offset = getDoubleChild(dataElement, "offset", 0);
            unlessLessThan = getDoubleChild(dataElement, "unlessLessThan", 0);
            offset2 = getDoubleChild(dataElement, "offset2", 0);
        }

        if (taxaElement != null) {
            for (Object ts : taxaElement.getChildren("taxonSet")) {
                Element taxonSetElement = (Element) ts;

                String id = getStringChild(taxonSetElement, XMLParser.ID, "");
                final Taxa taxonSet = new Taxa(id);

                Boolean enforceMonophyly = Boolean.valueOf(getStringChild(taxonSetElement, "enforceMonophyly", "false"));
                for (Object o : taxonSetElement.getChildren(TaxonParser.TAXON)) {
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
            if (nodeHeightPrior == NONE) nodeHeightPrior = YULE;

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
                    statistic = new Parameter(taxonSet, "tMRCA for taxon set ");
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
        }
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

    public void guessDates() {

        for (int i = 0; i < originalAlignment.getTaxonCount(); i++) {
            java.util.Date origin = new java.util.Date(0);

            double d = 0.0;

            try {
                if (guessDateFromOrder) {
                    d = guessDateFromOrder(originalAlignment.getTaxonId(i), order, fromLast);
                } else {
                    d = guessDateFromPrefix(originalAlignment.getTaxonId(i), prefix);
                }

            } catch (GuessDatesException gfe) {
                //
            }

            if (offset > 0) {
                if (unlessLessThan > 0) {
                    if (d < unlessLessThan) {
                        d += offset2;
                    } else {
                        d += offset;
                    }
                } else {
                    d += offset;
                }
            }

            Date date = Date.createTimeSinceOrigin(d, Units.Type.YEARS, origin);
            originalAlignment.getTaxon(i).setAttribute("date", date);
        }

        // adjust the dates to the current timescale...
        timeScaleChanged();
    }

    public double guessDateFromOrder(String label, int order, boolean fromLast) throws GuessDatesException {

        String field;

        if (fromLast) {
            int count = 0;
            int i = label.length() - 1;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c) && c != '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                if (i < 0) throw new GuessDatesException("Missing number field in taxon label, " + label);

                int j = i + 1;

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i--;
                    if (i < 0) break;
                    c = label.charAt(i);
                }

                field = label.substring(i + 1, j);

                count++;

            } while (count <= order);

        } else {
            int count = 0;
            int i = 0;

            char c = label.charAt(i);

            do {
                // first find a part of a number
                while (!Character.isDigit(c) && c != '.') {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }
                int j = i;

                if (i == label.length()) throw new GuessDatesException("Missing number field in taxon label, " + label);

                // now find the beginning of the number
                while (Character.isDigit(c) || c == '.') {
                    i++;
                    if (i == label.length()) break;
                    c = label.charAt(i);
                }

                field = label.substring(j, i);

                count++;

            } while (count <= order);
        }

        return Double.parseDouble(field);
    }

    public double guessDateFromPrefix(String label, String prefix) throws GuessDatesException {

        int i = label.indexOf(prefix);

        if (i == -1) throw new GuessDatesException("Missing prefix in taxon label, " + label);

        i += prefix.length();
        int j = i;

        // now find the beginning of the number
        char c = label.charAt(i);
        while (i < label.length() - 1 && (Character.isDigit(c) || c == '.')) {
            i++;
            c = label.charAt(i);
        }

        if (i == j) throw new GuessDatesException("Missing field after prefix in taxon label, " + label);

        String field = label.substring(j, i + 1);

        double d;

        try {
            d = Double.parseDouble(field);
        } catch (NumberFormatException nfe) {
            throw new GuessDatesException("Badly formated date in taxon label, " + label);
        }

        return d;
    }

    private void timeScaleChanged() {

        for (int i = 0; i < alignment.getTaxonCount(); i++) {
            Date date = alignment.getTaxon(i).getDate();
            double d = date.getTimeValue();

            Date newDate = createDate(d, units, datesDirection == BACKWARDS, 0.0);

            alignment.getTaxon(i).setDate(newDate);
        }

    }

    private Date createDate(double timeValue, Units.Type units, boolean backwards, double origin) {
        if (backwards) {
            return Date.createTimeAgoFromOrigin(timeValue, units, origin);
        } else {
            return Date.createTimeSinceOrigin(timeValue, units, origin);
        }
    }

    public class Parameter {

        /**
         * A constructor for "special" parameters which are not user-configurable
         *
         * @param name        the name
         * @param description the description
         */
        public Parameter(String name, String description) {
            this.name = name;
            this.description = description;
            this.scale = NONE;
            this.isNodeHeight = false;
            this.isStatistic = false;
            this.taxa = null;
            this.priorType = PriorType.NONE;
            this.initial = Double.NaN;
            this.lower = Double.NaN;
            this.upper = Double.NaN;
        }

        public Parameter(String name, String description, int scale,
                         double initial, double lower, double upper) {
            this.name = name;
            this.description = description;
            this.initial = initial;
            this.isNodeHeight = false;
            this.isStatistic = false;

            this.taxa = null;

            this.priorType = PriorType.UNIFORM_PRIOR;
            this.scale = scale;
            this.priorEdited = false;
            this.lower = lower;
            this.upper = upper;

            uniformLower = lower;
            uniformUpper = upper;
        }

        public Parameter(TaxonList taxa, String description) {
            this.taxa = taxa;
            this.name = null;
            this.description = description;

            this.isNodeHeight = true;
            this.isStatistic = true;
            this.priorType = PriorType.NONE;
            this.scale = TIME_SCALE;
            this.priorEdited = false;
            this.lower = 0.0;
            this.upper = Double.MAX_VALUE;

            uniformLower = lower;
            uniformUpper = upper;
        }

        public Parameter(String name, String description, boolean isDiscrete) {
            this.taxa = null;

            this.name = name;
            this.description = description;

            this.isNodeHeight = false;
            this.isStatistic = true;
            this.isDiscrete = isDiscrete;
            this.priorType = PriorType.UNIFORM_PRIOR;
            this.scale = NONE;
            this.priorEdited = false;
            this.initial = Double.NaN;
            this.lower = Double.NaN;
            this.upper = Double.NaN;
        }

        public Parameter(String name, String description, double lower, double upper) {
            this.taxa = null;

            this.name = name;
            this.description = description;

            this.isNodeHeight = false;
            this.isStatistic = true;
            this.isDiscrete = false;
            this.priorType = PriorType.UNIFORM_PRIOR;
            this.scale = NONE;
            this.priorEdited = false;
            this.initial = Double.NaN;
            this.lower = lower;
            this.upper = upper;

            uniformLower = lower;
            uniformUpper = upper;
        }

        public Parameter(String name, String description, boolean isNodeHeight,
                         double initial, double lower, double upper) {
            this.name = name;
            this.description = description;
            this.initial = initial;

            this.taxa = null;

            this.isNodeHeight = isNodeHeight;
            this.isStatistic = false;
            this.priorType = PriorType.NONE;
            this.scale = TIME_SCALE;
            this.priorEdited = false;
            this.lower = lower;
            this.upper = upper;

            uniformLower = lower;
            uniformUpper = upper;
        }

        public String getName() {
            if (taxa != null) {
                return "tmrca(" + taxa.getId() + ")";
            } else {
                return name;
            }
        }

        public String getXMLName() {
            if (taxa != null) {
                return "tmrca_" + taxa.getId();
            } else {
                return name;
            }
        }

        public String getDescription() {
            if (taxa != null) {
                return description + taxa.getId();
            } else {
                return description;
            }
        }

        private final String name;
        private final String description;
        public double initial;

        public final TaxonList taxa;

        public boolean isDiscrete = false;

        public boolean isFixed = false;
        public final boolean isNodeHeight;
        public final boolean isStatistic;

        public PriorType priorType;
        public boolean priorEdited;
        public final int scale;
        public double lower;
        public double upper;

        public double uniformUpper = 0.0;
        public double uniformLower = 0.0;
        public double exponentialMean = 1.0;
        public double exponentialOffset = 0.0;
        public double normalMean = 1.0;
        public double normalStdev = 1.0;
        public double logNormalMean = 0.0;
        public double logNormalStdev = 1.0;
        public double logNormalOffset = 0.0;
        public double gammaAlpha = 1.0;
        public double gammaBeta = 1.0;
        public double gammaOffset = 0.0;
        public double poissonMean = 1.0;
        public double poissonOffset = 0.0;
    }

    public class Operator {
        public Operator(String name, String description, Parameter parameter, String operatorType, double tuning, double weight) {
            this.name = name;
            this.description = description;
            this.parameter1 = parameter;
            this.parameter2 = null;

            this.type = operatorType;
            this.tuningEdited = false;
            this.tuning = tuning;
            this.weight = weight;

            this.inUse = true;
        }

        public Operator(String name, String description,
                        Parameter parameter1, Parameter parameter2,
                        String operatorType, double tuning, double weight) {
            this.name = name;
            this.description = description;
            this.parameter1 = parameter1;
            this.parameter2 = parameter2;

            this.type = operatorType;
            this.tuningEdited = false;
            this.tuning = tuning;
            this.weight = weight;

            this.inUse = true;
        }

        public String getDescription() {
            if (description == null || description.length() == 0) {
                String prefix = "";
                if (type.equals(SCALE)) {
                    prefix = "Scales the ";
                } else if (type.equals(RANDOM_WALK)) {
                    prefix = "A random-walk on the ";
                }
                return prefix + parameter1.getDescription();
            } else {
                return description;
            }
        }

        public boolean isTunable() {
            return tuning > 0;
        }

        public final String name;
        public final String description;

        public final String type;
        public boolean tuningEdited;
        public double tuning;
        public double weight;
        public boolean inUse;

        public final Parameter parameter1;
        public final Parameter parameter2;

    }

    public static final String version = "1.4";
    public static final int YEARS = 0;
    public static final int MONTHS = 1;
    public static final int DAYS = 2;
    public static final int FORWARDS = 0;
    public static final int BACKWARDS = 1;
    public static final int NONE = -1;

    public static final int JC = 0;
    public static final int HKY = 1;
    public static final int GTR = 2;

    public static final int BLOSUM_62 = 0;
    public static final int DAYHOFF = 1;
    public static final int JTT = 2;
    public static final int MT_REV_24 = 3;
    public static final int CP_REV_45 = 4;
    public static final int WAG = 5;

    public static final int BIN_SIMPLE = 0;
    public static final int BIN_COVARION = 1;

    public static final int ESTIMATED = 0;
    public static final int EMPIRICAL = 1;
    public static final int ALLEQUAL = 2;

    public static final int CONSTANT = 0;
    public static final int EXPONENTIAL = 1;
    public static final int LOGISTIC = 2;
    public static final int EXPANSION = 3;
    public static final int SKYLINE = 4;
    public static final int EXTENDED_SKYLINE = 5;
    public static final int YULE = 6;
    public static final int BIRTH_DEATH = 7;

    public static final int STRICT_CLOCK = 0;
    public static final int UNCORRELATED_EXPONENTIAL = 1;
    public static final int UNCORRELATED_LOGNORMAL = 2;
    public static final int RANDOM_LOCAL_CLOCK = 3;

    public static final int GROWTH_RATE = 0;
    public static final int DOUBLING_TIME = 1;
    public static final int CONSTANT_SKYLINE = 0;
    public static final int LINEAR_SKYLINE = 1;

    public static final int TIME_SCALE = 0;
    public static final int GROWTH_RATE_SCALE = 1;
    public static final int BIRTH_RATE_SCALE = 2;
    public static final int SUBSTITUTION_RATE_SCALE = 3;
    public static final int LOG_STDEV_SCALE = 4;
    public static final int SUBSTITUTION_PARAMETER_SCALE = 5;
    public static final int T50_SCALE = 6;
    public static final int UNITY_SCALE = 7;

    public static final String SCALE = "scale";
    public static final String RANDOM_WALK = "randomWalk";
    public static final String INTEGER_RANDOM_WALK = "integerRandomWalk";
    public static final String UP_DOWN = "upDown";
    public static final String SCALE_ALL = "scaleAll";
    public static final String CENTERED_SCALE = "centeredScale";
    public static final String DELTA_EXCHANGE = "deltaExchange";
    public static final String INTEGER_DELTA_EXCHANGE = "integerDeltaExchange";
    public static final String SWAP = "swap";
    public static final String BITFLIP = "bitFlip";
    public static final String TREE_BIT_MOVE = "treeBitMove";
    public static final String SAMPLE_NONACTIVE = "sampleNoneActiveOperator";
    public static final String SCALE_WITH_INDICATORS = "scaleWithIndicators";

    public static final String UNIFORM = "uniform";
    public static final String INTEGER_UNIFORM = "integerUniform";
    public static final String SUBTREE_SLIDE = "subtreeSlide";
    public static final String NARROW_EXCHANGE = "narrowExchange";
    public static final String WIDE_EXCHANGE = "wideExchange";
    public static final String WILSON_BALDING = "wilsonBalding";
    public String fileNameStem = "untitled";
    public String logFileName = null;
    public String treeFileName = null;
    public boolean mapTreeLog = false;
    public String mapTreeFileName = null;
    public boolean substTreeLog = false;
    public String substTreeFileName = null;

    // Data options
    public int dataType = DataType.NUCLEOTIDES;

    public TaxonList taxonList = null;
    public SimpleAlignment originalAlignment = null;
    public List<Taxa> taxonSets = new ArrayList<Taxa>();
    public Map<Taxa, Boolean> taxonSetsMono = new HashMap<Taxa, Boolean>();
    public Alignment alignment = null;
    public Tree tree = null;
    public boolean alignmentReset = true;
    public double meanDistance = 1.0;
    public int datesUnits = YEARS;
    public int datesDirection = FORWARDS;
    public double maximumTipHeight = 0.0;
    public int translation = 0;
    public boolean userTree = false;

    public boolean guessDates = false;
    public boolean guessDateFromOrder = true;
    public boolean fromLast = false;
    public int order = 0;
    public String prefix;
    public double offset = 0.0;
    public double unlessLessThan = 0.0;
    public double offset2 = 0.0;

    // Model options
    public int partitionCount = 1;
    public int nucSubstitutionModel = HKY;
    public int aaSubstitutionModel = BLOSUM_62;
    public int binarySubstitutionModel = BIN_SIMPLE;

    public int frequencyPolicy = ESTIMATED;

    public boolean gammaHetero = false;
    public int gammaCategories = 4;
    public boolean invarHetero = false;
    public String codonHeteroPattern = null;
    public double meanSubstitutionRate = 1.0;
    public boolean unlinkedSubstitutionModel = false;
    public boolean unlinkedHeterogeneityModel = false;
    public boolean unlinkedFrequencyModel = false;
    public int nodeHeightPrior = CONSTANT;
    public int parameterization = GROWTH_RATE;
    public int skylineGroupCount = 10;
    public int skylineModel = CONSTANT_SKYLINE;
    public String extendedSkylineModel = VariableDemographicModel.Type.LINEAR.toString();
    public double birthDeathSamplingProportion = 1.0;
    public boolean fixedTree = false;
    public Units.Type units = Units.Type.SUBSTITUTIONS;
    public boolean fixedSubstitutionRate = false;
    public boolean hasSetFixedSubstitutionRate = false;
    public int clockModel = STRICT_CLOCK;

    // MCMC options
    public boolean upgmaStartingTree = false;
    public int chainLength = 10000000;
    public int logEvery = 1000;
    public int echoEvery = 1000;
    public int burnIn = 100000;
    public String fileName = null;
    public boolean autoOptimize = true;
    public boolean performTraceAnalysis = false;
    public boolean generateCSV = true;  // until/if a button
    public boolean samplePriorOnly = false;

    public HashMap<String, Parameter> parameters = new HashMap<String, Parameter>();
    public HashMap<TaxonList, Parameter> statistics = new HashMap<TaxonList, Parameter>();
    public HashMap<String, Operator> operators = new HashMap<String, Operator>();
    public Parameter localClockRateChangesStatistic = null;
    public Parameter localClockRatesStatistic = null;
}
