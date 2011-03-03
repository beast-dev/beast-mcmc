package test.dr.calibration;


import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxa;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.CoalescentSimulator;
import dr.evomodel.coalescent.ConstantPopulationModel;
import dr.evomodel.coalescent.DemographicModel;
import dr.evomodel.operators.SubtreeSlideOperator;
import dr.evomodel.speciation.BirthDeathGernhard08Model;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciationModel;
import dr.evomodel.tree.*;
import dr.evomodelxml.coalescent.ConstantPopulationModelParser;
import dr.evomodelxml.tree.TreeModelParser;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.loggers.ArrayLogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.mcmc.MCMCOptions;
import dr.inference.model.BooleanLikelihood;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.*;
import dr.inference.trace.ArrayTraceList;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inferencexml.model.CompoundLikelihoodParser;
import dr.math.distributions.LogNormalDistribution;
import dr.util.NumberFormatter;
import dr.xml.XMLParseException;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.List;

/**
 * @auther Walter Xie
 */

public class TestCalibratedYuleModel {
    protected static final String TL = "TL";
    protected static final String TREE_HEIGHT = TreeModel.TREE_MODEL + "." + TreeModelParser.ROOT_HEIGHT;
    //    private final int treeSize;
    private final BufferedWriter out;
    Taxa taxa;

    public TestCalibratedYuleModel(int treeSize, double S, BufferedWriter out) throws Exception {
//        this.treeSize = treeSize;
        this.out = out;
        out.write(Integer.toString(treeSize) + "\t");

        TreeModel treeModel = createTreeModel(treeSize);

        Parameter brParameter = new Parameter.Default("birthRate", 2.0, 0.0, 100.0);

        OperatorSchedule schedule = new SimpleOperatorSchedule();
        MCMCOperator operator = new SubtreeSlideOperator(treeModel, 10, 1, true, false, false, false, CoercionMode.COERCION_ON);
        schedule.addOperator(operator);

        operator = new ScaleOperator(brParameter, 0.5);
        operator.setWeight(1.0);
        schedule.addOperator(operator);

        yuleTester(treeModel, schedule, brParameter, S);

    }


    protected TreeModel createTreeModel(int treeSize) throws Exception {
        taxa = new Taxa();
        for (int i = 0; i < treeSize; i++) {
            taxa.addTaxon(new Taxon("T" + Integer.toString(i)));
        }
        Taxa taxaSubSet = new Taxa();
        for (int i = 0; i < treeSize / 2; i++) {
            taxaSubSet.addTaxon(new Taxon("T" + Integer.toString(i)));
        }
        System.out.println("taxaSubSet_size = " + taxaSubSet.getTaxonCount());

        Parameter popSize = new Parameter.Default(treeSize);
        popSize.setId(ConstantPopulationModelParser.POPULATION_SIZE);
        ConstantPopulationModel startingTree = new ConstantPopulationModel(popSize, Units.Type.YEARS);
//        ConstantPopulation constant = (ConstantPopulation) startingTree.getDemographicFunction();
//        CoalescentSimulator simulator = new CoalescentSimulator();
//        Tree tree = simulator.simulateTree(taxonList, constant);

        Tree tree = calibration(taxa, startingTree, taxaSubSet);

        return new TreeModel(tree);//treeModel
    }

    private Tree calibration(final TaxonList taxa, DemographicModel demoModel, TaxonList taxaSubSet) throws Exception {
        CoalescentSimulator simulator = new CoalescentSimulator();

//        DemographicModel demoModel = (DemographicModel) xo.getChild(DemographicModel.class);
        List<TaxonList> taxonLists = new ArrayList<TaxonList>();
        List<Tree> subtrees = new ArrayList<Tree>();

        double rootHeight = -1.0;

        // should have one child that is node
//        for (int i = 0; i < xo.getChildCount(); i++) {
//            final Object child = xo.getChild(i);

        // AER - swapped the order of these round because Trees are TaxonLists...
//            if (child instanceof Tree) {
//                subtrees.add((Tree) child);
//            } else if (child instanceof TaxonList) {
//                taxonLists.add((TaxonList) child);
//            } else if (xo.getChildName(i).equals(CONSTRAINED_TAXA)) {
//                rootHeight = -1; // ignore it? should we errror?

//                XMLObject constrainedTaxa = (XMLObject) child;

        // all taxa
//                final TaxonList taxa = (TaxonList) constrainedTaxa.getChild(TaxonList.class);

        List<dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint> constraints
                = new ArrayList<dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint>();
        final String setsNotCompatibleMessage = "taxa sets not compatible";

        // pick up all constraints. order in partial order, where taxa_1 @in taxa_2 implies
        // taxa_1 is before taxa_2.


//                for (int nc = 0; nc < constrainedTaxa.getChildCount(); ++nc) {
//
//                    final Object object = constrainedTaxa.getChild(nc);
//                    if (object instanceof XMLObject) {
//                        final XMLObject constraint = (XMLObject) object;
//
//                        if (constraint.getName().equals(TMRCA_CONSTRAINT)) {
//                            TaxonList taxaSubSet = (TaxonList) constraint.getChild(TaxonList.class);
//                            ParametricDistributionModel dist =
//                                    (ParametricDistributionModel) constraint.getChild(ParametricDistributionModel.class);
        boolean isMono = true;

        final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint taxaConstraint
                = new dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint(taxaSubSet, null, isMono);
        int insertPoint;
        for (insertPoint = 0; insertPoint < constraints.size(); ++insertPoint) {
            // if new <= constraints[insertPoint] insert before insertPoint

            final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint iConstraint = constraints.get(insertPoint);
            if (iConstraint.isMonophyletic) {
                if (!taxaConstraint.isMonophyletic) {
                    continue;
                }

                final TaxonList taxonsip = iConstraint.taxons;
                final int nIn = simulator.sizeOfIntersection(taxonsip, taxaSubSet);
                if (nIn == taxaSubSet.getTaxonCount()) {
                    break;
                }
                if (nIn > 0 && nIn != taxonsip.getTaxonCount()) {
                    throw new XMLParseException(setsNotCompatibleMessage);
                }
            } else {
                // reached non mono area
                if (!taxaConstraint.isMonophyletic) {
                    if (iConstraint.upper >= taxaConstraint.upper) {
                        break;
                    }
                } else {
                    break;
                }
            }
        }
        constraints.add(insertPoint, taxaConstraint);
//                        }
//                    }
//                }
        final int nConstraints = constraints.size();
        System.out.println("nConstraints = " + nConstraints);
        if (nConstraints == 0) {
            if (taxa != null) {
                taxonLists.add(taxa);
            }
        } else {
            for (int nc = 0; nc < nConstraints; ++nc) {
                dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                if (!cnc.isMonophyletic) {
                    for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                        dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint cnc1 = constraints.get(nc1);
                        int x = simulator.sizeOfIntersection(cnc.taxons, cnc1.taxons);
                        if (x > 0) {
                            Taxa combinedTaxa = new Taxa(cnc.taxons);
                            combinedTaxa.addTaxa(cnc1.taxons);
                            cnc = new dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint(combinedTaxa, cnc.lower, cnc.upper, cnc.isMonophyletic);
                            constraints.set(nc, cnc);
                        }
                    }
                }
            }
            // determine upper bound for each set.
            double[] upper = new double[nConstraints];
            for (int nc = nConstraints - 1; nc >= 0; --nc) {
                final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                if (cnc.realLimits()) {
                    upper[nc] = cnc.upper;
                } else {
                    upper[nc] = Double.POSITIVE_INFINITY;
                }
            }

            for (int nc = nConstraints - 1; nc >= 0; --nc) {
                final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint cnc = constraints.get(nc);
                if (upper[nc] < Double.POSITIVE_INFINITY) {
                    for (int nc1 = nc - 1; nc1 >= 0; --nc1) {
                        final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint cnc1 = constraints.get(nc1);
                        if (simulator.contained(cnc1.taxons, cnc.taxons)) {
                            upper[nc1] = Math.min(upper[nc1], upper[nc]);
                            if (cnc1.realLimits() && cnc1.lower > upper[nc1]) {
                                throw new XMLParseException(setsNotCompatibleMessage);
                            }
                            break;
                        }
                    }
                }
            }
            // collect subtrees here
            List<Tree> st = new ArrayList<Tree>();
            for (int nc = 0; nc < constraints.size(); ++nc) {
                final dr.evomodel.coalescent.CoalescentSimulator.TaxaConstraint nxt = constraints.get(nc);
                // collect all previously built subtrees which are a subset of taxa set to be added
                List<Tree> subs = new ArrayList<Tree>();
                Taxa newTaxons = new Taxa(nxt.taxons);
                for (int k = 0; k < st.size(); ++k) {
                    final Tree stk = st.get(k);
                    int x = simulator.sizeOfIntersection(stk, nxt.taxons);
                    if (x == st.get(k).getTaxonCount()) {
                        final Tree tree = st.remove(k);
                        --k;
                        subs.add(tree);
                        newTaxons.removeTaxa(tree);
                    }
                }

                SimpleTree tree = simulator.simulateTree(newTaxons, demoModel);
                final double lower = nxt.realLimits() ? nxt.lower : 0;
                if (upper[nc] < Double.MAX_VALUE) {
                    simulator.attemptToScaleTree(tree, (lower + upper[nc]) / 2);
                }
                if (subs.size() > 0) {
                    if (tree.getTaxonCount() > 0) subs.add(tree);
                    double h = -1;
                    if (upper[nc] < Double.MAX_VALUE) {
                        for (Tree t : subs) {
                            h = Math.max(h, t.getNodeHeight(t.getRoot()));
                        }
                        h = (h + upper[nc]) / 2;
                    }
                    tree = simulator.simulateTree(subs.toArray(new Tree[subs.size()]), demoModel, h, true);
                }
                st.add(tree);

            }

            // add a taxon list for remaining taxa
            if (taxa != null) {
                final Taxa list = new Taxa();
                for (int j = 0; j < taxa.getTaxonCount(); ++j) {
                    Taxon taxonj = taxa.getTaxon(j);
                    for (Tree aSt : st) {
                        if (aSt.getTaxonIndex(taxonj) >= 0) {
                            taxonj = null;
                            break;
                        }
                    }
                    if (taxonj != null) {
                        list.addTaxon(taxonj);
                    }
                }
                if (list.getTaxonCount() > 0) {
                    taxonLists.add(list);
                }
            }
            if (st.size() > 1) {
                final Tree t = simulator.simulateTree(st.toArray(new Tree[st.size()]), demoModel, -1, false);
                subtrees.add(t);
            } else {
                subtrees.add(st.get(0));
            }
        }
//            }
//        }
//        System.out.println("subtrees = " + subtrees);
        if (taxonLists.size() == 0) {
            if (subtrees.size() == 1) {
                return subtrees.get(0);
            }
            throw new Exception("Expected at least one taxonList or two subtrees");
        }

        try {
            Tree[] trees = new Tree[taxonLists.size() + subtrees.size()];
            // simulate each taxonList separately
            for (int i = 0; i < taxonLists.size(); i++) {
                trees[i] = simulator.simulateTree(taxonLists.get(i), demoModel);
                System.out.println("trees[" + i + "] = " + trees[i]);
            }
            // add the preset trees
            for (int i = 0; i < subtrees.size(); i++) {
                trees[i + taxonLists.size()] = subtrees.get(i);
                System.out.println("trees[" + (i + taxonLists.size()) + "] = " + trees[i + taxonLists.size()]);
            }
            System.out.println("taxonLists.size() = " + taxonLists.size() + ";  "
                    + "subtrees.size() = " + subtrees.size() + ";  rootHeight = " + rootHeight);
            return simulator.simulateTree(trees, demoModel, rootHeight, trees.length != 1);

        } catch (IllegalArgumentException iae) {
            throw new Exception(iae.getMessage());
        }
    }

    private void yuleTester(TreeModel treeModel, OperatorSchedule schedule, Parameter brParameter, double S)
            throws IOException, Tree.MissingTaxonException {

        MCMC mcmc = new MCMC("mcmc1");
        MCMCOptions options = new MCMCOptions();
        options.setChainLength(20000000);
        options.setUseCoercion(true);
        options.setCoercionDelay(options.getChainLength() / 100);
        options.setTemperature(1.0);
        options.setFullEvaluationCount(2000);

        TreeLengthStatistic tls = new TreeLengthStatistic(TL, treeModel);
        TreeHeightStatistic rootHeight = new TreeHeightStatistic(TREE_HEIGHT, treeModel);

        SpeciationModel speciationModel = new BirthDeathGernhard08Model("yule", brParameter, null, null,
                BirthDeathGernhard08Model.TreeType.UNSCALED, Units.Type.SUBSTITUTIONS, false);
        Likelihood speciationLikelihood = new SpeciationLikelihood(treeModel, speciationModel, "yule.like");

        TMRCAStatistic tmrca = new TMRCAStatistic("tmrca(halfTaxa)", treeModel, taxa, false, false);
        DistributionLikelihood logNormalLikelihood = new DistributionLikelihood(
                new LogNormalDistribution(1.0, S), 0); // meanInRealSpace="false"
        logNormalLikelihood.addData(tmrca);

        MonophylyStatistic monophylyStatistic = new MonophylyStatistic("monophyly(halfTaxa)", treeModel, taxa, null);
        BooleanLikelihood booleanLikelihood = new BooleanLikelihood();
        booleanLikelihood.addData(monophylyStatistic);

        //CompoundLikelihood
        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        likelihoods.add(speciationLikelihood);
        likelihoods.add(logNormalLikelihood);
        likelihoods.add(booleanLikelihood);
        Likelihood prior = new CompoundLikelihood(0, likelihoods);
        prior.setId(CompoundLikelihoodParser.PRIOR);

        ArrayLogFormatter logformatter = new ArrayLogFormatter(false);

        MCLogger[] loggers = new MCLogger[1];
        loggers[0] = new MCLogger(logformatter, options.getChainLength() / 10000, false);
        loggers[0].add(speciationLikelihood);
        loggers[0].add(rootHeight);
        loggers[0].add(tls);
        loggers[0].add(brParameter);

//        loggers[1] = new MCLogger(new TabDelimitedFormatter(System.out), options.getChainLength() / 100000, false);
//        loggers[1].add(speciationLikelihood);
//        loggers[1].add(rootHeight);
//        loggers[1].add(tls);

        mcmc.setShowOperatorAnalysis(false);

        mcmc.init(options, prior, schedule, loggers);
        mcmc.run();

        List<Trace> traces = logformatter.getTraces();
        ArrayTraceList traceList = new ArrayTraceList("yuleModelTest", traces, 1000);

        for (int i = 1; i < traces.size(); i++) {
            traceList.analyseTrace(i);
        }

        NumberFormatter formatter = new NumberFormatter(8);

        TraceCorrelation tlStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TL));
        TraceCorrelation treeHeightStats = traceList.getCorrelationStatistics(traceList.getTraceIndex(TREE_HEIGHT));
        out.write(formatter.format(treeHeightStats.getMean()));
        out.write("\t");

        double expectedRootHeight = Math.pow(Math.E, (1 + (Math.pow(0.2, 2) / 2)));
        out.write(formatter.format(expectedRootHeight));
        out.write("\t");

        double error = Math.abs((treeHeightStats.getMean() - expectedRootHeight) / expectedRootHeight);
        NumberFormat percentFormatter = NumberFormat.getNumberInstance();
        percentFormatter.setMinimumFractionDigits(5);
        percentFormatter.setMinimumFractionDigits(5);
        out.write(percentFormatter.format(error));
        out.write("\t");
        out.write(Double.toString(tlStats.getESS()));

        System.out.println("rootHeight = " + formatter.format(treeHeightStats.getMean())
                + ";  expectation = " + formatter.format(expectedRootHeight)
                + ";  error = " + percentFormatter.format(error)
                + ";  tl.ess = " + tlStats.getESS());

    }


    public static void main(String[] args) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter("TestCalibratedYuleModel.txt"));
            out.write("treeSize\trootHeight\texpectation\terror");
            out.newLine();

            int[] taxaSchedule = new int[]{4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128};
            double[] S_Schedule = new double[]{0.05, 0.1, 0.2, 0.4};

            for (double S : S_Schedule) {
                for (int i : taxaSchedule) {
                    System.out.print("treeSize = " + i + "\t" + "S = " + S + "\t");
                    out.write("treeSize = " + i + "\t" + "S = " + S + "\t");
                    TestCalibratedYuleModel testCalibratedYuleModel = new TestCalibratedYuleModel(i, S, out);
                    out.newLine();
                }
            }
            out.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
