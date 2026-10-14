/*
 * MascotMappingFixtures.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent.mascot;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.coalescent.basta.AbstractPopulationSizeModel;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.PiecewiseConstantPopulationSizeModel;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/** Small fixed trees and parameter values for testing MASCOT stochastic mapping, plus matching references. */
final class MascotMappingFixtures {

    private MascotMappingFixtures() {
    }

    static final class Spec {
        final String name;
        final int K;
        final String newick;
        final Map<String, String> tipCodes = new HashMap<>();   // taxon -> state code ("0", "1", or "X" etc.)
        final Map<String, String[]> ambiguities = new HashMap<>(); // code -> allowed states
        double[] epochTimes = new double[0];
        double[][] migration;                  // [epoch][K K], off-diagonal backward rates
        double[][] popSizes;                   // [epoch][K]
        double strictClock = Double.NaN;
        final Map<String, Double> cladeRates = new HashMap<>(); // clade key ("A" or "A,B") -> branch rate

        Spec(String name, int K, String newick) {
            this.name = name;
            this.K = K;
            this.newick = newick;
        }

        Spec tips(String... taxonEqualsCode) {
            for (String entry : taxonEqualsCode) {
                String[] parts = entry.split("=");
                tipCodes.put(parts[0], parts[1]);
            }
            return this;
        }

        Spec ambiguity(String code, String... states) {
            ambiguities.put(code, states);
            return this;
        }

        Spec epochs(double... times) {
            epochTimes = times;
            return this;
        }

        Spec migration(double[]... rates) {
            migration = rates;
            return this;
        }

        Spec popSizes(double[]... sizes) {
            popSizes = sizes;
            return this;
        }

        Spec strictClock(double rate) {
            strictClock = rate;
            return this;
        }

        Spec branchRate(String clade, double rate) {
            cladeRates.put(clade, rate);
            return this;
        }

        Built build(double maxStep) throws Exception {
            return new Built(this, maxStep);
        }
    }

    static final class Built {
        final Spec spec;
        final double maxStep;
        final TreeModel tree;
        final GeneralDataType dataType;
        final MascotLikelihood likelihood;
        final Parameter migrationRates;
        final Parameter popSizes;
        final double[][] tipPartials;   // by node number; null for internal nodes
        final double[] clockRates;      // by node number; null when there is no branch-rate model

        private Built(Spec spec, double maxStep) throws Exception {
            this.spec = spec;
            this.maxStep = maxStep;
            int K = spec.K;
            int epochCount = spec.epochTimes.length + 1;

            String[] codes = new String[K];
            for (int s = 0; s < K; s++) {
                codes[s] = Integer.toString(s);
            }
            dataType = new GeneralDataType(codes);
            for (Map.Entry<String, String[]> entry : spec.ambiguities.entrySet()) {
                dataType.addAmbiguity(entry.getKey(), entry.getValue());
            }

            tree = new DefaultTreeModel("treeModel", new NewickImporter(spec.newick).importTree(null));
            Patterns patterns = new Patterns(dataType, tree);
            int[] pattern = new int[tree.getTaxonCount()];
            for (int i = 0; i < pattern.length; i++) {
                pattern[i] = dataType.getState(spec.tipCodes.get(tree.getTaxonId(i)));
            }
            patterns.addPattern(pattern);

            tipPartials = new double[tree.getNodeCount()][];
            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                NodeRef tip = tree.getExternalNode(i);
                int state = dataType.getState(spec.tipCodes.get(tree.getNodeTaxon(tip).getId()));
                boolean[] allowed = dataType.getStateSet(state);
                double[] partial = new double[K];
                for (int s = 0; s < K; s++) {
                    partial[s] = allowed[s] ? 1.0 : 0.0;
                }
                tipPartials[tip.getNumber()] = partial;
            }

            double[] rates = new double[epochCount * K * (K - 1)];
            double[] sizes = new double[epochCount * K];
            int index = 0;
            for (int e = 0; e < epochCount; e++) {
                for (int a = 0; a < K; a++) {
                    for (int b = 0; b < K; b++) {
                        if (a != b) {
                            rates[index++] = spec.migration[e][a * K + b];
                        }
                    }
                    sizes[e * K + a] = spec.popSizes[e][a];
                }
            }
            migrationRates = new Parameter.Default("migrationRates", rates);
            popSizes = new Parameter.Default("popSizes", sizes);
            Parameter epochTimes = spec.epochTimes.length == 0 ? null :
                    new Parameter.Default("gridPoints", spec.epochTimes.clone());
            AbstractPopulationSizeModel populationSizeModel = epochTimes == null ?
                    new ConstantPopulationSizeModel("constant", popSizes, K, -1) :
                    new PiecewiseConstantPopulationSizeModel("piecewiseConstant", popSizes, epochTimes, K, -1);

            BranchRateModel branchRateModel = null;
            if (!Double.isNaN(spec.strictClock)) {
                branchRateModel = new StrictClockBranchRates(new Parameter.Default("clock", spec.strictClock));
            } else if (!spec.cladeRates.isEmpty()) {
                ArbitraryBranchRates arbitrary = new ArbitraryBranchRates(tree,
                        new Parameter.Default("branchRates", new double[tree.getNodeCount() - 1]),
                        new ArbitraryBranchRates.BranchRateTransform.None(), false);
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    NodeRef node = tree.getNode(i);
                    if (!tree.isRoot(node)) {
                        Double rate = spec.cladeRates.get(cladeKey(node));
                        arbitrary.setBranchRate(tree, node, rate == null ? 1.0 : rate);
                    }
                }
                branchRateModel = arbitrary;
            }
            if (branchRateModel != null) {
                clockRates = new double[tree.getNodeCount()];
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    NodeRef node = tree.getNode(i);
                    if (!tree.isRoot(node)) {
                        clockRates[node.getNumber()] = branchRateModel.getBranchRate(tree, node);
                    }
                }
            } else {
                clockRates = null;
            }

            likelihood = new MascotLikelihood(null, tree, patterns, migrationRates, null, branchRateModel,
                    populationSizeModel, K, epochTimes, maxStep, false,
                    MascotLikelihood.DEFAULT_ANCESTRAL_STATE_TAG_NAME, MascotLikelihood.DEFAULT_MODE_STATE_TAG_NAME);
        }

        MascotSurrogateOracle oracle() {
            return oracle(maxStep);
        }

        MascotSurrogateOracle oracle(double step) {
            return new MascotSurrogateOracle(tree, tipPartials, spec.K, spec.epochTimes, spec.migration,
                    spec.popSizes, clockRates, step);
        }

        ExactStructuredCoalescent exact() {
            int n = tree.getNodeCount();
            int[] parent = new int[n];
            int[] child1 = new int[n];
            int[] child2 = new int[n];
            double[] height = new double[n];
            String[] label = new String[n];
            double[] clock = new double[n];
            for (int i = 0; i < n; i++) {
                NodeRef node = tree.getNode(i);
                int v = node.getNumber();
                parent[v] = tree.isRoot(node) ? -1 : tree.getParent(node).getNumber();
                child1[v] = tree.isExternal(node) ? -1 : tree.getChild(node, 0).getNumber();
                child2[v] = tree.isExternal(node) ? -1 : tree.getChild(node, 1).getNumber();
                height[v] = tree.getNodeHeight(node);
                label[v] = tree.isExternal(node) ? tree.getNodeTaxon(node).getId() : null;
                clock[v] = clockRates == null ? 1.0 : clockRates[v];
            }
            double[] epochStarts = new double[spec.epochTimes.length + 1];
            System.arraycopy(spec.epochTimes, 0, epochStarts, 1, spec.epochTimes.length);
            ExactStructuredCoalescent.Model model = new ExactStructuredCoalescent.Model(spec.K, epochStarts,
                    spec.migration, spec.popSizes);
            return new ExactStructuredCoalescent(model, new ExactStructuredCoalescent.TreeSpec(parent, child1,
                    child2, height, tipPartials, clock, label));
        }

        String cladeKey(NodeRef node) {
            List<String> taxa = new ArrayList<>();
            collectTaxa(node, taxa);
            Collections.sort(taxa);
            return String.join(",", taxa);
        }

        private void collectTaxa(NodeRef node, List<String> taxa) {
            if (tree.isExternal(node)) {
                taxa.add(tree.getNodeTaxon(node).getId());
            } else {
                for (int c = 0; c < tree.getChildCount(node); c++) {
                    collectTaxa(tree.getChild(node, c), taxa);
                }
            }
        }
    }

    /** K=2, three heterochronous tips, strongly asymmetric migration, unequal N. */
    static Spec c1() {
        return new Spec("c1", 2, "((A:0.8,B:0.5):1.2,C:1.6);")
                .tips("A=0", "B=1", "C=1")
                .migration(new double[]{0, 0.05, 1.5, 0})
                .popSizes(new double[]{0.3, 3.0});
    }

    /** K=3, four tips (two tied), very unequal N. */
    static Spec c2() {
        return new Spec("c2", 3, "((A:0.5,B:0.5):0.7,(C:0.4,D:0.9):0.6);")
                .tips("A=0", "B=2", "C=1", "D=1")
                .migration(new double[]{0, 0.3, 0.1, 0.2, 0, 0.6, 1.0, 0.05, 0})
                .popSizes(new double[]{0.1, 1.0, 10.0});
    }

    /** K=3, five tips, one ambiguous between states 0 and 2. */
    static Spec c3() {
        return new Spec("c3", 3, "(((A:0.4,B:0.6):0.5,C:1.0):0.6,(D:0.7,E:1.3):0.9);")
                .ambiguity("X", "0", "2")
                .tips("A=0", "B=1", "C=X", "D=2", "E=1")
                .migration(new double[]{0, 0.5, 0.2, 0.3, 0, 0.4, 0.1, 0.6, 0})
                .popSizes(new double[]{0.8, 1.5, 0.5});
    }

    /** K=2, four tips, two epochs of migration rates and population sizes. */
    static Spec c4() {
        return new Spec("c4", 2, "((A:0.6,B:0.55):0.9,(C:0.5,D:1.2):1.1);")
                .tips("A=0", "B=1", "C=1", "D=0")
                .epochs(1.3)
                .migration(new double[]{0, 0.4, 0.8, 0}, new double[]{0, 1.2, 0.3, 0})
                .popSizes(new double[]{0.5, 2.0}, new double[]{2.0, 0.5});
    }

    /** K=3, four tips, strict clock rate 2.5 (scales migration only). */
    static Spec c5() {
        return new Spec("c5", 3, "(((A:0.3,B:0.3):0.4,C:0.9):0.5,D:1.2);")
                .tips("A=0", "B=1", "C=2", "D=0")
                .migration(new double[]{0, 0.2, 0.05, 0.1, 0, 0.3, 0.25, 0.15, 0})
                .popSizes(new double[]{1.0, 0.5, 2.0})
                .strictClock(2.5);
    }

    /** K=2, five tips, two epochs and per-branch rates. */
    static Spec c6() {
        return new Spec("c6", 2, "(((A:0.3,B:0.5):0.4,C:0.6):0.8,(D:0.9,E:1.1):0.6);")
                .tips("A=0", "B=1", "C=1", "D=0", "E=1")
                .epochs(1.0)
                .migration(new double[]{0, 0.6, 0.9, 0}, new double[]{0, 0.2, 1.4, 0})
                .popSizes(new double[]{0.7, 1.2}, new double[]{1.5, 0.4})
                .branchRate("A", 0.5).branchRate("B", 1.5).branchRate("A,B", 0.8)
                .branchRate("C", 1.2).branchRate("D", 0.6).branchRate("E", 1.4).branchRate("D,E", 0.9)
                .branchRate("A,B,C", 1.1);
    }

    /**
     * K=3, only the branch above C migrates (every other branch has rate 0 and a known state), across two epochs.
     * At any time at most one lineage is mobile and the others are in known states, so Pi_h, the continuous
     * frozen-environment surrogate and the exact structured coalescent coincide.
     */
    static Spec singleMobile() {
        return new Spec("singleMobile", 3, "((A:1.0,B:1.0):0.8,C:1.5);")
                .ambiguity("X", "0", "2")
                .tips("A=1", "B=1", "C=X")
                .epochs(1.2)
                .migration(new double[]{0, 0.7, 0.3, 0.4, 0, 0.9, 0.5, 0.2, 0},
                        new double[]{0, 0.2, 1.1, 0.6, 0, 0.3, 0.8, 0.4, 0})
                .popSizes(new double[]{0.6, 0.9, 1.3}, new double[]{1.1, 0.4, 0.7})
                .branchRate("A", 0.0).branchRate("B", 0.0).branchRate("A,B", 0.0).branchRate("C", 1.3);
    }

    /**
     * The state of tests/TestXML/testMascotHighPressureLogging.xml, whose complexSubstitutionModel lists the rates
     * as (0,1), (0,2), (1,2), (1,0), (2,0), (2,1).
     */
    static Spec highPressure() {
        return new Spec("highPressure", 3,
                "(((A:2.8877280032784474,B:2.8877280032784474):1.5358558603633852," +
                        "(C:3.3883756235523608,D:3.3883756235523608):1.0352082400894718):1.4805952691315536," +
                        "(E:0.5493739204766552,F:0.5493739204766552):5.354805212296731);")
                .tips("A=0", "B=1", "C=2", "D=0", "E=1", "F=2")
                .migration(new double[]{0, 0.022840096758276696, 0.9991617072287613,
                        0.006329302680772018, 0, 0.6039584418543463,
                        0.6532680506509443, 0.004136765472939712, 0})
                .popSizes(new double[]{0.15729183324413404, 0.20647504338270894, 0.03191953596439705});
    }

    static List<Spec> oracleCases() {
        List<Spec> cases = new ArrayList<>();
        cases.add(c1());
        cases.add(c2());
        cases.add(c3());
        cases.add(c4());
        cases.add(c5());
        cases.add(c6());
        return cases;
    }
}
