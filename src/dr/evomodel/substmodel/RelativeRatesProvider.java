/*
 * RelativeRatesProvider.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.substmodel.nucleotide.TN93;
import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.List;

/**
 * Rate provider for a reversible model whose relative rates r_ij (in BaseSubstitutionModel pair order) are
 * constant on the MRCA classes of the character tree: HKY, TN93 and GTR with tied rates. With the dense
 * normalisation norm = 2 sum_{i<j} r_ij pi_i pi_j, the clade exchangeability is s(m) = r_ij / norm for any pair
 * with MRCA m, and lambda_root = s(root), lambda_m = pi_m (s(m) - s(parent(m))).
 * <p>
 * The generator is valid iff s(m) >= 0 at every internal node (docs/generator-constraints). Since norm > 0,
 * that is r_ij >= 0, so the bounds put on the parameters are the exact constraints: kappa >= 0 for HKY,
 * kappa1, kappa2 >= 0 for TN93, every rate >= 0 for GTR. kappa < 1 gives negative clade rates inside the valid
 * cone; only the hierarchical BEAGLE plugin, which needs non-negative clade rates, excludes them.
 *
 * @author Xiang Ji
 */
public class RelativeRatesProvider extends AbstractHierarchicalRateProvider {

    public interface RelativeRates {
        void setup(double[] rates);
    }

    public static final double DEFAULT_TOLERANCE = 1e-10;

    private final RelativeRates relativeRates;
    private final boolean normalise;
    private final double tolerance;
    private final double[] rates;
    private final double[] sigma;
    private final int[] representativePair;
    private final int[][] classPairs;

    public RelativeRatesProvider(String name, TreeModel tree, FrequencyModel frequency, RelativeRates relativeRates,
                                 boolean normalise, double tolerance, Variable... variables) {
        super(name, tree, frequency);
        this.relativeRates = relativeRates;
        this.normalise = normalise;
        this.tolerance = tolerance;
        this.rates = new double[stateCount * (stateCount - 1) / 2];
        this.sigma = new double[tree.getNodeCount()];
        this.representativePair = new int[tree.getNodeCount()];
        this.classPairs = new int[tree.getNodeCount()][];
        for (int n : internalPreOrder) {
            final NodeRef node = tree.getNode(n);
            final List<Integer> pairs = new ArrayList<Integer>();
            for (int a = 0; a < tree.getChildCount(node); a++) {
                final int[] left = statesBelow[tree.getChild(node, a).getNumber()];
                for (int b = a + 1; b < tree.getChildCount(node); b++) {
                    final int[] right = statesBelow[tree.getChild(node, b).getNumber()];
                    for (int x : left) {
                        for (int y : right) {
                            pairs.add(pairIndex(x, y, stateCount));
                        }
                    }
                }
            }
            classPairs[n] = toArray(pairs);
            representativePair[n] = classPairs[n][0];
        }
        for (Variable variable : variables) {
            addVariable(variable);
        }
    }

    public static int pairIndex(int i, int j, int stateCount) {
        if (i > j) {
            final int t = i;
            i = j;
            j = t;
        }
        return (2 * stateCount - i - 1) * i / 2 + (j - i - 1);
    }

    public static RelativeRatesProvider hky(Parameter kappa, FrequencyModel frequency) {
        return hky(kappa, frequency, characterTree(NUCLEOTIDE_NEWICK));
    }

    public static RelativeRatesProvider hky(final Parameter kappa, FrequencyModel frequency, TreeModel tree) {
        checkNucleotides(frequency, "kappa, kappa1 and kappa2");
        kappa.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        final RelativeRatesProvider provider = new RelativeRatesProvider("hkyRates", tree, frequency, new RelativeRates() {
            public void setup(double[] rates) {
                final double k = kappa.getParameterValue(0);
                rates[0] = 1.0;
                rates[1] = k;
                rates[2] = 1.0;
                rates[3] = 1.0;
                rates[4] = k;
                rates[5] = 1.0;
            }
        }, true, DEFAULT_TOLERANCE, kappa);
        provider.checkPurinePyrimidineClades();
        return provider;
    }

    public static RelativeRatesProvider tn93(Parameter kappa1, Parameter kappa2, FrequencyModel frequency) {
        return tn93(kappa1, kappa2, frequency, characterTree(NUCLEOTIDE_NEWICK));
    }

    public static RelativeRatesProvider tn93(final Parameter kappa1, final Parameter kappa2, FrequencyModel frequency,
                                             TreeModel tree) {
        checkNucleotides(frequency, "kappa, kappa1 and kappa2");
        kappa1.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        kappa2.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        final RelativeRatesProvider provider = new RelativeRatesProvider("tn93Rates", tree, frequency, new RelativeRates() {
            public void setup(double[] rates) {
                rates[0] = 1.0;
                rates[1] = kappa1.getParameterValue(0);
                rates[2] = 1.0;
                rates[3] = 1.0;
                rates[4] = kappa2.getParameterValue(0);
                rates[5] = 1.0;
            }
        }, true, DEFAULT_TOLERANCE, kappa1, kappa2);
        provider.checkPurinePyrimidineClades();
        return provider;
    }

    public static RelativeRatesProvider gtr(final Parameter sixRates, FrequencyModel frequency, TreeModel tree,
                                            double tolerance) {
        checkNucleotides(frequency, "gtr rates");
        if (sixRates.getDimension() != 6) {
            throw new IllegalArgumentException("gtr rates need 6 values (AC, AG, AT, CG, CT, GT); found "
                    + sixRates.getDimension());
        }
        sixRates.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 6));
        final RelativeRatesProvider provider = new RelativeRatesProvider("gtrRates", tree, frequency, new RelativeRates() {
            public void setup(double[] rates) {
                for (int k = 0; k < 6; k++) {
                    rates[k] = sixRates.getParameterValue(k);
                }
            }
        }, true, tolerance, sixRates);
        provider.getNormalization();
        return provider;
    }

    /**
     * The provider for an existing substitution model that is a special case of the hierarchical model: HKY and
     * TN93 for now (other models are added with their own providers). The kappas are read through the model's
     * getters and the provider listens to the model, so every change of the model reaches the hierarchical
     * rates; the provider shares the model's FrequencyModel.
     */
    public static RelativeRatesProvider of(final SubstitutionModel dense, TreeModel tree) {
        final FrequencyModel frequency = dense.getFrequencyModel();
        final RelativeRates relativeRates;
        final String name;
        if (dense instanceof HKY) {
            final HKY hky = (HKY) dense;
            name = "hkyRates";
            relativeRates = new RelativeRates() {
                public void setup(double[] rates) {
                    final double k = hky.getKappa();
                    rates[0] = 1.0;
                    rates[1] = k;
                    rates[2] = 1.0;
                    rates[3] = 1.0;
                    rates[4] = k;
                    rates[5] = 1.0;
                }
            };
        } else if (dense instanceof TN93) {
            final TN93 tn93 = (TN93) dense;
            name = "tn93Rates";
            relativeRates = new RelativeRates() {
                public void setup(double[] rates) {
                    rates[0] = 1.0;
                    rates[1] = tn93.getKappa1();
                    rates[2] = 1.0;
                    rates[3] = 1.0;
                    rates[4] = tn93.getKappa2();
                    rates[5] = 1.0;
                }
            };
        } else {
            throw new IllegalArgumentException("ratesFrom supports HKYModel and TN93Model; found "
                    + dense.getModelName());
        }
        checkNucleotides(frequency, "ratesFrom models");
        final RelativeRatesProvider provider = new RelativeRatesProvider(name, tree, frequency, relativeRates, true,
                DEFAULT_TOLERANCE);
        provider.checkPurinePyrimidineClades();
        provider.addModel(dense);
        return provider;
    }

    private static void checkNucleotides(FrequencyModel frequency, String what) {
        if (!(frequency.getDataType() instanceof Nucleotides)) {
            throw new IllegalArgumentException(what + " need a nucleotide frequency model; found "
                    + frequency.getDataType().getDescription());
        }
    }

    private void checkPurinePyrimidineClades() {
        final List<String> found = new ArrayList<String>();
        boolean ag = false, ct = false;
        for (int n : internalPreOrder) {
            if (parentOf[n] >= 0) {
                found.add(cladeName[n]);
                ag |= cladeName[n].equals("AG");
                ct |= cladeName[n].equals("CT");
            }
        }
        if (!ag || !ct || found.size() != 2) {
            throw new IllegalArgumentException("the character tree must have the clades {A,G} and {C,T}; found "
                    + describeClades());
        }
    }

    public static double normalisation(double[] rates, double[] pi) {
        final int stateCount = pi.length;
        double norm = 0.0;
        int k = 0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = i + 1; j < stateCount; j++) {
                norm += rates[k++] * pi[i] * pi[j];
            }
        }
        return 2.0 * norm;
    }

    @Override
    protected double computeLambdas(double[] lambda) {
        relativeRates.setup(rates);
        final double norm = normalise ? normalisation(rates, frequency.getFrequencies()) : 1.0;
        for (int n : internalPreOrder) {
            final double s = rates[representativePair[n]] / norm;
            for (int k : classPairs[n]) {
                final double other = rates[k] / norm;
                if (Math.abs(other - s) > tolerance * Math.max(1.0, Math.abs(s))) {
                    throw new IllegalArgumentException("exchangeabilities are not tree-structured on the character tree at clade "
                            + cladeName[n] + ": " + other + " vs " + s);
                }
            }
            sigma[n] = s;
            lambda[n] = parentOf[n] < 0 ? s : piNode[n] * (s - sigma[parentOf[n]]);
        }
        return norm;
    }

    public double getSigma(String sortedTipCodes) {
        getNormalization();
        return sigma[nodeOfClade(sortedTipCodes)];
    }
}
