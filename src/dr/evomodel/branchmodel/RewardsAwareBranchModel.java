/*
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

package dr.evomodel.branchmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.tree.TreeModel;
import dr.inference.markovjumps.SericolaSeriesMarkovRewardFastModel;
import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * TransitionMatrixProviderBranchModel that uses SericolaSeriesMarkovRewardFastModelRho
 * to compute reward-conditioned transition matrices on each branch.
 *
 * Design goals:
 *  - BranchModel only:
 *      (i) packs branch rewards X and times,
 *      (ii) provides output buffers W[nodeNr],
 *      (iii) caches "knownTransitionMatrices" at the branch-model level.
 *
 * @author Filippo Monti
 */
public class RewardsAwareBranchModel extends AbstractModel
        implements TransitionMatrixProviderBranchModel, Citable, Reportable {

    public static final String REWARDS_AWARE_BRANCH_MODEL = "RewardsAwareBranchModel";

    private final Parameter rewardRates;
    private final SubstitutionModel underlyingSubstitutionModel;
    private final TreeModel tree;
    private final ArbitraryBranchRates branchRateModel;
    private final Parameter indicator;                  // 0/1, same indexing as rho
    private final IndexedParameter indexedParameter;    // same indexing as rho

    private final int nstates;
    private final int dim2;

    // Output transition matrices in ORIGINAL state order, per node index.
    // Root row is unused but kept for direct nodeNr indexing.
    private final double[][] W;
    private final double[][] Watomic;

    // Packed inputs for Sericola: one entry per non-root branch (nodeCount - 1)
    private final double[] X;          // total reward
    private final double[] times;      // branchLength
    private final double[][] Wpacked;  // references into W[nodeNr], one per branch

    private final SericolaSeriesMarkovRewardFastModel sericola;

    // Cache flag at this branch-model layer
    private boolean knownTransitionMatrices = false;

    // For compatibility with SubstitutionModelDelegate
    private List<SubstitutionModel> substitutionModels;

    // Optional debug/testing
    private boolean DUMMYTESTING = false;
    private boolean DEBUG = false;

    private final int[] branchIndexToNodeNr;
    private final int[] nodeNrToBranchIndex;

    public RewardsAwareBranchModel(TreeModel tree,
                                   SubstitutionModel underlyingSubstitutionModel,
                                   Parameter rewardRates,
                                   Parameter indicator,
                                   IndexedParameter indexedParameter,
                                   ArbitraryBranchRates branchRateModel,
                                   Parameter extremeIndices,
                                   boolean conditional) {

        super(REWARDS_AWARE_BRANCH_MODEL);
        if (tree == null) throw new IllegalArgumentException("tree must be non-null");
        if (underlyingSubstitutionModel == null) {
            throw new IllegalArgumentException("RewardsAwareBranchModel must be provided with an underlying substitution model");
        }
        if (rewardRates == null) throw new IllegalArgumentException("rewardRates must be non-null");
        if (branchRateModel == null) throw new IllegalArgumentException("branchRateModel must be non-null");

        this.tree = tree;
        this.underlyingSubstitutionModel = underlyingSubstitutionModel;
        this.rewardRates = rewardRates;
        this.branchRateModel = branchRateModel;
        this.indicator = indicator;
        this.indexedParameter = indexedParameter;

        final int dim = branchRateModel.getRateParameter().getDimension();
        if (indicator.getDimension() != dim) {
            throw new IllegalArgumentException("indicator dimension must equal rho dimension (branchRateModel rate parameter).");
        }
        if (indexedParameter.getDimension() != dim) {
            throw new IllegalArgumentException("indexedParameter dimension must equal rho dimension (branchRateModel rate parameter).");
        }


        this.nstates = underlyingSubstitutionModel.getDataType().getStateCount();
        this.dim2 = nstates * nstates;

        final int nodeCount = tree.getNodeCount();
        final int branchCount = nodeCount - 1; // all non-root nodes

        this.W = new double[nodeCount][dim2];
        this.Watomic = new double[nodeCount][dim2];

        this.X = new double[branchCount];
        this.times = new double[branchCount];
        this.Wpacked = new double[branchCount][];

        final double epsilon = 1e-10;
        this.sericola = new SericolaSeriesMarkovRewardFastModel( //only for cts values
                underlyingSubstitutionModel,
                rewardRates,
                nstates,
                epsilon,
                conditional
        );

        addModel(tree);
        addModel(branchRateModel);
        addModel(sericola);
        addVariable(indicator);
        addVariable(indexedParameter.getIndicesParameter());
        addVariable(indexedParameter.getValuesParameter());
        addVariable(extremeIndices);

        final int nNodes = tree.getNodeCount();
        final int nBranches = nNodes - 1;

        branchIndexToNodeNr = new int[nBranches];
        nodeNrToBranchIndex = new int[nNodes];
        Arrays.fill(nodeNrToBranchIndex, -1);

        int k = 0;
        for (int i = 0; i < nNodes; i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            int nodeNr = node.getNumber();
            branchIndexToNodeNr[k] = nodeNr;
            nodeNrToBranchIndex[nodeNr] = k;
            k++;
        }
    }
    public int getNodeNumberForBranchIndex(int branchIndex) {
        return branchIndexToNodeNr[branchIndex];
    }
    private int getParameterIndexForNode(final int nodeNr) {
        final NodeRef node = tree.getNode(nodeNr);
        // This method name is from your commented code; if your version differs, swap accordingly.
        return branchRateModel.getParameterIndexFromNode(node);
    }

    public int getBranchIndexForNodeNumber(int nodeNr) {
        return nodeNrToBranchIndex[nodeNr];
    }
    public Parameter getIndicator() { return indicator; }
    public IndexedParameter getIndexedParameter() { return indexedParameter; }


    // -------------------- Basic accessors --------------------

    public FrequencyModel getRootFrequencyModel() { return underlyingSubstitutionModel.getFrequencyModel(); }

    @Override
    public SubstitutionModel getRootSubstitutionModel() { return underlyingSubstitutionModel; }

    public TreeModel getTree() { return tree; }

    public Parameter getRewardRates() { return rewardRates; }

    public BranchRateModel getRateBranchModel() { return branchRateModel; }

    public SericolaSeriesMarkovRewardFastModel getSericolaModel() { return sericola; }

    public double getUniformizationRate() {
        return sericola.getUniformizationRate();
    }

    // -------------------- Main API --------------------

    @Override
    public double[] getTransitionMatrix(NodeRef branch) {
        return getTransitionMatrix(branch.getNumber());
    }

    public double[] getTransitionMatrix(int nodeNr) {
        if (DUMMYTESTING) {
            NodeRef node = tree.getNode(nodeNr);
            double t = tree.getBranchLength(node);
            getRootSubstitutionModel().getTransitionProbabilities(t, W[nodeNr]);
            return W[nodeNr];
        }
        final int p = getParameterIndexForNode(nodeNr);
        final boolean atomicOn = isOne(indicator.getParameterValue(p));
        if (atomicOn) {
            return getTransitionMatrixAtomic(nodeNr);
        } else {
            return getTransitionMatrixCts(nodeNr);
        }
    }

    public double[] getTransitionMatrixCts(int nodeNr) {
        computeTransitionMatrices();
        if (DEBUG) {
            final double[] w = W[nodeNr];

            final double tol = 0.0; // or something like -1e-12 if you want numerical tolerance
            for (int i = 0; i < w.length; i++) {
                if (!(w[i] > tol)) {
                    throw new IllegalStateException(
                            "Transition matrix for node " + nodeNr +
                                    " contains non-positive entry at index " + i +
                                    ": value=" + w[i]
                    );
                }
            }
        }
        return W[nodeNr];
    }

    public double[] getTransitionMatrixAtomic(int nodeNr) {
        NodeRef node = tree.getNode(nodeNr);
        final double t = tree.getBranchLength(node);
        final double lambda = sericola.getLambda();
        final double scale = Math.exp(-lambda * t);
        final int paramIndex = getParameterIndexForNode(nodeNr);
        final int atomState = indexedParameter.getIndexValue(paramIndex);
        Arrays.fill(Watomic[nodeNr], 0.0);
        Watomic[nodeNr][atomState * nstates + atomState] = scale;
        return Watomic[nodeNr];
    }

    private static boolean isOne(final double x) {
        final long r = Math.round(x);
        return (Math.abs(x - r) <= 1e-9) && (r == 1L);
    }

    private void computeTransitionMatrices() {

        if (knownTransitionMatrices) return;
        if (DEBUG) System.out.println("RewardsAwareBranchModel: computeTransitionMatrices");

        final int nodeCount = tree.getNodeCount();


        if (DEBUG) {
//            int k = 0;
//            for (int i = 0; i < nodeCount; i++) {
//                final NodeRef node = tree.getNode(i);
//                if (tree.isRoot(node)) continue;
//
//                final int nodeNr = node.getNumber();
//                final double t = tree.getBranchLength(node);
//                final double rate = branchRateModel.getBranchRate(tree, node);
//
//                times[k] = t;
//                X[nodeNr] = t;
//
//                // Write results directly into W[nodeNr] (original order)
//                Wpacked[k] = W[nodeNr];
//                k++;
//            }
//                    System.err.print("X = ");
//        for (int i = 0; i < k; i++) {
//            System.err.print(X[i]);
//            if (i < k - 1) System.err.print(", ");
//        }
//            throw new IllegalStateException("DUMMY EXIT");
        }

//        double[] newValues = {2.17510446, 1.0, 1.1976619470000003, 3.025158366999996, 3.687673462999996, 4.0, 1.0, 4.0, 2.738096880999997, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 4.116476411000001, 1.0, 1.0, 3.6905935650000004, 5.0, 4.156313334000004, 5.0, 5.0, 3.197307416000001, 5.0, 4.938555213999997, 3.354191017000005, 5.0, 5.0, 4.467388521000004, 5.0, 4.3368921249999985, 1.0, 1.0, 5.0, 4.192920567999998, 1.9472911279999963, 4.517477322999994, 2.0602452609999986, 5.0, 4.0, 4.0, 3.3583967830000034, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.2885495749999976, 4.478665653, 5.0, 5.0, 2.481351459999999, 1.0, 1.0, 1.1848177740000025, 2.639734472999997, 4.0, 4.0, 4.950464521999997, 2.302737268999998, 4.0, 1.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 4.092127269999999, 4.072177705999998, 5.0, 4.0, 4.0, 1.498872034999998, 1.0, 1.0, 4.0, 4.0, 4.026512656000001, 5.0, 3.2332541090000007, 1.0, 1.0, 4.732089666, 4.940175623000002, 4.0, 4.0, 4.0, 2.857107292000002, 4.303362914999994, 4.0, 2.888708698000002, 4.0, 2.528267219, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.9276579819999995, 4.0, 1.0, 3.978308990000002, 4.0, 1.9957667160000057, 4.0, 4.0, 4.0, 2.247304284000002, 4.0, 4.979677209999998, 1.0, 1.0, 4.7234432190000035, 4.370942845999998, 5.0, 2.6411178140000047, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 1.0, 1.9122807570000049, 3.909500688999998, 3.6012112459999983, 2.285114275000005, 4.495522447999996, 5.0, 5.0, 4.9230427259999985, 3.9943087509999984, 5.0, 5.0, 5.0, 5.0, 2.564895012000008, 4.156215953, 1.000000000000007, 1.000000000000007, 2.811844958000002, 1.0, 1.0, 4.788114549999996, 4.0, 4.0, 4.0, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.408012491000001, 4.0, 5.0, 5.0, 1.0, 1.0, 1.0, 2.5069966730000033, 4.765514087999996, 1.6547086830000026, 4.0, 4.348457635999999, 1.9572366200000033, 2.1842416269999987, 4.297788584999999, 3.2892549869999996, 4.122712053999997, 3.496593279999999, 4.0, 4.0, 4.0, 4.8692309249999965, 4.417875105999997, 3.3970025469999996, 3.7588856229999976, 2.493398499999998, 1.0, 1.0, 4.834267717000003, 1.0, 1.1814947590000031, 1.7369990689999995, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 1.4129041440000023, 5.000000000000007, 4.278261620000002, 5.0, 5.0, 4.790492442999998, 1.0, 1.0, 1.0, 3.944705501999998, 1.632951718000001, 4.0, 1.7180738369999986, 5.0, 1.0, 1.0, 1.0, 1.0, 4.0, 4.0, 2.4813434949999973, 1.0, 4.1377083870000035, 3.914694589, 5.0, 5.0, 5.0, 5.0, 4.0, 4.0, 1.0, 5.0, 5.0, 1.0, 3.525796061000001, 1.0, 1.0, 1.0, 1.0, 1.5438018500000013, 3.1049545660000035, 1.0, 4.0, 4.0, 4.0, 4.699216808999999, 1.0, 2.5450417309999978, 5.0, 5.0, 2.0780470249999965, 1.0, 1.0, 5.0, 2.307770472999998, 5.0, 5.0, 4.012763737, 4.497399731999998, 5.0, 5.0, 4.6189770060000015, 1.0, 1.0, 1.0555149070000027, 2.9608575589999973, 4.0, 4.0, 4.372870311, 1.0, 3.1442411800000016, 4.0, 2.4983181599999966, 4.005384270999997, 4.0, 3.335627926000001, 2.933060306999998, 1.0, 3.6231112580000016, 1.0, 3.7861011930000004, 3.535606989999998, 1.0, 3.327503526000001, 3.5998157599999985, 4.0, 2.3172913340000036, 1.0, 1.0, 3.5024392289999966, 1.0, 1.0, 2.067716828000002, 4.0, 5.0, 4.0, 4.0, 4.0, 4.224263470000004, 5.0, 5.0, 2.741279720999998, 5.0, 4.0, 4.0, 2.6578477479999947, 2.835420597999999, 3.463229650999992, 1.0, 2.442566453000005, 5.0, 5.0, 2.6726882509999967, 5.0, 5.0, 3.7485276620000008, 1.0, 2.109939008000005, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.6881604800000076, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 1.314759256000002, 1.0, 1.0, 1.0, 3.2371001590000077, 4.0, 4.0, 1.0, 1.0, 3.6640951640000026, 1.0, 3.2239830580000017, 1.0, 2.802222971000006, 1.1902254360000057, 1.0, 1.5134443640000086, 1.0, 1.0, 1.0, 1.0, 1.0, 2.813791124000005, 1.994877945000006, 2.6254002720000074, 1.0, 1.0, 1.0, 3.506649685999996, 4.274445258, 1.0, 4.0, 4.0, 4.0, 3.1785986249999993, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.990632235999996, 4.625932239999997, 5.0, 5.0, 5.0, 4.051218147, 2.9290768449999973, 1.0, 1.0, 1.0, 1.839633206000002, 4.946140276000008, 4.143956343999996, 4.0, 2.1254908880000016, 4.0, 5.0, 5.0, 5.0, 5.0, 4.007189147999995, 1.5606179460000078, 3.4238191079999964, 5.0, 2.2812037079999996, 5.0, 4.544288444000003, 1.901372718999994, 1.0, 1.0, 4.0, 3.9092583719999965, 1.0, 1.0, 2.3426929230000013, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.424814769000001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0, 4.605194885000003, 2.2215265450000032, 1.0, 1.0, 4.478348295999993, 5.0, 4.831013135999996, 4.308030147000004, 4.0, 3.2481775009999936, 5.0, 3.542813226, 4.0, 4.442159719999999, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.7133067850000003, 1.0, 1.0, 2.241186855999999, 4.462588572000001, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.410405710999996, 4.0, 3.2949944930000044, 5.0, 4.0, 4.0, 4.4402345310000015, 5.0, 5.0, 5.0, 5.0, 5.0, 2.355008386999998, 3.9921064970000018, 2.7157906810000014, 5.0, 5.0, 5.0, 4.511831967000006, 5.0, 5.0, 5.0, 4.388779009000004, 4.0, 2.617407078999996, 5.0, 4.602153080000001, 4.287858096999997, 2.2102736800000002, 1.0, 1.0, 1.000000000000007, 4.0, 4.0, 1.3641576500000099, 1.0, 2.075727256999997, 4.293322576000001, 3.3599211109999914, 3.6383145339999956, 3.6078501779999996, 2.7220940839999983, 4.172063551999997, 4.083835688999997, 1.0, 2.603969448000001, 5.0, 3.192112991000002, 3.275256845999998, 4.0, 4.291442461999999, 4.0, 4.558713505, 4.1573153430000005, 4.2850716129999995, 4.0, 4.871713298000003, 2.767315916000001, 5.0, 3.3295863019999956, 4.0, 4.648197197999998, 3.672656691999997, 1.0, 1.0, 3.2630166550000013, 3.7390196409999987, 5.0, 1.4859382790000026, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 2.2134231830000033, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0265610579999986, 4.0, 4.746200268000003, 5.0, 1.7276993500000017, 1.0, 3.597978947999998, 1.0, 1.0, 2.243924053999997, 1.0, 1.0, 1.0, 2.3879684260000005, 4.101883588, 5.0, 2.306895951999998, 1.0, 1.0, 3.164731334999999, 1.0, 1.0, 1.0, 2.8739454939999973, 1.0, 1.0, 1.0, 1.0, 1.0, 1.060576367000003, 1.0, 2.497577810000003, 1.0, 5.0, 3.453973433999998, 1.0, 1.0, 1.0, 1.0, 1.0, 4.045361385, 1.0, 3.507001068000001, 1.0, 1.0, 1.0, 5.0, 4.269122152999998, 4.0, 4.0, 1.0, 2.0452696210000028, 4.0, 4.0, 4.0, 4.0, 3.259220880000001, 1.0, 3.378752193000004, 4.880000241999994, 3.8233784179999972, 5.0, 5.0, 4.822342354, 4.0, 4.0, 4.0, 3.1403317380000004, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.7966224529999977, 3.2311437059999975, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.606617235999998, 4.0, 4.0, 4.885441958999998, 5.0, 5.0, 5.0, 4.506123582000001, 4.0, 4.406390926, 4.322611705, 3.1891804400000012, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 2.2646077899999995, 5.0, 4.018150379000005, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5759997810000073, 1.0, 1.0, 1.0, 3.8938309220000065, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.979660527000007, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.036749490999995, 3.5901265939999973, 3.4428472750000054};
//        double[] newValues = {2.17510446, 1.0, 1.1976619470000003, 3.025158366999996, 3.687673462999996, 4.0, 1.0, 4.0, 2.738096880999997, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 4.116476411000001, 1.0, 1.0, 3.6905935650000004, 5.0, 4.156313334000004, 5.0, 5.0, 3.197307416000001, 5.0, 4.938555213999997, 3.354191017000005, 5.0, 5.0, 4.467388521000004, 5.0, 4.3368921249999985, 1.0, 1.0, 5.0, 4.192920567999998, 1.9472911279999963, 4.517477322999994, 2.0602452609999986, 5.0, 4.0, 4.0, 3.3583967830000034, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 3.2885495749999976, 4.478665653, 5.0, 5.0, 2.481351459999999, 1.0, 1.0, 1.1848177740000025, 2.639734472999997, 4.0, 4.0, 4.950464521999997, 2.302737268999998, 4.0, 1.0, 1.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 1.0, 1.0, 1.0, 4.092127269999999, 4.072177705999998, 5.0, 4.0, 4.0, 1.498872034999998, 1.0, 1.0, 4.0, 4.0, 4.026512656000001, 5.0, 3.2332541090000007, 1.0, 1.0, 4.732089666, 4.940175623000002, 4.0, 4.0, 4.0, 2.857107292000002, 4.303362914999994, 4.0, 2.888708698000002, 4.0, 2.528267219, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.9276579819999995, 4.0, 1.0, 3.978308990000002, 4.0, 1.9957667160000057, 4.0, 4.0, 4.0, 2.247304284000002, 4.0, 4.979677209999998, 1.0, 1.0, 4.7234432190000035, 4.370942845999998, 5.0, 2.6411178140000047, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 5.0, 5.0, 5.0, 1.0, 1.9122807570000049, 3.909500688999998, 3.6012112459999983, 2.285114275000005, 4.495522447999996, 5.0, 5.0, 4.9230427259999985, 3.9943087509999984, 5.0, 5.0, 5.0, 5.0, 2.564895012000008, 4.156215953, 1.000000000000007, 1.000000000000007, 2.811844958000002, 1.0, 1.0, 4.788114549999996, 4.0, 4.0, 4.0, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.408012491000001, 4.0, 5.0, 5.0, 1.0, 1.0, 1.0, 2.5069966730000033, 4.765514087999996, 1.6547086830000026, 4.0, 4.348457635999999, 1.9572366200000033, 2.1842416269999987, 4.297788584999999, 3.2892549869999996, 4.122712053999997, 3.496593279999999, 4.0, 4.0, 4.0, 4.8692309249999965, 4.417875105999997, 3.3970025469999996, 3.7588856229999976, 2.493398499999998, 1.0, 1.0, 4.834267717000003, 1.0, 1.1814947590000031, 1.7369990689999995, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 1.4129041440000023, 5.000000000000007, 4.278261620000002, 5.0, 5.0, 4.790492442999998, 1.0, 1.0, 1.0, 3.944705501999998, 1.632951718000001, 4.0, 1.7180738369999986, 5.0, 1.0, 1.0, 1.0, 1.0, 4.0, 4.0, 2.4813434949999973, 1.0, 4.1377083870000035, 3.914694589, 5.0, 5.0, 5.0, 5.0, 4.0, 4.0, 1.0, 5.0, 5.0, 1.0, 3.525796061000001, 1.0, 1.0, 1.0, 1.0, 1.5438018500000013, 3.1049545660000035, 1.0, 4.0, 4.0, 4.0, 4.699216808999999, 1.0, 2.5450417309999978, 5.0, 5.0, 2.0780470249999965, 1.0, 1.0, 5.0, 2.307770472999998, 5.0, 5.0, 4.012763737, 4.497399731999998, 5.0, 5.0, 4.6189770060000015, 1.0, 1.0, 1.0555149070000027, 2.9608575589999973, 4.0, 4.0, 4.372870311, 1.0, 3.1442411800000016, 4.0, 2.4983181599999966, 4.005384270999997, 4.0, 3.335627926000001, 2.933060306999998, 1.0, 3.6231112580000016, 1.0, 3.7861011930000004, 3.535606989999998, 1.0, 3.327503526000001, 3.5998157599999985, 4.0, 2.3172913340000036, 1.0, 1.0, 3.5024392289999966, 1.0, 1.0, 2.067716828000002, 4.0, 5.0, 4.0, 4.0, 4.0, 4.224263470000004, 5.0, 5.0, 2.741279720999998, 5.0, 4.0, 4.0, 2.6578477479999947, 2.835420597999999, 3.463229650999992, 1.0, 2.442566453000005, 5.0, 5.0, 2.6726882509999967, 5.0, 5.0, 3.7485276620000008, 1.0, 2.109939008000005, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.6881604800000076, 4.0, 4.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 1.314759256000002, 1.0, 1.0, 1.0, 3.2371001590000077, 4.0, 4.0, 1.0, 1.0, 3.6640951640000026, 1.0, 3.2239830580000017, 1.0, 2.802222971000006, 1.1902254360000057, 1.0, 1.5134443640000086, 1.0, 1.0, 1.0, 1.0, 1.0, 2.813791124000005, 1.994877945000006, 2.6254002720000074, 1.0, 1.0, 1.0, 3.506649685999996, 4.274445258, 1.0, 4.0, 4.0, 4.0, 3.1785986249999993, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.990632235999996, 4.625932239999997, 5.0, 5.0, 5.0, 4.051218147, 2.9290768449999973, 1.0, 1.0, 1.0, 1.839633206000002, 4.946140276000008, 4.143956343999996, 4.0, 2.1254908880000016, 4.0, 5.0, 5.0, 5.0, 5.0, 4.007189147999995, 1.5606179460000078, 3.4238191079999964, 5.0, 2.2812037079999996, 5.0, 4.544288444000003, 1.901372718999994, 1.0, 1.0, 4.0, 3.9092583719999965, 1.0, 1.0, 2.3426929230000013, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.424814769000001, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.0, 4.605194885000003, 2.2215265450000032, 1.0, 1.0, 4.478348295999993, 5.0, 4.831013135999996, 4.308030147000004, 4.0, 3.2481775009999936, 5.0, 3.542813226, 4.0, 4.442159719999999, 5.0, 5.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0, 3.7133067850000003, 1.0, 1.0, 2.241186855999999, 4.462588572000001, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 4.410405710999996, 4.0, 3.2949944930000044, 5.0, 4.0, 4.0, 4.4402345310000015, 5.0, 5.0, 5.0, 5.0, 5.0, 2.355008386999998, 3.9921064970000018, 2.7157906810000014, 5.0, 5.0, 5.0, 4.511831967000006, 5.0, 5.0, 5.0, 4.388779009000004, 4.0, 2.617407078999996, 5.0, 4.602153080000001, 4.287858096999997, 2.2102736800000002, 1.0, 1.0, 1.000000000000007, 4.0, 4.0, 1.3641576500000099, 1.0, 2.075727256999997, 4.293322576000001, 3.3599211109999914, 3.6383145339999956, 3.6078501779999996, 2.7220940839999983, 4.172063551999997, 4.083835688999997, 1.0, 2.603969448000001, 5.0, 3.192112991000002, 3.275256845999998, 4.0, 4.291442461999999, 4.0, 4.558713505, 4.1573153430000005, 4.2850716129999995, 4.0, 4.871713298000003, 2.767315916000001, 5.0, 3.3295863019999956, 4.0, 4.648197197999998, 3.672656691999997, 1.0, 1.0, 3.2630166550000013, 3.7390196409999987, 5.0, 1.4859382790000026, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 2.2134231830000033, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0265610579999986, 4.0, 4.746200268000003, 5.0, 1.7276993500000017, 1.0, 3.597978947999998, 1.0, 1.0, 2.243924053999997, 1.0, 1.0, 1.0, 2.3879684260000005, 4.101883588, 5.0, 2.306895951999998, 1.0, 1.0, 3.164731334999999, 1.0, 1.0, 1.0, 2.8739454939999973, 1.0, 1.0, 1.0, 1.0, 1.0, 1.060576367000003, 1.0, 2.497577810000003, 1.0, 5.0, 3.453973433999998, 1.0, 1.0, 1.0, 1.0, 1.0, 4.045361385, 1.0, 3.507001068000001, 1.0, 1.0, 1.0, 5.0, 4.269122152999998, 4.0, 4.0, 1.0, 2.0452696210000028, 4.0, 4.0, 4.0, 4.0, 3.259220880000001, 1.0, 3.378752193000004, 4.880000241999994, 3.8233784179999972, 5.0, 5.0, 4.822342354, 4.0, 4.0, 4.0, 3.1403317380000004, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.7966224529999977, 3.2311437059999975, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 4.606617235999998, 4.0, 4.0, 4.885441958999998, 5.0, 5.0, 5.0, 4.506123582000001, 4.0, 4.406390926, 4.322611705, 3.1891804400000012, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 5.0, 5.0, 5.0, 2.2646077899999995, 5.0, 4.018150379000005, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.5759997810000073, 1.0, 1.0, 1.0, 3.8938309220000065, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.979660527000007, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 3.036749490999995, 3.5901265939999973, 3.4428472750000054};
        //        System.arraycopy(newValues, 0, X, 0, X.length);
        double minRate = Double.POSITIVE_INFINITY;
        double maxRate = Double.NEGATIVE_INFINITY;
        double sumRate = 0.0;
        int count = 0;
        int k = 0;
        for (int i = 0; i < nodeCount; i++) {
            final NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            final int nodeNr = node.getNumber();
            final double t = tree.getBranchLength(node);
            final double rate = branchRateModel.getBranchRate(tree, node);

            times[k] = t;
            X[k] = rate;

            // Write results directly into W[nodeNr] (original order)
            Wpacked[k] = W[nodeNr];
            k++;
        }

        if (DEBUG) {
            System.err.print("times = ");
            for (int i = 0; i < k; i++) {
                System.err.print(times[i] - 1);
                if (i < k - 1) System.err.print(" ");
            }
            System.err.println("]");
            throw new RuntimeException("finished printing");
        }

        // Sericola handles all sorting/lazy caches and writes into ORIGINAL order by design
//        sericola.computePdfIntoY(X, times, true, Wpacked);
        sericola.computePdfInto(X, times, true, Wpacked);
        knownTransitionMatrices = true;
    }
    public double[] getWPacked(int i) {
        return Wpacked[i];
    }

    // -------------------- Branch model mapping --------------------

    @Override
    public Mapping getBranchModelMapping(NodeRef node) {
        final double[] weights = new double[]{1.0};
        final int[] order = new int[]{node.getNumber()};

        return new Mapping() {
            @Override
            public int[] getOrder() { return order; }

            @Override
            public double[] getWeights() { return weights; }
        };
    }

    @Override
    public boolean requiresMatrixConvolution() {
        return false;
    }

    // -------------------- Model events --------------------

    private boolean ignoreModelChangedEvent = false;

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (ignoreModelChangedEvent) return;
        knownTransitionMatrices = false;
        fireModelChanged();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        knownTransitionMatrices = false;
        fireModelChanged();
    }

    // -------------------- MCMC store/restore --------------------

    private boolean storedKnownTransitionMatrices;

    @Override
    protected void storeState() {
        storedKnownTransitionMatrices = knownTransitionMatrices;
    }

    @Override
    protected void restoreState() {
        knownTransitionMatrices = false;
    }

    @Override
    protected void acceptState() {
        // nothing
//        System.out.println("Accepting state inside RewardsAwareBranchModel");
    }

    public int getStateCount() {
        return nstates;
    }

    // -------------------- Citable / Reportable --------------------

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Rewards Aware Branch model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(
                new Citation(new Author[]{new Author("F", "Monti"),
                        new Author("MA", "Suchard")},
                        "Dependencies between CTMCs",
                        2026,
                        "TOBE",
                        1,
                        1,
                        1,
                        Citation.Status.IN_PRESS));
    }

    @Override
    public String getReport() {
        if (!DUMMYTESTING) computeTransitionMatrices();
        StringBuilder sb = new StringBuilder();
        sb.append("W matrix: ");
        for (int nodeNr = 0; nodeNr < tree.getNodeCount(); nodeNr++) {
            NodeRef node = tree.getNode(nodeNr);
            if (tree.isRoot(node)) continue;
            for (double val : W[nodeNr]) sb.append(val).append(" ");
        }
        return sb.toString();
    }

    // -------------------- Compatibility: SubstitutionModelDelegate --------------------

    @Override
    public List<SubstitutionModel> getSubstitutionModels() {
        if (substitutionModels == null) {
            buildSubstitutionModels();
        }
        return substitutionModels;
    }

    protected void buildSubstitutionModels() {

        substitutionModels = new ArrayList<>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            if (tree.isRoot(node)) continue;

            ignoreModelChangedEvent = true;
            SubstitutionModel substitutionModel = new TransitionMatrixProvider(
                    "RewardsAwareSubstitutionModel",
                    underlyingSubstitutionModel.getDataType(),
                    underlyingSubstitutionModel.getFrequencyModel(),
                    W[node.getNumber()]
            );
            ignoreModelChangedEvent = false;

            substitutionModels.add(substitutionModel);
        }
    }

    class TransitionMatrixProvider extends ComplexSubstitutionModel {
        private final double[] transitionMatrixRef;

        public TransitionMatrixProvider(String name, DataType dataType,
                                        FrequencyModel freqModel, double[] transitionMatrixRef) {
            super(name, dataType, freqModel, null);
            this.transitionMatrixRef = transitionMatrixRef;
        }

        @Override
        public void getTransitionProbabilities(double distance, double[] matrix) {
            // Ensure current before exposing.
            if (!DUMMYTESTING) computeTransitionMatrices();
            System.arraycopy(transitionMatrixRef, 0, matrix, 0, transitionMatrixRef.length);
            throw new UnsupportedOperationException("This provider has not yet implemented the atomic case");
        }

        @Override protected void handleModelChangedEvent(Model model, Object object, int index) {}
        @Override protected void frequenciesChanged() {}
        @Override protected void ratesChanged() {}
        @Override protected void setupRelativeRates(double[] rates) {}
    }
}
