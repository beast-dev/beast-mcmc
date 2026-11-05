package dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.treedatalikelihood.AugmentedNodeRegistry;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EpochEvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.xml.Reportable;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public final class DiscretePreOrderReport implements Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final PreOrderMessageProvider discreteDelegate;
    private final TreeTrait treeTrait;
    private final double roundingTolerance;
    private final int roundingScale;
    private final DiscretePartialsType displayType;
    private final DiscretePartialsType beagleType;
    private final EpochEvolutionaryProcessDelegate epochProcessDelegate; // null unless there are augmented epoch nodes

    public DiscretePreOrderReport(TreeDataLikelihood treeDataLikelihood) {
        this(treeDataLikelihood, DiscretePartialsType.TOP,0.0);
    }

    public DiscretePreOrderReport(TreeDataLikelihood treeDataLikelihood,
                                  double roundingTolerance) {
        this(treeDataLikelihood, DiscretePartialsType.TOP, roundingTolerance);
    }

    public DiscretePreOrderReport(TreeDataLikelihood treeDataLikelihood,
                                  DiscretePartialsType type,
                                  double roundingTolerance) {
        this.treeDataLikelihood = treeDataLikelihood;
        this.roundingTolerance = roundingTolerance;
        this.roundingScale = roundingTolerance > 0.0 ?
                Math.max(0, BigDecimal.valueOf(roundingTolerance).stripTrailingZeros().scale()) : -1;
        this.displayType = type;
        if (type == DiscretePartialsType.TOP_SPECTRAL) {
            this.beagleType = DiscretePartialsType.TOP;
        } else if (type == DiscretePartialsType.BOTTOM_SPECTRAL) {
            this.beagleType = DiscretePartialsType.BOTTOM;
        } else {
            this.beagleType = displayType;
        }

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (delegate instanceof PreOrderMessageProvider) {
            this.discreteDelegate = (PreOrderMessageProvider) delegate;
            this.treeTrait = null;
            this.epochProcessDelegate = null;
        } else if (delegate instanceof BeagleDataLikelihoodDelegate) {

            BeagleDataLikelihoodDelegate likelihoodDelegate = (BeagleDataLikelihoodDelegate) delegate;
            this.epochProcessDelegate = (likelihoodDelegate.getEvolutionaryProcessDelegate()
                    instanceof EpochEvolutionaryProcessDelegate) ?
                    (EpochEvolutionaryProcessDelegate) likelihoodDelegate.getEvolutionaryProcessDelegate() : null;
            String name = "preorder.partials";

            if (treeDataLikelihood.getTreeTrait(name) == null) {

                ProcessSimulationDelegate preOrderDelegate = new AbstractBeagleGradientDelegate(
                        "test", treeDataLikelihood.getTree(), likelihoodDelegate) {

                    @Override
                    protected DiscretePartialsType getPreOrderType() {
                        return beagleType;
                    }

                    @Override
                    protected int getGradientLength() {
                        return 0;
                    }

                    @Override
                    protected void getNodeDerivatives(Tree tree, double[] first, double[] second) {

                    }

                    @Override
                    public void getPreorderPartials(int nodeNumber, DiscretePartialsType type, double[] out) {

                        assert out.length == stateCount * categoryCount * patternCount;

                        if (type == DiscretePartialsType.TOP || displayType == DiscretePartialsType.TOP_SPECTRAL) {

                            if (displayType == DiscretePartialsType.TOP_SPECTRAL) {
                                rotatePreorderPartials(nodeNumber, DiscretePartialsType.TOP, false, out);
                            } else {
                                super.getPreorderPartials(nodeNumber, DiscretePartialsType.TOP, out);
                            }

                        } else if (type == DiscretePartialsType.BOTTOM || displayType == DiscretePartialsType.BOTTOM_SPECTRAL) {

                            if (displayType == DiscretePartialsType.BOTTOM_SPECTRAL) {
                                rotatePreorderPartials(nodeNumber, DiscretePartialsType.BOTTOM, true, out);
                            } else {
                                super.getPreorderPartials(nodeNumber, DiscretePartialsType.BOTTOM, out);
                            }
                        } else {
                            super.getPreorderPartials(nodeNumber, type, out);
                        }
                    }

                    /**
                     * Rotates the raw (state-space) pre-order partial of nodeNumber -- an original or, with
                     * augmented epoch nodes, a degree-2 node -- into the eigenbasis of the substitution model that
                     * applies on the branch immediately above it (the model that makes that branch's transition
                     * matrix diagonal in this basis; see SpectralBeagleCrossProductDelegate.getRotatedPartial).
                     * Root has no branch above it, so its partial (the root state frequencies) is left unrotated.
                     * <p>
                     * A branch spanning more than one epoch (i.e., requiring convolution of several substitution
                     * models, rather than the augmented degree-2 nodes that give each segment its own single model)
                     * has no single diagonalizing basis; getSubstitutionModelForBranch() throws in that case, which
                     * is the correct signal here too.
                     */
                    private void rotatePreorderPartials(int nodeNumber, DiscretePartialsType rawType, boolean bottom,
                                                        double[] out) {

                        double[] intermediate = new double[out.length];
                        super.getPreorderPartials(nodeNumber, rawType, intermediate);

                        if (nodeNumber == tree.getRoot().getNumber()) {
                            System.arraycopy(intermediate, 0, out, 0, out.length);
                            return;
                        }

                        EigenDecomposition ed = likelihoodDelegate.getEvolutionaryProcessDelegate()
                                .getSubstitutionModelForBranch(nodeNumber).getEigenDecomposition();
                        EigenDecomposition edT = ed.transpose();
                        double[] rotation = bottom ? edT.getEigenVectors() : edT.getInverseEigenVectors();

                        int offset = 0;
                        for (int c = 0; c < categoryCount; ++c) {
                            for (int p = 0; p < patternCount; ++p) {
                                multiplyMatrixVector(rotation, intermediate, offset, out, offset, stateCount);
                                offset += stateCount;
                            }
                        }
                    }

                    @Override
                    protected void constructTraits(Helper treeTraitHelper) {

                        treeTraitHelper.addTrait(new TreeTrait.DA() {
                            @Override
                            public String getTraitName() {
                                return name;
                            }

                            @Override
                            public Intent getIntent() {
                                return Intent.WHOLE_TREE;
                            }

                            @Override
                            public double[] getTrait(Tree tree, NodeRef node) {
                                return getGradient(node);
                            }
                        });
                    }
                };

                TreeTraitProvider traitProvider = new ProcessSimulation(treeDataLikelihood, preOrderDelegate);
                treeDataLikelihood.addTraits(traitProvider.getTreeTraits());

                this.discreteDelegate = (PreOrderMessageProvider) preOrderDelegate;
                this.treeTrait = treeDataLikelihood.getTreeTrait("preorder.partials");

            } else {
                throw new RuntimeException("Unknown error");
            }
        } else {
            throw new IllegalArgumentException("Unknown likelihood");
        }
    }

    @Override
    public String getReport() {

        treeDataLikelihood.getLogLikelihood();
        Tree tree = treeDataLikelihood.getTree();

        if (treeTrait != null) {
            treeTrait.getTrait(tree, null);
        }

        int categoryCount = discreteDelegate.getCategoryCount();
        int patternCount = discreteDelegate.getPatternCount();
        int stateCount = discreteDelegate.getStateCount();

        double[] allStart = new double[stateCount * categoryCount * patternCount];
        double[] allEnd = new double[stateCount * categoryCount * patternCount];

        double[] start = new double[stateCount];
        double[] end = new double[stateCount];

        StringBuilder sb = new StringBuilder();
        sb.append("Discrete pre-order partials\n");
        sb.append("  treeDataLikelihood = ").append(treeDataLikelihood.getId()).append('\n');
        sb.append("  nodes = ").append(tree.getNodeCount()).append('\n');
        sb.append("  patterns = ").append(patternCount).append('\n');
        sb.append("  categories = ").append(categoryCount).append('\n');
        sb.append("  states = ").append(stateCount).append('\n');

        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            int nodeNumber = node.getNumber();
            sb.append("node ").append(nodeNumber);
            if (tree.isRoot(node)) {
                sb.append(" root");
            }
            if (tree.isExternal(node)) {
                Taxon taxon = tree.getNodeTaxon(node);
                sb.append(" taxon=").append(taxon == null ? "null" : taxon.getId());
            }
            sb.append('\n');

            appendPartialRows(sb, nodeNumber, allStart, allEnd, start, end, categoryCount, patternCount, stateCount);
        }

        // With augmented epoch nodes, each branch also carries a chain of degree-2 nodes, one per epoch
        // transition time it spans, each with its own pre-order partial.
        if (epochProcessDelegate != null) {
            AugmentedNodeRegistry registry = epochProcessDelegate.getAugmentedNodeRegistry();
            int[] chain = new int[registry.getBoundaryCount()];

            for (int branch = 0; branch < tree.getNodeCount(); branch++) {
                if (tree.isRoot(tree.getNode(branch))) {
                    continue;
                }

                int count = registry.copyChain(branch, chain);
                for (int j = 0; j < count; j++) {
                    int id = chain[j];

                    sb.append("node ").append(id).append(" augmented, branch ").append(branch)
                            .append(", epoch ").append(registry.getMatrixEpoch(id)).append('\n');

                    appendPartialRows(sb, id, allStart, allEnd, start, end, categoryCount, patternCount, stateCount);
                }
            }
        }

        return sb.toString();
    }

    private void appendPartialRows(StringBuilder sb, int nodeNumber, double[] allStart, double[] allEnd,
                                   double[] start, double[] end, int categoryCount, int patternCount,
                                   int stateCount) {

        discreteDelegate.getPreorderPartials(nodeNumber, DiscretePartialsType.TOP, allStart);
        discreteDelegate.getPreorderPartials(nodeNumber, DiscretePartialsType.BOTTOM, allEnd);

        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {

                int offset = (c * patternCount + p) * stateCount;
                System.arraycopy(allStart, offset, start, 0, stateCount);
                System.arraycopy(allEnd, offset, end, 0, stateCount);

                sb.append("  category ").append(c)
                        .append(" pattern ").append(p)
                        .append(" start=").append(formatArray(start))
                        .append(" end=").append(formatArray(end))
                        .append('\n');
            }
        }
    }

    private static void multiplyMatrixVector(double[] matrix, double[] vector, int vectorOffset, double[] out, int outOffset, int dim) {
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            final int rowBase = i * dim;
            for (int j = 0; j < dim; j++) {
                sum += matrix[rowBase + j] * vector[vectorOffset + j];
            }
            out[outOffset + i] = sum;
        }
    }

    private String formatArray(double[] values) {
        if (roundingTolerance <= 0.0) {
            return Arrays.toString(values);
        }

        StringBuilder sb = new StringBuilder("[");
        for (int i = 0; i < values.length; i++) {
            if (i > 0) {
                sb.append(", ");
            }
            BigDecimal rounded = BigDecimal.valueOf(values[i]).setScale(roundingScale, RoundingMode.HALF_UP);
            sb.append(rounded.stripTrailingZeros().toPlainString());
        }
        sb.append(']');
        return sb.toString();
    }
}
