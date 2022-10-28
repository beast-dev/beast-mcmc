package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treelikelihood.AncestralStateBeagleTreeLikelihood;
import dr.inference.distribution.BinomialLikelihood;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Statistic;
import dr.math.Binomial;
import dr.math.MathUtils;
import dr.math.UnivariateFunction;
import dr.math.UnivariateMinimum;
import dr.math.distributions.GammaDistribution;
import dr.xml.Reportable;

import java.util.Set;

/**
 * @author Andy Magee
 */
public class ASRBranchSubstitutionMatrix extends Statistic.Abstract implements Reportable {

    public ASRBranchSubstitutionMatrix(String name,
                                       AncestralStateBeagleTreeLikelihood asrLike,
                                       DataType dataType,
                                       Boolean doublets,
                                       TaxonList mrcaTaxaDescendant

    ) throws TreeUtils.MissingTaxonException {
        this.name = name;
        this.asrLikelihood = asrLike;
        this.tree = asrLikelihood.getTreeModel();
        this.dataType = dataType;
        if (dataType == null) { throw new RuntimeException("DataType not found;");}
        this.leafSetDescendant = (mrcaTaxaDescendant != null) ? TreeUtils.getLeavesForTaxa(tree, mrcaTaxaDescendant) : null;
        // Check validity of paired data model

        this.useDoublets = doublets;
        int d = dataType.getStateCount() * dataType.getStateCount();
        if ( useDoublets ) {
            d *= d;
        }
        this.dim = d;
        this.dimNames = constructDimNames();
    }

    public int getDimension() {
        return dim;
    }
    
    public String getDimensionName(int i) {
        return dimNames[i];
    }

    private int[] singlesFromDoublet(int i, int nStates) {
        int[] singlets = new int[2];
        singlets[0] = Math.floorDiv(i, nStates);
        singlets[1] = i - (singlets[0] * nStates);
        return singlets;
    }

    public String[] constructDimNames() {
        int n = dataType.getStateCount();
        String[] dimNames = new String[dim];
        if ( useDoublets ) {
            int nSq = n * n;
            for (int i = 0; i < dim; i++) {
                int fromDoublet = Math.floorDiv(i, nSq);
                int toDoublet = i - fromDoublet * nSq;
                System.err.println("from " + fromDoublet + " to " + toDoublet);
                int[] from = singlesFromDoublet(fromDoublet, n);
                int[] to = singlesFromDoublet(toDoublet, n);
                System.err.println("  from[0] " + from[0] + " from[1] " + from[1]);
                System.err.println("  to[0] " + to[0] + " to[1] " + to[1]);
//                dimNames[i] = name + "[" + String.valueOf(dataType.getChar(from[0])) + String.valueOf(dataType.getChar(from[1])) + "," + String.valueOf(dataType.getChar(to[0])) + String.valueOf(dataType.getChar(to[1])) + "]";
                dimNames[i] = name + "[" + dataType.getChar(from[0]) + dataType.getChar(from[1]) + "," + dataType.getChar(to[0]) + dataType.getChar(to[1]) + "]";
            }
        } else {
            for (int i = 0; i < dim; i++) {
                int from = Math.floorDiv(i, n);
                int to = i - from * n;
                dimNames[i] = name + "[" + dataType.getChar(from) + "," + dataType.getChar(to) + "]";
            }
        }
        return dimNames;
    }

    public String getStatisticName() {
        return NAME;
    }

    private NodeRef getNode(Set<String> leafSet) {
        return (leafSet != null) ? TreeUtils.getCommonAncestorNode(tree, leafSet) : tree.getRoot();
    }

    private double[] getSingletCountMatrix(int[] nodeStatesAncestor, int[] nodeStatesDescendant) {
        int n = dataType.getStateCount();
        double[] countMatrix = new double[n * n];
        for (int s = 0; s < nodeStatesAncestor.length - 1; s++) {
            int from =  nodeStatesAncestor[s];
            int to =  nodeStatesDescendant[s];
            countMatrix[n * from + to] += 1.0;
        }
        return countMatrix;
    }


    private double[] getDoubletCountMatrix(int[] nodeStatesAncestor, int[] nodeStatesDescendant) {
        int n = dataType.getStateCount();
        int nSq = n * n;
        double[] doubletCountMatrix = new double[nSq * nSq];
        for (int s = 0; s < nodeStatesAncestor.length - 1; s++) {
            int from =  n * nodeStatesAncestor[s] + nodeStatesAncestor[s + 1];
            int to =  n * nodeStatesDescendant[s] + nodeStatesDescendant[s + 1];
            doubletCountMatrix[nSq * from + to] += 1.0;
        }
        return doubletCountMatrix;
    }

    private int getDoublet(int first, int second, int nStates) {
        return first * nStates + second;
    }

    private int[] getTargetsForDim(int dim) {
        int[] targets = new int[2];
        int n = dataType.getStateCount();

        if ( useDoublets ) {
            int nSq = n * n;
            targets[0] = Math.floorDiv(dim, nSq);
            targets[1] = dim - targets[0] * nSq;
        } else {
            targets[0] = Math.floorDiv(dim, n);
            targets[1] = dim - targets[0] * n;
        }

        return targets;
    }

    public double getStatisticValue(int dim) {
        // TODO doing this separately nStates^2 times is a bit wasteful
        int n = dataType.getStateCount();
        int nSq = n * n;

        double count = 0.0;

        int targets[] = getTargetsForDim(dim);

        NodeRef nodeDescendant = getNode(leafSetDescendant);
        NodeRef nodeAncestor = tree.getParent(nodeDescendant);

        int[] ancestorStates = asrLikelihood.getStatesForNode(tree, nodeAncestor);
        int[] descendantStates = asrLikelihood.getStatesForNode(tree, nodeDescendant);

        int last = ancestorStates.length;
        if ( useDoublets ) {
            last--;
        }

        for (int s = 0; s < last; s++) {
            int from;
            int to;
            if (useDoublets) {
                from = getDoublet(ancestorStates[s], ancestorStates[s + 1], n);
                to = getDoublet(descendantStates[s], descendantStates[s + 1], n);
            } else {
                from = ancestorStates[s];
                to = descendantStates[s];
            }
            if (from == targets[0] && to == targets[1]) {
                count += 1.0;
            }
        }
        return count;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder("ASRBranchSubstitutionMatrix Report\n\n");
        for (int i = 0; i < dim; i++) {
            sb.append(dimNames[i] + ": " + getStatisticValue(i) + "; ");
        }
        sb.append(getStatisticValue(0));
        return sb.toString();
    }

    private AncestralStateBeagleTreeLikelihood asrLikelihood;
    private final Set<String> leafSetDescendant;
    private final Tree tree;
    private final String name;
    private final boolean useDoublets;
    private final int dim;
    private final DataType dataType;
    private final String[] dimNames;
}
