package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.CompoundSymmetricMatrix;
import dr.inference.model.Model;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class SampledMultivariateTraitLikelihood extends AbstractMultivariateTraitLikelihood {

    public SampledMultivariateTraitLikelihood(String traitName,
                                              TreeModel treeModel,
                                              MultivariateDiffusionModel diffusionModel,
                                              CompoundParameter traitParameter,
                                              List<Integer> missingIndices,
                                              boolean cacheBranches, boolean scaleByTime, boolean useTreeLength,
                                              BranchRateModel rateModel, Model samplingDensity,
                                              boolean reportAsMultivariate,
                                              boolean reciprocalRates) {
        super(traitName, treeModel, diffusionModel, traitParameter, missingIndices, cacheBranches, scaleByTime,
                useTreeLength, rateModel, samplingDensity, reportAsMultivariate, reciprocalRates);
    }

    protected String extraInfo() {
        return "\tSampling internal trait values: true\n";
    }

    /**
     * Calculate the log likelihood of the current state.
     *
     * @return the log likelihood.
     */
    public double calculateLogLikelihood() {

        double logLikelihood;

        if (!cacheBranches)
            logLikelihood = traitLogLikelihood(null, treeModel.getRoot());
        else
            logLikelihood = traitCachedLogLikelihood(null, treeModel.getRoot());
        if (logLikelihood > maxLogLikelihood) {
            maxLogLikelihood = logLikelihood;
        }
        return logLikelihood;
    }

    protected  double calculateAscertainmentCorrection(int taxonIndex) {
        throw new RuntimeException("Ascertainment correction not yet implemented for sampled trait likelihoods");
    }

    public final double getLogDataLikelihood() {
        double logLikelihood = 0;
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef tip = treeModel.getExternalNode(i); // TODO Do not include integrated tips; how to check???

            if (cacheBranches && validLogLikelihoods[tip.getNumber()])
                logLikelihood += cachedLogLikelihoods[tip.getNumber()];
            else {
                NodeRef parent = treeModel.getParent(tip);

                double[] tipTrait = treeModel.getMultivariateNodeTrait(tip, traitName);
                double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, traitName);
                double time = getRescaledBranchLength(tip);

                logLikelihood += diffusionModel.getLogLikelihood(parentTrait, tipTrait, time);
            }
        }
        return logLikelihood;
    }


    private double traitCachedLogLikelihood(double[] parentTrait, NodeRef node) {

        double logL = 0.0;
        double[] childTrait = null;
        final int nodeNumber = node.getNumber();

        if (!treeModel.isRoot(node)) {

            if (!validLogLikelihoods[nodeNumber]) { // recompute

                childTrait = treeModel.getMultivariateNodeTrait(node, traitName);
                double time = getRescaledBranchLength(node);
                if (parentTrait == null)
                    parentTrait = treeModel.getMultivariateNodeTrait(treeModel.getParent(node), traitName);
                logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
                cachedLogLikelihoods[nodeNumber] = logL;
                validLogLikelihoods[nodeNumber] = true;
            } else
                logL = cachedLogLikelihoods[nodeNumber];
        }

        int childCount = treeModel.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            logL += traitCachedLogLikelihood(childTrait, treeModel.getChild(node, i));
        }

        return logL;
    }

    private double traitLogLikelihood(double[] parentTrait, NodeRef node) {

        double logL = 0.0;
        double[] childTrait = treeModel.getMultivariateNodeTrait(node, traitName);

        if (parentTrait != null) {

            double time = getRescaledBranchLength(node);
            logL = diffusionModel.getLogLikelihood(parentTrait, childTrait, time);
            if (new Double(logL).isNaN()) {
                System.err.println("AbstractMultivariateTraitLikelihood: likelihood is undefined");
                System.err.println("time = " + time);
                System.err.println("parent trait value = " + new Vector(parentTrait));
                System.err.println("child trait value = " + new Vector(childTrait));

                double[][] precisionMatrix = diffusionModel.getPrecisionmatrix();
                if (precisionMatrix != null) {
                    System.err.println("precision matrix = " + new Matrix(diffusionModel.getPrecisionmatrix()));
                    if (diffusionModel.getPrecisionParameter() instanceof CompoundSymmetricMatrix) {
                        CompoundSymmetricMatrix csMatrix = (CompoundSymmetricMatrix) diffusionModel.getPrecisionParameter();
                        System.err.println("diagonals = " + new Vector(csMatrix.getDiagonals()));
                        System.err.println("off diagonal = " + csMatrix.getOffDiagonal());
                    }
                }
            }
        }
        int childCount = treeModel.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            logL += traitLogLikelihood(childTrait, treeModel.getChild(node, i));
        }

        if (new Double(logL).isNaN()) {
            System.err.println("logL = " + logL);
//            System.err.println(new Matrix(diffusionModel.getPrecisionmatrix()));
            System.exit(-1);
        }

        return logL;
    }


    public double[] getTraitForNode(Tree treeModel, NodeRef node, String traitName) {
        return ((TreeModel)treeModel).getMultivariateNodeTrait(node,traitName);
    }

}