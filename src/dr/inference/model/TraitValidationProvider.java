package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inferencexml.model.TraitValidationProviderParser;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * @author Gabriel Hassler
 */

public class TraitValidationProvider implements CrossValidationProvider, Reportable {

    private final Parameter trueTraits;
    private final int[] missingInds;
    private final String[] dimNames;
    private final String sumName;
    private final int dimTrait;
    private final ContinuousExtensionDelegate extensionDelegate;


    public TraitValidationProvider(Parameter trueTraits,
                                   ContinuousTraitPartialsProvider dataModel,
                                   Tree treeModel,
                                   String id,
                                   Parameter missingParameter,
                                   TreeDataLikelihood treeLikelihood,
                                   String inferredValuesName,
                                   List<Integer> trueMissingIndices) {


        this.trueTraits = trueTraits;
        this.dimTrait = dataModel.getPartialDimension();

        this.missingInds = setupMissingInds(dataModel, missingParameter, trueMissingIndices);
        int nMissing = missingInds.length;


        TreeTrait treeTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + inferredValuesName);

        if (dataModel instanceof ModelExtensionProvider) {
            this.extensionDelegate = ((ModelExtensionProvider) dataModel).getExtensionDelegate(
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate(),
                    treeTrait,
                    treeModel);
        } else { //Simply returns the tree traits
            this.extensionDelegate = new ContinuousExtensionDelegate(
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate(),
                    treeTrait,
                    treeModel);
        }

        this.dimNames = new String[nMissing];
        setupDimNames(treeModel, id);
        sumName = getSumName(id);

    }

    private int[] setupMissingInds(ContinuousTraitPartialsProvider dataModel, Parameter missingParameter,
                                   List<Integer> trueMissing) {
        int[] missingInds;
        int nMissing = 0;
        if (missingParameter == null) {
            List<Integer> missingList = dataModel.getMissingIndices();
            missingList.removeAll(trueMissing);
            nMissing = missingList.size();

            missingInds = new int[nMissing];

            for (int i = 0; i < nMissing; i++) {
                missingInds[i] = missingList.get(i);
            }

        } else {


            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0 && !trueMissing.contains(i)) {
                    //TODO: search more efficiently through the `trueMissing` array
                    nMissing += 1;
                }
            }

            missingInds = new int[nMissing];
            int counter = 0;

            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0 && !trueMissing.contains(i)) {
                    //TODO: (see above)
                    missingInds[counter] = i;
                    counter += 1;
                }
            }

        }

        return missingInds;

    }

    private String getId(String id) {
        if (id == null) return TraitValidationProviderParser.TRAIT_VALIDATION_PROVIDER;
        return id;
    }

    private void setupDimNames(Tree treeModel, String id) {
        id = getId(id);
        int dim = 0;
        for (int i : missingInds) {
            int taxonInd = i / dimTrait;
            int traitInd = i - taxonInd * dimTrait;
            String taxonName = treeModel.getTaxonId(taxonInd);
            dimNames[dim] = id + "." + taxonName + "." + (traitInd + 1);
            dim += 1;
        }
    }

    private String getSumName(String id) {
        id = getId(id);
        return id + ".TotalSum";
    }

    @Override
    public double[] getTrueValues() {
        return trueTraits.getParameterValues();
    }

    @Override
    public double[] getInferredValues() {
        return extensionDelegate.getExtendedValues();
    }

//    private void updateTraitsFromTree() {
//        double[] tipValues = (double[]) treeTrait.getTrait(treeLikehoood.getTree(), null);
//        assert (tipValues.length == inferredTraits.getDimension());
//        for (int i : missingInds) {
//            inferredTraits.setParameterValueQuietly(i, tipValues[i]);
//        }
//        inferredTraits.fireParameterChangedEvent();
//    }

    @Override
    public int[] getRelevantDimensions() {
        return missingInds;
    }

    @Override
    public String getName(int dim) {
        return dimNames[dim];
    }

    @Override
    public String getNameSum(int dim) {
        return sumName;
    }

    @Override
    public String getReport() {
        int iterations = 1000000;
        int nMissing = missingInds.length;

        double[] sums = new double[nMissing];
        double[][] sumSquares = new double[nMissing][nMissing];

        for (int iteration = 0; iteration < iterations; iteration++) {

            double[] extendedVals = extensionDelegate.getExtendedValues();

            for (int i = 0; i < nMissing; i++) {

                double vi = extendedVals[missingInds[i]];

                sums[i] += vi;
                sumSquares[i][i] += vi * vi;

                for (int j = 0; j < i; j++) {

                    double vj = extendedVals[missingInds[j]];
                    sumSquares[i][j] += vi * vj;
                }
            }
        }

        for (int i = 0; i < nMissing; i++) {
            sums[i] = sums[i] / iterations;
            sumSquares[i][i] = sumSquares[i][i] / iterations - sums[i] * sums[i];
            for (int j = 0; j < i; j++) {
                sumSquares[i][j] = sumSquares[i][j] / iterations - sums[i] * sums[j];
                sumSquares[j][i] = sumSquares[i][j];
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Mean:\t");
        sb.append(new Vector(sums));
        sb.append("\n");
        sb.append("Covariance:\t");
        sb.append(new Matrix(sumSquares));


        return sb.toString();
    }
}
