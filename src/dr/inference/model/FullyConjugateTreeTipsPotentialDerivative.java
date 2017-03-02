package dr.inference.model;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;

/**
 * @author Max Tolkoff
 */
public class FullyConjugateTreeTipsPotentialDerivative implements PotentialDerivativeInterface {
    FullyConjugateMultivariateTraitLikelihood treeLikelihood;

    public FullyConjugateTreeTipsPotentialDerivative(FullyConjugateMultivariateTraitLikelihood treeLikelihood){
        this.treeLikelihood = treeLikelihood;
    }

    @Override
    public double[] getDerivative() {
        Parameter traitParameter = treeLikelihood.getTraitParameter();
        int dimTraits = treeLikelihood.getDimTrait() * treeLikelihood.getNumData();
        int ntaxa = traitParameter.getDimension() / dimTraits;
        double[] derivative = new double[traitParameter.getDimension()];
        double[][] mean = treeLikelihood.getConditionalMeans();
        double[] precfactor = treeLikelihood.getPrecisionFactors(); //TODO need to be flexible to handle diagonal case and dense case. Right now just diagonal.


        for (int i = 0; i < dimTraits; i++) {
            for (int j = 0; j < ntaxa; j++) {
                derivative[j * dimTraits + i] += (traitParameter.getParameterValue(j * dimTraits + i) - mean[j][i]) * precfactor[j];
            }

        }


        return derivative;
    }
}
