package dr.evomodel.operators;

import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.AbstractHamiltonianMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;
import dr.math.MathUtils;

/**
 * Created by max on 12/2/15.
 */
public class LatentFactorHamiltonianMC extends AbstractHamiltonianMCOperator{
    private LatentFactorModel lfm;
    private FullyConjugateMultivariateTraitLikelihood tree;
    private MatrixParameter factors;
    private MatrixParameter loadings;
    private MatrixParameter Precision;
    private int nfac;
    private int ntaxa;
    private int ntraits;
    private double stepSize;
    private int nSteps;


    public LatentFactorHamiltonianMC(LatentFactorModel lfm, FullyConjugateMultivariateTraitLikelihood tree, double weight, CoercionMode mode, double stepSize, int nSteps, double momentumSd){
        super(mode, momentumSd);
        setWeight(weight);
        this.lfm=lfm;
        this.tree=tree;
        this.factors=lfm.getFactors();
        this.loadings=lfm.getLoadings();
        this.Precision=lfm.getColumnPrecision();
        nfac=lfm.getFactorDimension();
        ntaxa=lfm.getFactors().getColumnDimension();
        ntraits=Precision.getRowDimension();
        this.stepSize=stepSize;
        this.nSteps=nSteps;
    }



    @Override
    public double getCoercableParameter() {
        return 0;
    }

    @Override
    public void setCoercableParameter(double value) {

    }

    @Override
    public double getRawParameter() {
        return 0;
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public String getOperatorName() {
        return "Latent Factor Hamiltonian Monte Carlo";
    }

    private double[] getMatrix(int nfac, double[] residual){
        double answer[]=new double[this.nfac];
        for (int i = 0; i <this.nfac ; i++) {
            for (int j = 0; j < ntraits; j++) {
                answer[i] +=loadings.getParameterValue(i,j)*Precision.getParameterValue(j,j)*
                        residual[j*ntaxa+nfac];
            }
        }
        return answer;
    }

    private double[] getGradient(int randel, double[] mean, double[][] prec){
        double[] residual=lfm.getResidual();
        double[] derivative=getMatrix(randel, residual);


        for (int i = 0; i <mean.length ; i++) {
            double sumi=0;
            for (int j = 0; j <mean.length ; j++) {
                sumi+=prec[i][j]*(factors.getParameterValue(j, randel)-mean[j]);
            }
            derivative[i]-=sumi;
        }
        return derivative;
    }

    @Override
    public double doOperation() throws OperatorFailedException {
        int randel = MathUtils.nextInt(ntaxa);



        double[] mean=tree.getConditionalMean(randel);
        double[][] prec=tree.getConditionalVariance(randel);

        double[] derivative=getGradient(randel, mean, prec);
        drawMomentum(lfm.getFactorDimension());

        double prop=0;
        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            prop+=momentum[i]*momentum[i]/2;
        }

        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            momentum[i] = momentum[i] - stepSize / 2 * derivative[i];
        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                factors.setParameterValueQuietly(j, randel, factors.getParameterValue(j,randel)+stepSize*momentum[j]);
                factors.fireParameterChangedEvent(factors.getRowDimension()*randel+j,null);
            }


            if(i!=nSteps){
                derivative=getGradient(randel,mean,prec);

                for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                    momentum[j] = momentum[j] - stepSize * derivative[j];
                }
            }
        }

        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            momentum[i] = momentum[i] - stepSize / 2 * derivative[i];
        }

        double res=0;
        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            res+=momentum[i]*momentum[i]/2;
        }

        return res-prop;
    }
}
