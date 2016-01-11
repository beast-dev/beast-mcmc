package dr.evomodel.operators;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameter;
import dr.inference.operators.AbstractHamiltonianMCOperator;
import dr.inference.operators.CoercionMode;
import dr.inference.operators.OperatorFailedException;

/**
 * Created by max on 1/11/16.
 */
public class LoadingsHamiltonianMC extends AbstractHamiltonianMCOperator {
    private LatentFactorModel lfm;
    private MomentDistributionModel prior;
    private MatrixParameter factors;
    private MatrixParameter loadings;
    private MatrixParameter Precision;
    private int nfac;
    private int ntaxa;
    private int ntraits;
    private double stepSize;
    private int nSteps;


    public LoadingsHamiltonianMC(LatentFactorModel lfm, MomentDistributionModel prior, double weight, CoercionMode mode, double stepSize, int nSteps, double momentumSd){
        super(mode, momentumSd);
        setWeight(weight);
        this.lfm=lfm;
        this.prior=prior;
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
        return "LoadingsHamiltonianMC";
    }

    @Override
    public double doOperation() throws OperatorFailedException {


        double[][] derivative=getGradient();
        drawMomentum(lfm.getFactorDimension()*ntraits);

        double prop=0;
        for (int i = 0; i <momentum.length ; i++) {
            prop+=momentum[i]*momentum[i]/(2*getMomentumSd()*getMomentumSd());
        }

        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            for (int j = 0; j <ntraits ; j++) {
                momentum[i*ntraits+j] = momentum[i*ntraits+j] - stepSize / 2 * derivative[i][j];
            }

        }

        for (int i = 0; i <nSteps ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                for (int k = 0; k <ntraits ; k++) {
                    loadings.setParameterValueQuietly(j,k, loadings.getParameterValue(j,k)+stepSize*momentum[j*ntraits+k]);
                }

            }
            loadings.fireParameterChangedEvent();


            if(i!=nSteps){
                derivative=getGradient();

                for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                    for (int k = 0; k <ntraits ; k++) {
                        momentum[j*ntraits+k] = momentum[j*ntraits+k] - stepSize * derivative[j][k];
                    }

                }
            }
        }

        derivative=getGradient();
        for (int i = 0; i <lfm.getFactorDimension() ; i++) {
            for (int j = 0; j <ntraits ; j++) {
                momentum[i*ntraits+j] = momentum[i*ntraits+j] - stepSize / 2 * derivative[i][j];
            }

        }
        double res=0;
        for (int i = 0; i <momentum.length ; i++) {
            res+=momentum[i]*momentum[i]/(2*getMomentumSd()*getMomentumSd());
        }
        return prop-res;
    }



    private double[][] getLFMDerivative(){
        double[] residual=lfm.getResidual();
        double[][] answer= new double[ntraits][lfm.getFactorDimension()];
        for (int i = 0; i <ntaxa; i++) {
            for (int j = 0; j <ntraits; j++) {
                for (int k = 0; k <lfm.getFactorDimension() ; k++) {
                    answer[j][k]-=residual[i*ntaxa+j]*factors.getParameterValue(i,k);
                }
            }

        }
        for (int i = 0; i <ntraits ; i++) {
            for (int j = 0; j <lfm.getFactorDimension() ; j++) {
                answer[i][j]*=Precision.getParameterValue(i,i);
            }

        }
        return answer;
    }

    private double[][] getGradient(){
        double[][] answer=getLFMDerivative();
        for (int i = 0; i <loadings.getRowDimension() ; i++) {
            for (int j = 0; j <loadings.getColumnDimension() ; j++) {
                answer[i][j]-=2/loadings.getParameterValue(i,j)-(loadings.getParameterValue(i,j)-prior.getMean()[0])/prior.getScaleMatrix()[0][0];
            }

        }
        return answer;
    }

}
