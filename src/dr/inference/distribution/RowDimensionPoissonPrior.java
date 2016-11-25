package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.distributions.PoissonDistribution;

/**
 * Created by maxryandolinskytolkoff on 7/20/16.
 */
public class RowDimensionPoissonPrior extends AbstractModelLikelihood implements MatrixSizePrior {
    public RowDimensionPoissonPrior(double untruncatedMean, AdaptableSizeFastMatrixParameter parameter){
        super("");
        this.poisson = new PoissonDistribution(untruncatedMean);
        this.parameter = parameter;
        addVariable(parameter);
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
//        System.out.println(poisson.logPdf(parameter.getRowDimension()) - Math.log(1 - Math.exp(-poisson.mean())));
        return poisson.logPdf(parameter.getRowDimension()) - Math.log(1 - Math.exp(-poisson.mean()));
    }

    @Override
    public void makeDirty() {

    }

    @Override
    public double getSizeLogLikelihood() {
        return getLogLikelihood();
    }

    PoissonDistribution poisson;
    AdaptableSizeFastMatrixParameter parameter;
}
