package dr.inference.distribution;

import dr.inference.model.*;

/**
 * @author Max Tolkoff
 */
public class RowDimensionMultinomialPrior extends AbstractModelLikelihood implements MatrixSizePrior {
    AdaptableSizeFastMatrixParameter data;
    Parameter distribution;
    final boolean transpose;

    public RowDimensionMultinomialPrior(String id, AdaptableSizeFastMatrixParameter data, Parameter distribution, boolean transpose){
        super(id);
        addVariable(data);
        addVariable(distribution);
        this.data = data;
        this.distribution = distribution;
        this.transpose = transpose;
    }

    @Override
    public double getSizeLogLikelihood() {
        if(!transpose){
            return distribution.getParameterValue(data.getColumnDimension() - 1);
        }
        else
        {
            return distribution.getParameterValue(data.getRowDimension() - 1);
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return getSizeLogLikelihood();
    }

    @Override
    public void makeDirty() {
        getSizeLogLikelihood();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

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
}
