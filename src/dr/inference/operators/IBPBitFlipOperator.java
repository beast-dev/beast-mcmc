package dr.inference.operators;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.IndianBuffetProcessPrior;
import dr.inference.model.Likelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

/**
 * Created by maxryandolinskytolkoff on 11/21/16.
 */
public class IBPBitFlipOperator extends BitFlipOperator {
    public IBPBitFlipOperator(Parameter parameter, double weight, IndianBuffetProcessPrior IBP) {
        super(parameter, weight, true);
        this.sparsity = (MatrixParameterInterface) parameter;
        this.IBP = IBP;
    }

    @Override
    public double doOperation(Likelihood likelihood)
    { //throws OperatorFailedException {
        return super.doOperation(likelihood);
    }

    protected double sum(int pos) {
        double sum = 0;
        int column = pos / sparsity.getRowDimension();
        for (int i = 0; i < sparsity.getRowDimension(); i++) {
            sum += sparsity.getParameterValue(i, column);
        }


        return sum;
    }

    protected  int getDimension(){
        return sparsity.getRowDimension();
    }

    IndianBuffetProcessPrior IBP;
    MatrixParameterInterface sparsity;


}
