package dr.inference.operators;

import dr.inference.model.MatrixParameterInterface;
import dr.math.MathUtils;

public class LoadingsSparsityOperator extends SimpleMCMCOperator{
    public LoadingsSparsityOperator(double weight, LoadingsGibbsTruncatedOperator loadingsGibbs, MatrixParameterInterface sparse){

        setWeight(weight);
        this.loadingsGibbs = loadingsGibbs;
        this.sparse = sparse;
    }

    @Override
    public String getPerformanceSuggestion() {

        return null;
    }

    @Override
    public String getOperatorName() {
        return "LoadingsSparsityOperator";
    }

    @Override
    public double doOperation() {
        int row = MathUtils.nextInt(sparse.getRowDimension());
        int col = MathUtils.nextInt(sparse.getColumnDimension());

        double hastings = 0;
        hastings += loadingsGibbs.drawI(row, col, false);
        if(sparse.getParameterValue(row, col) == 0)
            sparse.setParameterValue(row, col, 1);
        else
            sparse.setParameterValue(row, col, 0);
        hastings -= loadingsGibbs.drawI(row, col, true);

        if(Double.isNaN(hastings)){
            System.out.println("is NaN");
        }
        if(Double.isInfinite(hastings))
        {
            return Double.NEGATIVE_INFINITY;
        }
        return hastings;
    }


    LoadingsGibbsTruncatedOperator loadingsGibbs;
    MatrixParameterInterface sparse;
}
