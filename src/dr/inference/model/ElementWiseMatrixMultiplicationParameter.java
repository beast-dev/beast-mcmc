package dr.inference.model;

/**
 * Created by max on 11/30/15.
 */
public class ElementWiseMatrixMultiplicationParameter extends MatrixParameter {
    private MatrixParameter[] matlist;

    public ElementWiseMatrixMultiplicationParameter(String name) {
        super(name);
    }

    public ElementWiseMatrixMultiplicationParameter(String name, MatrixParameter[] matList) {
        super(name);
        this.matlist=matList;
    }

    @Override
    public double getParameterValue(int dim) {
        double prod=1;
        for (int i = 0; i <matlist.length ; i++) {
            prod=prod*matlist[i].getParameterValue(dim);
        }
        return prod;
    }

    public double getParameterValue(int row, int col){
        double prod=1;
        for (int i = 0; i <matlist.length ; i++) {
            prod=prod*matlist[i].getParameterValue(row,col);
        }
        return prod;
    }

    @Override
    public int getDimension() {
        return matlist[0].getDimension();
    }

    @Override
    public int getColumnDimension() {
        return matlist[0].getColumnDimension();
    }

    public int getRowDimension(){
        return matlist[0].getRowDimension();
    }
}
