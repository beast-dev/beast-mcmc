package dr.inference.model;

/**
 * Created by max on 11/30/15.
 */
public class ElementWiseMatrixMultiplicationParameter extends MatrixParameter {
    private MatrixParameterInterface[] paramList;

    public ElementWiseMatrixMultiplicationParameter(String name) {
        super(name);
    }

    public ElementWiseMatrixMultiplicationParameter(String name, MatrixParameterInterface[] matList) {
        super(name);
        this.paramList =matList;
        for (MatrixParameterInterface mat : matList) {
            mat.addVariableListener(this);
        }
    }

    @Override
    public double getParameterValue(int dim) {
        double prod = 1;
        for (int i = 0; i < paramList.length ; i++) {
            prod = prod * paramList[i].getParameterValue(dim);
        }
        return prod;
    }

    public double getParameterValue(int row, int col){
        double prod=1;
        for (int i = 0; i < paramList.length ; i++) {
            prod=prod * paramList[i].getParameterValue(row, col);
        }
        return prod;
    }


    protected void storeValues() {
        for (Parameter p : paramList) {
            p.storeParameterValues();
        }
    }

    protected void restoreValues() {
        for (Parameter p : paramList) {
            p.restoreParameterValues();
        }
    }

    protected void acceptValues() {
        for (Parameter p : paramList) {
            p.acceptParameterValues();
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index, type);
    }

//    @Override
//    public void fireParameterChangedEvent(){
//        for (int i = 0; i < paramList.length ; i++) {
//            paramList[i].fireParameterChangedEvent();
//        }
//    }
//
//    @Override
//    public void fireParameterChangedEvent(int index, ChangeType type){
//        for (int i = 0; i < paramList.length; i++) {
//            paramList[i].fireParameterChangedEvent(index, type);
//        }
//    }

    @Override
    public int getDimension() {
        return paramList[0].getDimension();
    }

    @Override
    public int getColumnDimension() {
        return paramList[0].getColumnDimension();
    }

    public int getRowDimension(){
        return paramList[0].getRowDimension();
    }
}
