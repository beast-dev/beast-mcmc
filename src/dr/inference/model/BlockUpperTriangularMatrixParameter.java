package dr.inference.model;

/*
@author Max Tolkoff
*/

public class BlockUpperTriangularMatrixParameter extends MatrixParameter {
    private int rowDim;

    public TransposedBlockUpperTriangularMatrixParameter transposeBlock(){
        return TransposedBlockUpperTriangularMatrixParameter.recast(getVariableName(), this);
    }

    public BlockUpperTriangularMatrixParameter(String name, Parameter[] params) {
        super(name);

        int rowDim=params[params.length-1].getSize();
        int colDim=params.length;

        for(int i=0; i<colDim; i++){
            if(i<rowDim)
            {params[i].setDimension(i+1);
            this.addParameter(params[i]);}
            else
            {params[i].setDimension(rowDim);
                this.addParameter(params[i]);
//                System.err.print(colDim-rowDim+i+1);
//                System.err.print("\n");
            }
        }
        this.rowDim=rowDim;
//        System.err.print("Dimension Check\n");
//        System.err.print(rowDim);
//        System.err.print("\n");
//        System.err.print(colDim);
//        System.err.print("\n");

//
//        double[][] temp=getParameterAsMatrix();
//        for(int i=0; i<getRowDimension(); i++){
//            for(int j=0; j<getColumnDimension(); j++)
//            {
//                System.err.print(temp[i][j]);
//                System.err.print(" ");
//            }
//            System.err.print("\n");
//        }
    }

    @Override
    public int getRowDimension() {
        return rowDim;
    }

    public void setRowDimension(int rowDim){
        this.rowDim=rowDim;
    }

//    public double[][] getParameterAsMatrix(){
//        double[][] answer=new double[getRowDimension()][getColumnDimension()];
//        for(int i=0; i<getRowDimension(); i++){
//            for(int j=0; j<getColumnDimension(); j++){
//                if(i<=j){
////                    System.err.print(parameters[i].getSize());
////                    System.err.print("we get here\n");
//                    answer[i][j]=getParameter(j).getParameterValue(i);
//                }
//                else{
////                    System.err.print(i);
////                    System.err.print(" ");
////                    System.err.print(j);
////                    System.err.print("\n");
//                    answer[i][j]=0;
////                    System.err.print("getting here?\n");
//                }
//            }
//        }
//
//        return answer;
//    }

    public double getParameterValue(int row, int col) {
        if (row > col) {
            return 0.0;
        } else {
            return getParameter(col).getParameterValue(row);
        }
    }



}
