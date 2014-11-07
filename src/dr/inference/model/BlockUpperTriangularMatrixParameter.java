package dr.inference.model;

/*
@author Max Tolkoff
*/

public class BlockUpperTriangularMatrixParameter extends MatrixParameter {
    private int rowDim;
    private Bounds bounds = null;

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

    private int getRow(int PID){
        return  PID%getColumnDimension();
    }

    private int getColumn(int PID){
        return PID/getRowDimension();
    }

    public void setParameterValue(int PID, double value){
        int row=getRow(PID);
        int col=getColumn(PID);

        if(row<=col){
            setParameterValue(row, col, value);
        }
    }

    public double getParameterValue(int id){
        int row=getRow(id);
        int col=getColumn(id);

        if(row<=col){
            return getParameterValue(row, col);
        }
        else
        {
            return 0;
        }
    }

    public void addBounds(Bounds<Double> boundary) {

        if (bounds == null) {
            bounds = new BUTMPBounds();
//            return;
        } //else {
        IntersectionBounds newBounds = new IntersectionBounds(getDimension());
        newBounds.addBounds(bounds);

//        }
        ((IntersectionBounds) bounds).addBounds(boundary);
    }

    public Bounds<Double> getBounds() {

        if (bounds == null) {
            bounds = new BUTMPBounds();
        }
        return bounds;
    }


    private class BUTMPBounds implements Bounds<Double>{


            public Double getUpperLimit(int dim) {
                int row=getRow(dim);
                int col=getColumn(dim);

                if(row <= col)
                return getParameters().get(dim-row).getBounds().getUpperLimit(getPindex().get(dim-row));
                else
                    return 0.0;
            }

            public Double getLowerLimit(int dim) {
                int row=getRow(dim);
                int col=getColumn(dim);

                if(row <= col)
                return getParameters().get(dim-row).getBounds().getLowerLimit(getPindex().get(dim-row));
                else
                    return 0.0;
            }

            public int getBoundsDimension() {
                int nBlanks = 0;
                for (int i = 0; i <getColumnDimension() ; i++) {
                    nBlanks+=1;
                }

                return getDimension()-nBlanks;
        }
    }

    public int getDimension(){
        return getRowDimension()*getColumnDimension();
    }


}
