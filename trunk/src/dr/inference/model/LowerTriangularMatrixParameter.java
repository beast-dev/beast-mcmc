package dr.inference.model;

/*
@author Max Tolkoff
*/

public class LowerTriangularMatrixParameter extends MatrixParameter {
    private int rowDim;
    private int colDim;
    private Parameter[] parameters;


    public LowerTriangularMatrixParameter(String name, Parameter[] params, int rowDim, int colDim) {
        super(name, params);

        for(int i=0; i<rowDim; i++){
            if(colDim-rowDim+i<0)
            {params[i].setDimension(0);}
            else
            {params[i].setDimension(colDim-rowDim+i+1);
//                System.err.print(colDim-rowDim+i+1);
//                System.err.print("\n");
            }
        }

//        System.err.print("Dimension Check\n");
//        System.err.print(rowDim);
//        System.err.print("\n");
//        System.err.print(colDim);
//        System.err.print("\n");
        this.parameters=params;
        this.rowDim=rowDim;
        this.colDim=colDim;

//        double[][] temp=getParameterAsMatrix();
//        for(int i=0; i<rowDim; i++){
//            for(int j=0; j<colDim; j++)
//            {
//                System.err.print(temp[i][j]);
//                System.err.print(" ");
//            }
//            System.err.print("\n");
//        }
    }

    public double[][] getParameterAsMatrix(){
        double[][] answer=new double[rowDim][colDim];
        for(int i=0; i<rowDim; i++){
            for(int j=0; j<colDim; j++){
                if(parameters[i].getSize()>j){
//                    System.err.print(parameters[i].getSize());
//                    System.err.print("we get here\n");
                    answer[i][j]=parameters[i].getParameterValue(j);
                }
                else{
//                    System.err.print(i);
//                    System.err.print(" ");
//                    System.err.print(j);
//                    System.err.print("\n");
                    answer[i][j]=0;
//                    System.err.print("getting here?\n");
                }
            }
        }

        return answer;
    }
    public int getRowDim()
    {return rowDim;}

    public int getColDim()
    {return colDim;}


    //Which parameter is changed is a function of the length of the parameter
    public void setParameterValue(double[] New){
        int len=New.length;
        int paramNum=colDim-rowDim+len-1;
        for(int i=0; i<paramNum; i++){
            parameters[paramNum].setParameterValueQuietly(i, New[i]);
        }
        parameters[paramNum].fireParameterChangedEvent();
    }

    public double[] getParameterValues(int j){
        double[] answer=parameters[j].getParameterValues();
        return answer;
    }
}
