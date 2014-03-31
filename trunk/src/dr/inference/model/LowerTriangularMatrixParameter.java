package dr.inference.model;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 3/27/14
 * Time: 3:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class LowerTriangularMatrixParameter extends CompoundParameter {
    int rowDim;
    int colDim;
    Parameter[] parameters;


    public LowerTriangularMatrixParameter(String name, Parameter[] params, int rowDim, int colDim) {
        super(name, params);

        for(int i=0; i<params.length; i++){
            if(colDim-rowDim+i<0)
            {params[i].setDimension(0);}
            else
            {params[i].setDimension(colDim-rowDim+i+1);}
        }


        this.parameters=params;
        this.rowDim=rowDim;
        this.colDim=colDim;
    }
}
