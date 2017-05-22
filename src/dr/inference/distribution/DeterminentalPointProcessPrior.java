package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;

import java.util.List;
import java.util.Vector;

/**
 * Created by max on 4/6/16.
 */
public class DeterminentalPointProcessPrior extends AbstractModelLikelihood{

    double theta;
    MatrixParameterInterface data;
    boolean likelihoodKnown = false;
    boolean storedLikelihoodKnown;
    double logLikelihood;
    double storedLogLikelihood;
    double[][] relationshipList;
    double[][] storedRelationshipList;
    int size;
    int storedSize;
    Vector<Integer> changedList;
    Vector<Integer> storedChangedList;

    public DeterminentalPointProcessPrior(String name, double theta, MatrixParameterInterface data) {
        super(name);
        this.theta = theta;
        this.data = data;
        addVariable(data);
        relationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        storedRelationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        size = data.getColumnDimension();
        for (int i = 0; i < data.getRowDimension(); i++) {
            for (int j = 0; j < data.getColumnDimension(); j++) {
                if(i % (j+1) == 0) {
                    data.setParameterValueQuietly(i, j, 0);
                }
            }
        }

        reset();
        changedList = new Vector<Integer>();
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
        storedSize = size;
//        System.out.println("first");
//        for (int i = 0; i < relationshipList.length; i++) {
//            for (int j = 0; j < relationshipList.length ; j++) {
//                System.out.println(i + " + " + j + ": " + relationshipList[i][j]);
//            }
//        }


//        if(relationshipList.length != storedRelationshipList.length)
//            storedRelationshipList = new double[size][size];
//        for(int i = 0; i < relationshipList.length; i++)
//            storedRelationshipList[i] = relationshipList[i].clone();



//        System.out.println("stored");
//        for (int i = 0; i < relationshipList.length; i++) {
//            for (int j = 0; j < relationshipList.length ; j++) {
//                System.out.println(i + " + " + j + ": " + storedRelationshipList[i][j]);
//            }
//        }
//        System.arraycopy(relationshipList, 0, storedRelationshipList, 0, relationshipList.length);
        storedChangedList = (Vector<Integer>) changedList.clone();

    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        double[][] relationshipListTemp = relationshipList;
        relationshipList = storedRelationshipList;
        storedRelationshipList = relationshipListTemp;
        size = storedSize;
        Vector<Integer> changedListTemp = changedList;
        changedList = storedChangedList;
        storedChangedList = changedListTemp;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
        if(index == -1)
            reset();
        else {
            changedList.add(index);
        }
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        makeDirty();
            if(!likelihoodKnown) {
                logLikelihood = computeLogLikelihood();
                likelihoodKnown = true;
            }
        return logLikelihood;
    }

    public double computeLogLikelihood(){
        int newSize = data.getColumnDimension();
        while(!changedList.isEmpty()){
            int index = changedList.remove(0);
            int row = index % data.getRowDimension();
            int col = index / data.getRowDimension();
            for (int i = 0; i < data.getColumnDimension(); i++) {
                if(col != i){
                    if(data.getParameterValue(row, col) == data.getParameterValue(row, i)){
                        relationshipList[col][i] *= Math.exp(1 / (theta * theta));
                        relationshipList[i][col]=relationshipList[col][i];
                    }
                    else{
                        relationshipList[col][i] *= Math.exp(- 1 / (theta * theta));
                        relationshipList[i][col] = relationshipList[col][i];
                    }
                }
            }
        }

//        reset();

        if (newSize != size){
            size = newSize;
            relationshipList = new double[size][size];
                reset();
        }

        CholeskyDecomposition chol = null;
        try {
            chol = new CholeskyDecomposition(relationshipList);
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        double product = 0;
        for (int i = 0; i <newSize ; i++) {
            product += Math.log(chol.getL()[i][i]);
        }
        product *= 2;
        return product;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        changedList = new Vector<Integer>();
//        storedChangedList = new Vector<Integer>();
        reset();
    }

    public void reset(){
        if(relationshipList.length != data.getColumnDimension()){
            relationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        }
        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < i; j++) {
                int count = 0;
                for (int k = 0; k < data.getRowDimension(); k++) {
                    count += Math.abs(data.getParameterValue(k, i) - data.getParameterValue(k, j));
                }
                relationshipList[i][j] = Math.exp(- count / (theta * theta));
                relationshipList[j][i] = relationshipList[i][j];
            }
        }
        for (int i = 0; i <data.getColumnDimension() ; i++) {
                relationshipList[i][i] = 1;
        }
    }
}
