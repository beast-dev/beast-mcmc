package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;

import java.util.List;
import java.util.Vector;

/**
 * Created by max on 4/6/16.
 */
public class DeterminentalPointProcessPrior extends AbstractModelLikelihood{

    double theta;
    AdaptableSizeFastMatrixParameter data;
    boolean likelihoodKnown = false;
    boolean storedLikelihoodKnown;
    double logLikelihood;
    double storedLogLikelihood;
    double[][] relationshipList;
    double[][] storedRelationshipList;
    int size;
    int storedSize;
    Vector<Integer> changedList;

    public DeterminentalPointProcessPrior(String name, double theta, AdaptableSizeFastMatrixParameter data) {
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

        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getColumnDimension(); j++) {
                int count = 0;
                for (int k = 0; k < data.getRowDimension(); k++) {
                    count += Math.abs(data.getParameterValue(k, i) - data.getParameterValue(k, j));
                }
                relationshipList[i][j] = Math.exp(- count / (theta * theta));
                relationshipList[j][i] = relationshipList[i][j];
            }
        }
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
        if(relationshipList.length != storedRelationshipList.length)
            storedRelationshipList = new double[size][size];
        System.arraycopy(relationshipList, 0, storedRelationshipList, 0, relationshipList.length);

    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        double[][] relationshipListTemp = relationshipList;
        relationshipList = storedRelationshipList;
        storedRelationshipList = relationshipListTemp;
        size = storedSize;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
        changedList.add(index);
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
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
                        relationshipList[i][col]=relationshipList[i][col];
                    }
                    else{
                        relationshipList[col][i] *= Math.exp(- 1 / (theta * theta));
                        relationshipList[i][col] = relationshipList[col][i];
                    }
                }
                
            }
        }

        if (newSize != size){
            double[][] relationshipListTemp = new double[size][size];
            for (int i = 0; i < Math.min(newSize, size) ; i++) {
                for (int j = 0; j < Math.min(newSize, size); j++) {
                    relationshipListTemp[i][j] = relationshipList[i][j];
                }
            }
            if(relationshipList.length < relationshipListTemp.length){
                for (int i = 0; i < data.getColumnDimension(); i++) {
                    int count = 0;
                    for (int j = 0; j < data.getRowDimension(); j++) {
                        count += Math.abs(data.getParameterValue(i, j) - data.getParameterValue(newSize - 1, j));
                    }
                    relationshipListTemp[i][newSize - 1] = Math.exp(count / theta);
                    relationshipListTemp[newSize-1][i] = relationshipListTemp[i][newSize - 1];
                }
            }
            size = newSize;
            relationshipList = relationshipListTemp;
        }

        CholeskyDecomposition chol = null;
        try {
            chol = new CholeskyDecomposition(relationshipList);
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        double product = 1;
        for (int i = 0; i <newSize ; i++) {
            product *= chol.getL()[i][i];
        }
        logLikelihood = Math.log(product);
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getColumnDimension(); j++) {
                int count = 0;
                for (int k = 0; k < data.getRowDimension(); k++) {
                    count += Math.abs(data.getParameterValue(k, i) - data.getParameterValue(k, j));
                }
                relationshipList[i][j] = Math.exp(- count / (theta * theta));
            }
        }
    }
}
