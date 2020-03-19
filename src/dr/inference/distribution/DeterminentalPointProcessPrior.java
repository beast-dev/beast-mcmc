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
    boolean[] notZero;
    boolean[] storedNotZero;
    int size;
    int storedSize;
    Vector<Integer> changedList;
    Vector<Integer> storedChangedList;
    Parameter normalizingConstants;
    boolean pathSampling;
    int sum;
    int storedSum;

    public DeterminentalPointProcessPrior(String name, double theta, MatrixParameterInterface data, Parameter normalizingConstants, boolean noZeros, boolean pathSampling, boolean resetData) {
        super(name);
        this.normalizingConstants = normalizingConstants;
        this.theta = theta;
        this.data = data;
        addVariable(data);
        relationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        storedRelationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        size = data.getColumnDimension();
        if(resetData){
            for (int i = 0; i < data.getRowDimension(); i++) {
                for (int j = 0; j < data.getColumnDimension(); j++) {
                    if(i % (j+1) == 0) {
                        data.setParameterValueQuietly(i, j, 0);
                    }
                }
            }
        }
        if(noZeros){
            notZero = new boolean[data.getColumnDimension()];
            storedNotZero = new boolean[data.getColumnDimension()];
        }
        reset();
        this.pathSampling = pathSampling;
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
//        System.out.println(notZero);
//        if(notZero != null && notZero.length != storedNotZero.length){
//            storedNotZero = new boolean[notZero.length];
//        }
//        if(notZero != null)
//            System.arraycopy(notZero, 0, storedNotZero, 0, storedNotZero.length);
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
        storedSum = sum;

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
//        boolean[] temp = notZero;
//        notZero = storedNotZero;
//        storedNotZero = temp;
        sum = storedSum;
    }

    @Override
    protected void acceptState() {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        likelihoodKnown = false;
        if(index == -1)
            reset();
//        else {
//            changedList.add(index);
//        }
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
                if(normalizingConstants != null) {
                    if(sum != 0)
                        logLikelihood += normalizingConstants.getParameterValue(sum - 1);
                    else
                        logLikelihood += normalizingConstants.getParameterValue(0);
                }
                likelihoodKnown = true;
            }
        return logLikelihood;
    }

    public double computeLogLikelihood(){
//        int newSize = data.getColumnDimension();
//        while(!changedList.isEmpty()){
//            int index = changedList.remove(0);
//            int row = index % data.getRowDimension();
//            int col = index / data.getRowDimension();
//            for (int i = 0; i < data.getColumnDimension(); i++) {
//                if(col != i){
//                    if(data.getParameterValue(row, col) == data.getParameterValue(row, i)){
//                        relationshipList[col][i] *= Math.exp(1 / (theta * theta));
//                        relationshipList[i][col]=relationshipList[col][i];
//                    }
//                    else{
//                        relationshipList[col][i] *= Math.exp(- 1 / (theta * theta));
//                        relationshipList[i][col] = relationshipList[col][i];
//                    }
//                }
//            }
//        }
//
////        reset();
//
//        if (newSize != size){
//            size = newSize;
//            relationshipList = new double[size][size];
//                reset();
//        }

        if(pathSampling && notZero != null && sum != data.getColumnDimension()){
            return Double.NEGATIVE_INFINITY;
        }

        CholeskyDecomposition chol = null;
        try {
            chol = new CholeskyDecomposition(relationshipList);
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
        double product = 0;
        for (int i = 0; i <relationshipList.length ; i++) {
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

    private void findZeros(){
        sum = 0;
        notZero = new boolean[data.getColumnDimension()];
        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getRowDimension(); j++) {
                if(data.getParameterValue(j, i) == 1)
                {notZero[i] = true;
                    sum++;
                    break;}
            }
        }

    }

    public int getSum(){
        findZeros();
        return sum;
    }


    public void reset(){
        if(notZero != null){
            findZeros();
            if(relationshipList.length != sum && sum != 0){
                relationshipList = new double[sum][sum];
            }
            else if(sum == 0)
                relationshipList = new double[1][1];
        }
        else if(relationshipList.length != data.getColumnDimension() || relationshipList == null){
            relationshipList = new double[data.getColumnDimension()][data.getColumnDimension()];
        }
        int offset1 = 0;
        for (int i = 0; i < relationshipList.length; i++) {
            int offset2 = 0;
            if((notZero != null && notZero[i] == true) || notZero == null){
            for (int j = 0; j < i; j++) {
                if((notZero !=null && notZero[j] == true) || notZero == null){
                    int count = 0;
                    for (int k = 0; k < data.getRowDimension(); k++) {
                        count += Math.abs(data.getParameterValue(k, i) - data.getParameterValue(k, j));
                    }
                    relationshipList[i + offset1][j + offset2] = Math.exp(- count / (theta * theta));
                    relationshipList[j + offset2][i + offset1] = relationshipList[i + offset1][j + offset2];
                }
                else
                    offset2 --;
            }
            }
            else
                offset1 --;
        }
        for (int i = 0; i <relationshipList.length ; i++) {
                relationshipList[i][i] = 1;
        }
    }
}
