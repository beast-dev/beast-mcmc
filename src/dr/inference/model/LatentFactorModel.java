/*
 * LatentFactorModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

import dr.inference.distribution.LatentFactorModelInterface;
import dr.math.matrixAlgebra.Matrix;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;


/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class LatentFactorModel extends AbstractModelLikelihood implements Citable, LatentFactorModelInterface {
//    private Matrix data;
//    private Matrix factors;
//    private Matrix loadings;

    private final MatrixParameterInterface data;
    private final MatrixParameterInterface factors;
    private final MatrixParameterInterface loadings;
    private MatrixParameterInterface sData;
    private final DiagonalMatrix rowPrecision;
    private final DiagonalMatrix colPrecision;
    private final Parameter continuous;

    private final boolean scaleData;

    private final int dimFactors;
    private final int dimData;
    private final int nTaxa;

    private boolean newModel;

    private boolean likelihoodKnown = false;
    private boolean isDataScaled = false;
    private boolean storedLikelihoodKnown;
    private boolean residualKnown = false;
    private boolean dataKnown = false;
    private boolean storedDataKnown;
    private boolean LxFKnown = false;
    private boolean storedResidualKnown = false;
    private boolean storedLxFKnown;
    private boolean traceKnown = false;
    private boolean storedTraceKnown;
    private boolean logDetColKnown = false;
    private boolean storedLogDetColKnown;
    private double trace;
    private double storedTrace;
    private double logLikelihood;
    private double storedLogLikelihood;
    private double logDetCol;
    private double storedLogDetCol;
    private boolean[][] changed;
    private boolean[][] storedChanged;
    private boolean RecomputeResiduals;
    private boolean RecomputeFactors;
    private boolean RecomputeLoadings;
    private Vector<Integer> changedValues;
    private Vector<Integer> storedChangedValues;
    private boolean factorsKnown = false;
    private boolean storedFactorsKnown = false;
    private boolean loadingsKnown = false;
    private boolean storedLoadingsKnown = false;
    private boolean totalRecompute = true;
    private boolean storedTotalRecompute = false;

    private double[] residual;
    private double[] LxF;
    private double[] storedResidual;
    private double[] storedLxF;

    private double pathParameter = 1.0;

    private final Parameter missingIndicator;
    private final int[] rowCount;
    private final int nmeasurements;

    public LatentFactorModel(MatrixParameterInterface data, MatrixParameterInterface factors, MatrixParameterInterface loadings,
                             DiagonalMatrix rowPrecision, DiagonalMatrix colPrecision,
                             Parameter missingIndicator,
                             boolean scaleData, Parameter continuous, boolean newModel, boolean recomputeResiduals,
                             boolean recomputeFactors, boolean recomputeLoadings
    ) {
        super("");
        this.RecomputeResiduals = recomputeResiduals;
        this.RecomputeFactors = recomputeFactors;
        this.RecomputeLoadings = recomputeLoadings;
        changedValues = new Vector<Integer>();
        for (int i = 0; i < data.getDimension(); i++) {
            changedValues.add(i);
        }
        storedChangedValues = new Vector<Integer>();
//        data = new Matrix(dataIn.getParameterAsMatrix());
//        factors = new Matrix(factorsIn.getParameterAsMatrix());
//        loadings = new Matrix(loadingsIn.getParameterAsMatrix());
        this.newModel = newModel;
        this.scaleData = scaleData;
        this.data = data;
        this.factors = factors;
        // Put default bounds on factors
        if(factors instanceof MatrixParameter){
            for (int i = 0; i < factors.getColumnDimension(); ++i) {
                Parameter p = factors.getParameter(i);
                System.err.println(p.getId() + " " + p.getDimension());
                p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, p.getDimension()));
            }
        }
        this.continuous=continuous;

        this.loadings = loadings;

        this.missingIndicator = missingIndicator;
        rowCount = new int[colPrecision.getRowDimension()];
        if(missingIndicator != null){
            for (int i = 0; i < data.getRowDimension(); i++) {
                for (int j = 0; j < data.getColumnDimension(); j++) {
                    if(missingIndicator.getParameterValue(j * data.getRowDimension() + i) == 0){
                        rowCount[i]++;
                    }
                }

            }
        }
        else{
            for (int i = 0; i < data.getRowDimension(); i++) {
                rowCount[i] = data.getColumnDimension();
            }

        }
        int nmeasurements = 0;
        for (int i = 0; i < rowCount.length ; i++) {
            nmeasurements += rowCount[i];
        }
        this.nmeasurements = nmeasurements;

//        storedData=new MatrixParameter(null);
//        for (int i = 0; i <continuous.getDimension(); i++) {
//            if(continuous.getParameterValue(i)==0)
//                storedData.addParameter(new Parameter.Default(data.getColumnDimension()));
//
//        }

        // Put default bounds on loadings
//        loadings.addBounds();


        changed = new boolean[loadings.getRowDimension()][factors.getColumnDimension()];
        storedChanged = new boolean[loadings.getRowDimension()][factors.getColumnDimension()];

        for (int i = 0; i < loadings.getRowDimension(); i++) {
            for (int j = 0; j < factors.getColumnDimension(); j++) {
                changed[i][j] = true;
            }
        }

        this.rowPrecision = rowPrecision;
        this.colPrecision = colPrecision;

        addVariable(data);
        addVariable(factors);
        addVariable(loadings);
        addVariable(rowPrecision);
        addVariable(colPrecision);

        dimFactors = factors.getRowDimension();
        dimData = loadings.getRowDimension();
//        nTaxa = factors.getParameterCount();
//        nTaxa = factors.getParameter(0).getDimension();
        nTaxa = factors.getColumnDimension();

//        System.out.print(nTaxa);
//        System.out.print("\n");
//        System.out.print(dimData);
//        System.out.print("\n");
//        System.out.println(dimFactors);
//        System.out.println(data.getDimension());
//        System.out.println(data.getRowDimension());
//        System.out.println(data.getColumnDimension());

//        System.out.println(new Matrix(data.getParameterAsMatrix()));
//        System.out.println(new Matrix(factors.getParameterAsMatrix()));


        if (nTaxa != data.getColumnDimension()) {
            throw new RuntimeException("DATA COLUMNS MUST HAVE THE SAME DIMENSION AS FACTOR COLUMNS\n");
//            System.exit(10);
        }
        if (dimData != data.getRowDimension()) {
            System.out.println(dimData);
            System.out.println(data.getRowDimension());
            System.out.println(loadings.getRowDimension());
            throw new RuntimeException("DATA ROWS MUST HAVE THE SAME DIMENSION AS LOADINGS ROWS\n");
//            System.exit(10);
        }
        if (factors.getRowDimension() != loadings.getColumnDimension()) {
            System.out.println(getModelName());
            throw new RuntimeException("LOADINGS AND FACTORS MUST HAVE THE SAME NUMBER OF FACTORS\n");
//            System.exit(10);
        }
        if (dimData < dimFactors) {
            throw new RuntimeException("MUST HAVE FEWER FACTORS THAN DATA POINTS\n");
        }

        residual = new double[loadings.getRowDimension() * factors.getColumnDimension()];
        LxF = new double[loadings.getRowDimension() * factors.getColumnDimension()];
        storedResidual = new double[residual.length];
        storedLxF = new double[LxF.length];

        if (!isDataScaled & !scaleData) {
            sData = this.data;
            isDataScaled = true;
        }
        if (!isDataScaled) {
            sData = computeScaledData();
            isDataScaled = true;
            for (int i = 0; i < sData.getRowDimension(); i++) {
                for (int j = 0; j < sData.getColumnDimension(); j++) {
                    this.data.setParameterValue(i, j, sData.getParameterValue(i, j));
//                    System.out.println(this.data.getParameterValue(i,j));
                }

            }
            data.fireParameterChangedEvent();
        }

        double sum = 0;
        for (int i = 0; i < sData.getRowDimension(); i++) {
            for (int j = 0; j < sData.getColumnDimension(); j++) {
                if (continuous.getParameterValue(i) == 1 && sData.getParameterValue(i, j) != 0) {
                    sum += -.5 * Math.log(2 * StrictMath.PI) - .5 * sData.getParameterValue(i, j) * sData.getParameterValue(i, j);
                }
            }
        }
        System.out.println("Constant Value for Path Sampling (normal 0,1): " + -1 * sum);

        computeResiduals();
//        System.out.print(new Matrix(residual.toComponents()));
//        System.out.print(calculateLogLikelihood());
    }


    //    public Matrix getData(){
//        Matrix ans=data;
//        return ans;
//    }
//
//    public Matrix getFactors(){
//        Matrix ans=factors;
//        return ans;
//    }
//
//    public Matrix getLoadings(){
//        Matrix ans=loadings;
//        return ans;
//    }
//
//    public Matrix getResidual(){
//        Matrix ans=residual;
//        return ans;
//    }
    public MatrixParameterInterface getFactors(){
        return factors;
    }

    public MatrixParameter getColumnPrecision() {
        return colPrecision;
    }

    public MatrixParameterInterface getLoadings(){
        return loadings;
    }

    public MatrixParameterInterface getData(){
        return data;
    }

    public Parameter returnIntermediate() {
        if (!residualKnown && checkLoadings()) {
            computeResiduals();
        }
        return data;
    }

//    public Parameter returnIntermediate(int PID)
//    {   //residualKnown=false;
//        if(!residualKnown && checkLoadings()){
//        computeResiduals();
//        }
//    return data.getParameter(PID);
//    }


    public MatrixParameterInterface getScaledData(){
        return data;
    }

    public Parameter getContinuous() {
        return continuous;
    }

    public int getFactorDimension() {
        return factors.getRowDimension();
    }

    public double[] getResidual() {
        computeResiduals();
        return residual;
    }

    public Parameter getMissingIndicator() {
        return missingIndicator;
    }

    public int getRowCount(int i){
        return rowCount[i];
    }

    private void Multiply(MatrixParameterInterface Left, MatrixParameterInterface Right, double[] answer){
        int dim = Left.getColumnDimension();
        int n = Left.getRowDimension();
        int p = Right.getColumnDimension();

        if(((factorsKnown == false && !RecomputeFactors) || (!dataKnown && !RecomputeResiduals) || (!loadingsKnown && !RecomputeLoadings)) && !totalRecompute){
            double sum;
            ListIterator<Integer> li = changedValues.listIterator();
            while (li.hasNext()) {
                int index = li.next();
                int i = index % n;
                int j = index / n;
                if(missingIndicator == null || missingIndicator.getParameterValue(j * n + i) != 1){
                    sum = 0;
                    for (int k = 0; k < dim; k++) {
//                System.out.println(data.getColumnDimension());
//                System.out.println(index);
                        sum += Left.getParameterValue(i, k) *
                                Right.getParameterValue(k, j);
                    }
                    answer[j * n + i] = sum;
                }
            }
        } else {
            for (int i = 0; i < answer.length; i++) {
                answer[i] = 0;
            }
            for (int i = 0; i < n; i++) {
                for (int k = 0; k < dim; k++){
                    if (Left.getParameterValue(i, k) != 0){
                        for (int j = 0; j < p; j++) {
                            if(missingIndicator == null || missingIndicator.getParameterValue(j * n + i) != 1){
                                if ((changed[i][j] == true && continuous.getParameterValue(i) != 0) || newModel) {
                                    double sum = 0;

                                    sum += Left.getParameterValue(i, k) * Right.getParameterValue(k, j);
                                    answer[j * n + i] += sum;
                                    //changed[i][j]=false;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void add(MatrixParameter Left, MatrixParameter Right, double[] answer) {
        int row = Left.getRowDimension();
        int col = Left.getColumnDimension();
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                answer[i * col + j] = Left.getParameterValue(i, j) + Right.getParameterValue(i, j);
            }

        }
    }

    private void subtract(MatrixParameterInterface Left, double[] Right, double[] answer) {
        int row = Left.getRowDimension();
        int col = Left.getColumnDimension();
        if(((!RecomputeResiduals && !dataKnown) || (!RecomputeFactors && !factorsKnown) || (!RecomputeLoadings && !loadingsKnown)) && !totalRecompute) {
            int size = changedValues.size();
            for (int i = 0; i < size; i++) {
                int id = changedValues.get(i);
                int tcol=id / row;
                int trow=id % row;
//                System.out.println(Left.getParameterValue(id)==Left.getParameterValue(tcol,trow));
                if(missingIndicator == null || missingIndicator.getParameterValue(tcol * row + trow) != 1)
                    answer[tcol * row + trow] = Left.getParameterValue(id) - Right[tcol * row + trow];
            }
        } else {
            for (int i = 0; i < row; i++) {
                if (continuous.getParameterValue(i) != 0 || newModel) {
                    for (int j = 0; j < col; j++) {
                        if(missingIndicator == null || missingIndicator.getParameterValue(j * row + i) != 1)
                            answer[j * row + i] = Left.getParameterValue(i, j) - Right[j * row + i];
                    }
                }
//              else{
//                  for (int j = 0; j <col; j++) {
//                        Left.setParameterValueQuietly(i,j, Right[i*col+j]);
//                  }
//                    containsDiscrete=true;
//                }

            }
        }
        changedValues.clear();
//        if(containsDiscrete){
//            Left.fireParameterChangedEvent();}
    }

    private double TDTTrace(double[] array, DiagonalMatrix middle) {
        int innerDim = middle.getRowDimension();
        int outerDim = array.length / innerDim;
        double sum = 0;
        for (int j = 0; j < innerDim; j++) {
            if (continuous.getParameterValue(j) != 0 || newModel) {
                for (int i = 0; i < outerDim; i++) {
                    if(missingIndicator == null || missingIndicator.getParameterValue(i * innerDim + j) != 1){
                        double s1 = array[i * innerDim + j];
                        double s2 = middle.getParameterValue(j, j);
                        sum += s1 * s1 * s2;
                    }
                }
            }
        }
        return sum;
    }


    private MatrixParameter computeScaledData() {
        MatrixParameter answer = new MatrixParameter(data.getParameterName() + ".scaled");
        answer.setDimensions(data.getRowDimension(), data.getColumnDimension());
        //       Matrix answer=new Matrix(data.getRowDimension(), data.getColumnDimension());
        double[][] aData = data.getParameterAsMatrix();
        double[] meanList = new double[data.getRowDimension()];
        double[] varList = new double[data.getRowDimension()];
        double[] count = new double[data.getRowDimension()];
        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getRowDimension(); j++) {
                if (data.getParameterValue(j, i) != 0) {
                    meanList[j] += data.getParameterValue(j, i);
                    count[j]++;
                }
            }
        }
        for (int i = 0; i < data.getRowDimension(); i++) {
            if (continuous.getParameterValue(i) == 1)
                meanList[i] = meanList[i] / count[i];
            else
                meanList[i] = 0;
        }

        double[][] answerTemp = new double[data.getRowDimension()][data.getColumnDimension()];
        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getRowDimension(); j++) {
                if (aData[j][i] != 0) {
                    answerTemp[j][i] = aData[j][i] - meanList[j];
                }
            }
        }
//        System.out.println(new Matrix(answerTemp));

        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getRowDimension(); j++) {
                varList[j] += answerTemp[j][i] * answerTemp[j][i];
            }
        }

        for (int i = 0; i < data.getRowDimension(); i++) {
            if (continuous.getParameterValue(i) == 1) {
                varList[i] = varList[i] / (count[i] - 1);
                varList[i] = StrictMath.sqrt(varList[i]);
            } else {
                varList[i] = 1;
            }
        }
//        System.out.println(data.getColumnDimension());
//        System.out.println(data.getRowDimension());

        for (int i = 0; i < data.getColumnDimension(); i++) {
            for (int j = 0; j < data.getRowDimension(); j++) {
                answer.setParameterValue(j, i, answerTemp[j][i] / varList[j]);
            }
        }
//        System.out.println(new Matrix(answerTemp));
//        computeResiduals();
        return answer;
    }

    private Matrix copy(CompoundParameter parameter, int dimMajor, int dimMinor) {
        return new Matrix(parameter.getParameterValues(), dimMajor, dimMinor);
    }

    private void computeResiduals() {
//    LxFKnown=false;


//        if(firstTime || (!factorVariablesChanged.empty() && !loadingVariablesChanged.empty())){
        if (!LxFKnown) {
            Multiply(loadings, factors, LxF);

        }
        subtract(data, LxF, residual);
        LxFKnown = true;
        residualKnown = true;
        factorsKnown = true;
        loadingsKnown = true;
        dataKnown = true;
        totalRecompute = false;
//        firstTime=false;}
//        else{
//            while(!factorVariablesChanged.empty()){
//
//            }
//            while(!loadingVariablesChanged.empty()){
//
//            }
//        }
//


    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
//        data.storeParameterValues();
//        loadings.storeValues();
//        factors.storeValues();
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogDetColKnown = logDetColKnown;
        storedLogDetCol = logDetCol;
        storedTrace = trace;
        storedTraceKnown = traceKnown;
        storedResidualKnown = residualKnown;
        storedLxFKnown = LxFKnown;
        storedFactorsKnown = factorsKnown;
        storedLoadingsKnown = loadingsKnown;
        storedDataKnown = dataKnown;
        storedTotalRecompute = totalRecompute;
        System.arraycopy(residual, 0, storedResidual, 0, residual.length);

        System.arraycopy(LxF, 0, storedLxF, 0, residual.length);
        System.arraycopy(changed, 0, storedChanged, 0, changed.length);
//        for (int i = 0; i <changedValues.size() ; i++) {
//            storedChangedValues.addElement(changedValues.elementAt(i));    ;
//        }
        storedChangedValues = (Vector<Integer>) changedValues.clone();

    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        changed = storedChanged;
//        data.restoreParameterValues();
//        loadings.restoreValues();
//        factors.restoreValues();
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        trace = storedTrace;
        traceKnown = storedTraceKnown;
        residualKnown = storedResidualKnown;
        LxFKnown = storedLxFKnown;
        double[] temp = residual;
        residual = storedResidual;
        storedResidual = temp;
        temp = LxF;
        LxF = storedLxF;
        storedLxF = temp;
        logDetCol = storedLogDetCol;
        logDetColKnown = storedLogDetColKnown;
        factorsKnown = storedFactorsKnown;
        loadingsKnown = storedLoadingsKnown;
        dataKnown = storedDataKnown;
        totalRecompute = storedTotalRecompute;
        changedValues = storedChangedValues;
//        changedValues=storedChangedValues;
//        storedChangedValues=new Vector<Integer>();

//        System.out.println(data.getParameterValue(10, 19));


//        int index=0;
//        for (int i = 0; i <continuous.getDimension() ; i++) {
//            if(continuous.getParameterValue(i)==0){
//                for (int j = 0; j <data.getParameter(i).getDimension() ; j++) {
//                    data.getParameter(i).setParameterValueQuietly(j, storedData.getParameter(index).getParameterValue(j));
//                }
//                index++;
//            }
//        }
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {
        // Do nothing
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == getScaledData()) {
            residualKnown = false;
            traceKnown = false;
            likelihoodKnown = false;
            if (!RecomputeResiduals) {
                if (index != -1 && !changedValues.contains(index))
                    changedValues.add(index);
                else{
                    totalRecompute = true;
                    changedValues.clear();
                }
                dataKnown = false;
            }
        }
        if (variable == factors) {


//            for (int i = 0; i <loadings.getRowDimension() ; i++) {
//                changed[i][index/factors.getRowDimension()]=true;
//            }

            if (!RecomputeFactors) {
//                System.out.println("index");
//                System.out.println(index);
                factorsKnown = false;
                int row = index / factors.getRowDimension();
                if (index != -1)
                    for (int i = 0; i < data.getRowDimension(); i++) {
                        if(!changedValues.contains(row * data.getRowDimension() + i))
                            changedValues.add(row * data.getRowDimension() + i);
                    }
                else{
                    totalRecompute = true;
                    changedValues.clear();
                }

            }


//            factorVariablesChanged.push(index);


            LxFKnown = false;
            residualKnown = false;
            traceKnown = false;
            likelihoodKnown = false;
        }
        if (variable == loadings) {
            if (!RecomputeLoadings) {
                loadingsKnown = false;
                int col = index % loadings.getRowDimension();
                if (index != -1) {
                    for (int i = 0; i < data.getColumnDimension(); i++) {
                        if(!changedValues.contains(i * data.getRowDimension() + col))
                            changedValues.add(i * data.getRowDimension() + col);
                    }
                }
                else{totalRecompute = true;
                    changedValues.clear();}
            }
//            System.out.println("Loadings Changed");
//            System.out.println(index);
//            System.out.println(index/loadings.getRowDimension());


//            for (int i = 0; i <factors.getColumnDimension(); i++) {
//                changed[index%loadings.getRowDimension()][i]=true;
//            }


//            factorVariablesChanged.push(index);


            LxFKnown = false;
            residualKnown = false;
            traceKnown = false;
            likelihoodKnown = false;
        }
        if (variable == colPrecision) {
            logDetColKnown = false;
            traceKnown = false;
            likelihoodKnown = false;
        }

    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Latent factor model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.CYBIS_2015_ASSESSING);
    }

    /**
     * Get the model.
     *
     * @return the model.
     */
    @Override
    public Model getModel() {
        return this;
    }

    /**
     * Get the log likelihood.
     *
     * @return the log likelihood.
     */
    @Override
    public double getLogLikelihood() {
        likelihoodKnown = false;
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    /**
     * Forces a complete recalculation of the likelihood next time getLikelihood is called
     */
    @Override
    public void makeDirty() {
        likelihoodKnown = false;

        factorsKnown = false;
        loadingsKnown = false;
        residualKnown = false;
        totalRecompute = true;
        changedValues.clear();
        LxFKnown = false;
        logDetColKnown = false;
        traceKnown = false;
        dataKnown = false;

    }

    private boolean checkLoadings() {
        for (int i = 0; i < StrictMath.min(loadings.getRowDimension(), loadings.getColumnDimension()); i++) {
            if (loadings.getParameterValue(i, i) < 0) {
                return false;
            }
        }
        return true;
    }

    private double calculateLogLikelihood() {
//         if(!checkLoadings()){
//             if(pathParameter==1)
//                return Double.NEGATIVE_INFINITY;
//            else{
//                return Math.log(1-pathParameter);}}

//        Matrix tRowPrecision= new Matrix(rowPrecision.getParameterAsMatrix());
//        Matrix tColPrecision= new Matrix(colPrecision.getParameterAsMatrix());


//        residualKnown=false;
        if (!residualKnown) {
            computeResiduals();
        }


//        expPart = residual.productInPlace(rowPrecision.productInPlace(residual.transposeThenProductInPlace(colPrecision, TResidualxC), RxTRxC), expPart);
//            logDetRow=StrictMath.log(rowPrecision.getDeterminant());
        //       logDetColKnown=false;
        if (!logDetColKnown) {
            logDetColKnown = true;
            logDetCol = 0;
            for (int i = 0; i < colPrecision.getRowDimension(); i++) {
                if (continuous.getParameterValue(i) != 0)
                    logDetCol += Math.log(colPrecision.getParameterValue(i, i)) * rowCount[i];
            }
        }
//            System.out.println(logDetCol);
//            System.out.println(logDetRow);
//        traceKnown=false;
        if (!traceKnown) {
            traceKnown = true;
            trace = TDTTrace(residual, colPrecision);
        }
//        if(expPart.getRowDimension()!=expPart.getColumnDimension())
//        {
//            System.err.print("Matrices are not conformable");
//            System.exit(0);
//        }

//        else{
//            for(int i=0; i<expPart.getRowDimension(); i++){
//                trace+=expPart.getParameterValue(i, i);
//            }
//        }
//        System.out.println(expPart);


        return -.5 * trace + .5 * logDetCol //+ .5 * data.getRowDimension()

                - .5 * nmeasurements * Math.log(2.0 * StrictMath.PI);
    }

//    public void setPathParameter(double beta){
//        pathParameter=beta;
//        data.product(pathParameter);
//    }

//    @Override
//    public double getLikelihoodCorrection() {
//        return 0;
//    }

    public double[] getLxF(){
        if(!LxFKnown)
            computeResiduals();
        return LxF;
    }
}
