/*
 * LatentFactorModel.java
 *
 * Copyright (c) 2002-2014 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.math.matrixAlgebra.Matrix;
import dr.util.Citable;
import dr.util.Citation;

import java.util.List;


/**
 * @author Max Tolkoff
 * @author Marc Suchard
 */

public class LatentFactorModel extends AbstractModelLikelihood implements Citable {
//    private Matrix data;
//    private Matrix factors;
//    private Matrix loadings;

    private final MatrixParameter data;
    private final MatrixParameter factors;
    private final MatrixParameter loadings;
    private MatrixParameter sData;
    private final DiagonalMatrix rowPrecision;
    private final DiagonalMatrix colPrecision;
    private final Parameter continuous;

    private final boolean scaleData;

    private final int dimFactors;
    private final int dimData;
    private final int nTaxa;

    private boolean firstTime=true;

    private boolean likelihoodKnown = false;
    private boolean isDataScaled=false;
    private boolean storedLikelihoodKnown;
    private boolean residualKnown=false;
    private boolean LxFKnown=false;
    private boolean storedResidualKnown=false;
    private boolean storedLxFKnown;
    private boolean traceKnown=false;
    private boolean storedTraceKnown;
    private boolean logDetColKnown=false;
    private boolean storedLogDetColKnown;
    private double trace;
    private double storedTrace;
    private double logLikelihood;
    private double storedLogLikelihood;
    private double logDetCol;
    private double storedLogDetCol;
    private boolean[][] changed;
    private boolean[][] storedChanged;

    private double[] residual;
    private double[] LxF;
    private double[] storedResidual;
    private double[] storedLxF;

    public LatentFactorModel(MatrixParameter data, MatrixParameter factors, MatrixParameter loadings,
                             DiagonalMatrix rowPrecision, DiagonalMatrix colPrecision,
                             boolean scaleData, Parameter continuous
    ) {
        super("");
//        data = new Matrix(dataIn.getParameterAsMatrix());
//        factors = new Matrix(factorsIn.getParameterAsMatrix());
//        loadings = new Matrix(loadingsIn.getParameterAsMatrix());
        this.scaleData=scaleData;
        this.data = data;
        this.factors = factors;
        // Put default bounds on factors
        for (int i = 0; i < factors.getParameterCount(); ++i) {
            Parameter p = factors.getParameter(i);
            System.err.println(p.getId() + " " + p.getDimension());
            p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, p.getDimension()));
        }
        this.continuous=continuous;

        this.loadings = loadings;

//        storedData=new MatrixParameter(null);
//        for (int i = 0; i <continuous.getDimension(); i++) {
//            if(continuous.getParameterValue(i)==0)
//                storedData.addParameter(new Parameter.Default(data.getColumnDimension()));
//
//        }

        // Put default bounds on loadings
//        loadings.addBounds();


        changed=new boolean[loadings.getRowDimension()][factors.getColumnDimension()];
        storedChanged=new boolean[loadings.getRowDimension()][factors.getColumnDimension()];

        for (int i = 0; i <loadings.getRowDimension() ; i++) {
            for (int j = 0; j <factors.getColumnDimension() ; j++) {
                changed[i][j]=false;
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


        if (nTaxa * dimData != data.getDimension()) {
            throw new RuntimeException("LOADINGS MATRIX AND FACTOR MATRIX MUST HAVE EXTERNAL DIMENSIONS WHOSE PRODUCT IS EQUAL TO THE NUMBER OF DATA POINTS\n");
//            System.exit(10);
        }
        if (dimData < dimFactors) {
            throw new RuntimeException("MUST HAVE FEWER FACTORS THAN DATA POINTS\n");
        }

        residual=new double[loadings.getRowDimension()*factors.getColumnDimension()];
        LxF=new double[loadings.getRowDimension()*factors.getColumnDimension()];
        storedResidual=new double[residual.length];
        storedLxF=new double[LxF.length];

        if(!isDataScaled & !scaleData){
            sData=this.data;
            isDataScaled=true;
        }
        if(!isDataScaled){
            sData = computeScaledData();
            isDataScaled=true;
            for (int i = 0; i <sData.getRowDimension() ; i++) {
                for (int j = 0; j <sData.getColumnDimension() ; j++) {
                        this.data.setParameterValue(i,j,sData.getParameterValue(i,j));
                    System.out.println(this.data.getParameterValue(i,j));
                }

            }
            data.fireParameterChangedEvent();
        }

//       computeResiduals();
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
    public MatrixParameter getFactors(){return factors;}

    public MatrixParameter getColumnPrecision(){return colPrecision;}

    public MatrixParameter getLoadings(){return loadings;}

    public MatrixParameter getData(){return data;}

    public Parameter returnIntermediate(){
        if(!residualKnown && checkLoadings()){
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


    public MatrixParameter getScaledData(){return data;}

    public Parameter getContinuous(){return continuous;}

    public int getFactorDimension(){return factors.getRowDimension();}

    private void Multiply(MatrixParameter Left, MatrixParameter Right, double[] answer){
        int dim=Left.getColumnDimension();
        int n=Left.getRowDimension();
        int p=Right.getColumnDimension();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < p; j++) {
                double sum = 0;
                for (int k = 0; k < dim; k++)
                    sum += Left.getParameterValue(i, k) * Right.getParameterValue(k,j);
                answer[i*p+j]=sum;
            }
        }
    }

    private void add(MatrixParameter Left, MatrixParameter Right, double[] answer){
        int row=Left.getRowDimension();
        int col=Left.getColumnDimension();
        for (int i = 0; i <row ; i++) {
            for (int j = 0; j < col; j++) {
                answer[i*col+j]=Left.getParameterValue(i,j)+Right.getParameterValue(i,j);
            }

        }
    }

    private void subtract(MatrixParameter Left, double[] Right, double[] answer){
        int row=Left.getRowDimension();
        int col=Left.getColumnDimension();
        boolean containsDiscrete=false;
        for (int i = 0; i <row ; i++) {
            if(continuous.getParameterValue(i)!=0){
                for (int j = 0; j < col; j++) {
                    answer[i*col+j]=Left.getParameterValue(i,j)-Right[i*col+j];
                }
            }
//            else{
//                for (int j = 0; j <col; j++) {
//                    Left.setParameterValueQuietly(i,j, Right[i*col+j]);
//                }
//                containsDiscrete=true;
//            }

        }
//        if(containsDiscrete){
//            Left.fireParameterChangedEvent();}
    }

    private double TDTTrace(double[] array, DiagonalMatrix middle){
        int innerDim=middle.getRowDimension();
        int outerDim=array.length/innerDim;
        double sum=0;
        for (int j = 0; j <innerDim ; j++){
            if(continuous.getParameterValue(j)!=0) {
                for (int i = 0; i < outerDim; i++) {
                        double s1 = array[j * outerDim + i];
                        double s2 = middle.getParameterValue(j, j);
                        sum += s1 * s1 * s2;
                }
            }
        }
        return sum;
    }



    private MatrixParameter computeScaledData(){
        MatrixParameter answer=new MatrixParameter(data.getParameterName() + ".scaled");
        answer.setDimensions(data.getRowDimension(), data.getColumnDimension());
 //       Matrix answer=new Matrix(data.getRowDimension(), data.getColumnDimension());
        double[][] aData=data.getParameterAsMatrix();
        double[] meanList=new double[data.getRowDimension()];
        double[] varList=new double[data.getRowDimension()];
        for(int i=0; i<data.getColumnDimension(); i++){
            for (int j=0; j<data.getRowDimension(); j++){
                meanList[j]+=data.getParameterValue(j,i);
            }
        }
        for(int i=0; i<data.getRowDimension(); i++){
            if(continuous.getParameterValue(i)==1)
                meanList[i]=meanList[i]/data.getColumnDimension();
            else
                meanList[i]=0;
        }

        double[][] answerTemp=new double[data.getRowDimension()][data.getColumnDimension()];
        for(int i=0; i<data.getColumnDimension(); i++){
            for(int j=0; j<data.getRowDimension(); j++){
                answerTemp[j][i]=aData[j][i]-meanList[j];
            }
        }
//        System.out.println(new Matrix(answerTemp));

        for(int i=0; i<data.getColumnDimension(); i++){
            for(int j=0; j<data.getRowDimension(); j++){
                varList[j]+=answerTemp[j][i]*answerTemp[j][i];
            }
        }

        for(int i=0; i<data.getRowDimension(); i++){
            if(continuous.getParameterValue(i)==1){
            varList[i]=varList[i]/(data.getColumnDimension()-1);
            varList[i]=StrictMath.sqrt(varList[i]);}
            else{
                varList[i]=1;
            }
        }
//        System.out.println(data.getColumnDimension());
//        System.out.println(data.getRowDimension());

        for(int i=0; i<data.getColumnDimension(); i++){
            for(int j=0; j<data.getRowDimension(); j++){
                answer.setParameterValue(j,i, answerTemp[j][i]/varList[j]);
            }
        }
//        System.out.println(new Matrix(answerTemp));
        computeResiduals();
        return answer;
    }

    private Matrix copy(CompoundParameter parameter, int dimMajor, int dimMinor) {
        return new Matrix(parameter.getParameterValues(), dimMajor, dimMinor);
    }

    private void computeResiduals() {
//    LxFKnown=false;



//        if(firstTime || (!factorVariablesChanged.empty() && !loadingVariablesChanged.empty())){
    if(!LxFKnown){
    Multiply(loadings, factors, LxF);
        LxFKnown=true;
    }
        subtract(data, LxF, residual);
        residualKnown=true;
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
        data.storeParameterValues();
        loadings.storeValues();
        factors.storeValues();
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogDetColKnown=logDetColKnown;
        storedLogDetCol=logDetCol;
        storedTrace=trace;
        storedTraceKnown=traceKnown;
        storedResidualKnown=residualKnown;
        storedLxFKnown=LxFKnown;
        System.arraycopy(residual, 0, storedResidual, 0, residual.length);

//        System.arraycopy(LxF, 0, storedLxF, 0, residual.length);
//        System.arraycopy(changed, 0, storedChanged, 0, residual.length);

//        System.out.println(data.getParameterValue(10, 19));


//        int index=0;
//        for (int i = 0; i <continuous.getDimension() ; i++) {
//            if(continuous.getParameterValue(i)==0){
//                for (int j = 0; j <data.getParameter(i).getDimension() ; j++) {
//                    storedData.getParameter(index).setParameterValueQuietly(j, getData().getParameterValue(i,j));
//                }
//                index++;
//            }
//        }
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        changed=storedChanged;
        data.restoreParameterValues();
        loadings.restoreValues();
        factors.restoreValues();
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        trace=storedTrace;
        traceKnown=storedTraceKnown;
        residualKnown=storedResidualKnown;
        LxFKnown=storedLxFKnown;
        residual=storedResidual;
        storedResidual=new double[residual.length];
        LxF=storedLxF;
        storedLxF=new double[LxF.length];
        logDetCol=storedLogDetCol;
        logDetColKnown=storedLogDetColKnown;

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
//        if(variable==data){
//            isDataScaled=false;
//            residualKnown=false;
//            traceKnown=false;
//        }
//        if(variable==loadings){
//        System.out.println(variable.getVariableName());
//        System.out.println(index);
//            System.out.println(loadings.getParameterValue(index));}
        if(variable==factors){



//            int column=index/factors.getColumnDimension();
//            for (int i = 0; i <factors.getColumnDimension() ; i++) {
//                changed[i][column]=true;
//
//            }




//            factorVariablesChanged.push(index);


            LxFKnown=false;
            residualKnown=false;
            traceKnown=false;
            likelihoodKnown = false;
        }
        if(variable==loadings){


//            factorVariablesChanged.push(index);


            LxFKnown=false;
            residualKnown=false;
            traceKnown=false;
            likelihoodKnown = false;
        }
        if(variable==colPrecision){
            logDetColKnown=false;
            traceKnown=false;
            likelihoodKnown = false;
        }

    }

    /**
     * @return a list of citations associated with this object
     */
    @Override
    public List<Citation> getCitations() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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
       likelihoodKnown=false;
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
    }

    private boolean checkLoadings(){
        for(int i=0; i<StrictMath.min(loadings.getRowDimension(),loadings.getColumnDimension()); i++)
        {
            if(loadings.getParameter(i).getParameterValue(i)<0)
            {
                return false;
            }
        }
        return true;
    }

    private double calculateLogLikelihood() {
         if(!checkLoadings())
            return Double.NEGATIVE_INFINITY;

//        Matrix tRowPrecision= new Matrix(rowPrecision.getParameterAsMatrix());
//        Matrix tColPrecision= new Matrix(colPrecision.getParameterAsMatrix());


//        residualKnown=false;
        if(!residualKnown){
            computeResiduals();
        }



//        expPart = residual.productInPlace(rowPrecision.productInPlace(residual.transposeThenProductInPlace(colPrecision, TResidualxC), RxTRxC), expPart);
//            logDetRow=StrictMath.log(rowPrecision.getDeterminant());
 //       logDetColKnown=false;
        if(!logDetColKnown){
            logDetColKnown=true;
            double product=1;
            for (int i = 0; i <colPrecision.getRowDimension() ; i++) {
                if (continuous.getParameterValue(i)!=0)
                    product*=colPrecision.getParameterValue(i,i);
            }

            logDetCol=StrictMath.log(product);
        }
//            System.out.println(logDetCol);
//            System.out.println(logDetRow);
//        traceKnown=false;
        if(!traceKnown){
            traceKnown=true;
            trace=TDTTrace(residual, colPrecision);
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


       return -.5*trace + .5*data.getColumnDimension()*logDetCol

               -.5*data.getRowDimension()*data.getColumnDimension()*Math.log(2.0 * StrictMath.PI);
    }
}
