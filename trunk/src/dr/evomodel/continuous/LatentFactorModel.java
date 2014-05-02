/*
 * LatentFactorModel.java
 *
 * Copyright (c) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.inference.model.*;
import dr.math.matrixAlgebra.IllegalDimension;
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

    private Matrix residual;

    private final MatrixParameter data;
    private final MatrixParameter factors;
    private final MatrixParameter loadings;
    private final DiagonalMatrix rowPrecision;
    private final DiagonalMatrix colPrecision;

    private final int dimFactors;
    private final int dimData;
    private final int nTaxa;

    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;

    public LatentFactorModel(MatrixParameter data, MatrixParameter factors, MatrixParameter loadings,
                             DiagonalMatrix rowPrecision, DiagonalMatrix colPrecision
    ) {
        super("");
//        data = new Matrix(dataIn.getParameterAsMatrix());
//        factors = new Matrix(factorsIn.getParameterAsMatrix());
//        loadings = new Matrix(loadingsIn.getParameterAsMatrix());
        this.data = data;
        this.factors = factors;
        // Put default bounds on factors
        for (int i = 0; i < factors.getParameterCount(); ++i) {
            Parameter p = factors.getParameter(i);
            System.err.println(p.getId() + " " + p.getDimension());
            p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, p.getDimension()));
        }


        this.loadings = loadings;

        // Put default bounds on loadings
//        loadings.addBounds();


        this.rowPrecision = rowPrecision;
        this.colPrecision = colPrecision;

        addVariable(data);
        addVariable(factors);
        addVariable(loadings);
        addVariable(rowPrecision);
        addVariable(colPrecision);

        dimFactors = factors.getRowDimension();
        dimData = loadings.getColumnDimension();
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

    private Matrix copy(CompoundParameter parameter, int dimMajor, int dimMinor) {
        return new Matrix(parameter.getParameterValues(), dimMajor, dimMinor);
    }

    private Matrix computeResiduals() {
//        Parameter[] dataTemp=new Parameter[nTaxa];
//        for(int i=0; i<nTaxa; i++)
//        {
//            dataTemp[i] = new Parameter.Default(dimData);
//            for(int j=0; j<dimData; j++)
//            {
//                dataTemp[i].setParameterValue(j, data.getParameterValue(i*dimData+j));
//            }
//
//        }
//        MatrixParameter dataMatrix=new MatrixParameter(null, dataTemp);

        Matrix residual = null;
        Matrix tLoadings = new Matrix(loadings.getParameterAsMatrix());
        Matrix tData = new Matrix(data.getParameterAsMatrix());
        Matrix tFactors = new Matrix(factors.getParameterAsMatrix());
        try {
            residual = tData.subtract(tLoadings.transpose().product(tFactors));
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        return residual;
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
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood; // TODO Possible error in store/restore -- when changed to 42, no error arises
        likelihoodKnown = storedLikelihoodKnown;
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
        likelihoodKnown = false;
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

    private double calculateLogLikelihood() {
        for(int i=0; i<StrictMath.min(loadings.getRowDimension(),loadings.getColumnDimension()); i++)
        {
            if(loadings.getParameter(i).getParameterValue(i)<0)
            {
                return Double.NEGATIVE_INFINITY;
            }
        }
        Matrix tRowPrecision= new Matrix(rowPrecision.getParameterAsMatrix());
        Matrix tColPrecision= new Matrix(colPrecision.getParameterAsMatrix());
        Matrix residual = computeResiduals();
//        computeResiduals();
        Matrix expPart=null;
        double logDetRow=0;
        double logDetCol=0;
        try{
        expPart = residual.product(tRowPrecision.product(residual.transpose())).product(tColPrecision);
            logDetRow=tRowPrecision.logDeterminant();
            logDetCol=tColPrecision.logDeterminant();
        } catch (IllegalDimension illegalDimension) {
        illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        double trace=0;
        if(!expPart.isSquare())
        {
            System.err.print("Matrices are not conformable");
            System.exit(0);
        }
        else{
            for(int i=0; i<expPart.rows(); i++){
                trace+=expPart.component(i,i);
            }
        }
       return -.5*trace - .5*tColPrecision.rows()*logDetCol-.5*tRowPrecision.rows()*logDetRow-.5*tRowPrecision.rows()*tColPrecision.rows()*StrictMath.log(StrictMath.PI);
    }
}
