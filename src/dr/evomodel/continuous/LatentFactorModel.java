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
    private final CompoundParameter factors;
    private final MatrixParameter loadings;
    private final MatrixParameter precision;

    private final int dimFactors;
    private final int dimData;
    private final int nTaxa;

    private boolean likelihoodKnown = false;
    private double logLikelihood;

    public LatentFactorModel(MatrixParameter data, CompoundParameter factors, MatrixParameter loadings,
                             MatrixParameter precision,
                             int numFactors) {
        super("");
//        data = new Matrix(dataIn.getParameterAsMatrix());
//        factors = new Matrix(factorsIn.getParameterAsMatrix());
//        loadings = new Matrix(loadingsIn.getParameterAsMatrix());
        this.data = data;
        this.factors = factors;
        this.loadings = loadings;
        this.precision = precision;

        addVariable(data);
        addVariable(factors);
        addVariable(loadings);

        dimFactors = numFactors;
        nTaxa = data.getColumnDimension();
        dimData = data.getRowDimension();

        // TODO Check dimensions of loadings (dimFactors x dimData)
        // TODO dimData >= dimFactors

        computeResiduals();
        System.out.print(new Matrix(residual.toComponents()));
        System.out.print(calculateLogLikelihood());
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

    private void computeResiduals() {

        Matrix tLoadings = new Matrix(loadings.getParameterAsMatrix());
        Matrix tData = new Matrix(data.getParameterAsMatrix());
        Matrix tFactors = copy(factors, nTaxa, dimFactors);
        try {
            residual = tData.subtract(tFactors.product(tLoadings).transpose());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }



    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    @Override
    protected void storeState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void restoreState() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * This call specifies that the current state is accept. Most models will not need to do anything.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    @Override
    protected void acceptState() {
        //To change body of implemented methods use File | Settings | File Templates.
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
        return null;  //To change body of implemented methods use File | Settings | File Templates.
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
        //To change body of implemented methods use File | Settings | File Templates.
    }

    private double calculateLogLikelihood() {
        //TODO return correct answer
        Matrix tPrecision= new Matrix(precision.getParameterAsMatrix());
        computeResiduals();
        Matrix expPart=null;
        double logDet=0;
        try{
        expPart = residual.product(tPrecision.product(residual.transpose()));
            logDet=tPrecision.logDeterminant();
        } catch (IllegalDimension illegalDimension) {
        illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
       return -.5*expPart.toComponents()[0][0] - .5*logDet;
    }
}
