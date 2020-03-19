/*
 * MomentDistributionModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;


import dr.inference.model.*;
import dr.inference.model.Parameter;
import dr.inferencexml.distribution.MomentDistributionModelParser;
import dr.math.distributions.RandomGenerator;

//@author Max Tolkoff
public class MomentDistributionModel extends AbstractModelLikelihood implements ParametricMultivariateDistributionModel, RandomGenerator {

    public MomentDistributionModel(String id, Parameter mean, Parameter precision, Parameter cutoff, Parameter data) {
        super(id);

        this.mean=mean;
        this.precision=precision;
//        this.mean = new DuplicatedParameter(mean);
//        this.mean.addDuplicationParameter(new Parameter.Default(cutoff.getDimension()));
//        DuplicatedParameter precTemp= new DuplicatedParameter(precision);
//        precTemp.addDuplicationParameter(new Parameter.Default(cutoff.getDimension()));
//        this.precision=new DiagonalMatrix(precTemp);
        addVariable(mean);
        mean.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, 1));
        addVariable(precision);
//        precision.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        this.cutoff=cutoff;
        if(cutoff!=null){
            addVariable(cutoff);
            int boundsDim;
            if(cutoff.getBounds()!=null)
                boundsDim = cutoff.getBounds().getBoundsDimension();
            else boundsDim = cutoff.getDimension();
            cutoff.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, boundsDim));}
        addVariable(data);
        this.data=data;
        untruncated=new NormalDistributionModel(mean, precision, true);
        sumKnown=false;
        untruncatedKnown=false;
    }

    private final Parameter mean;
    private final Parameter precision;
//    private final DuplicatedParameter mean;
//    private final DiagonalMatrix precision;
    private final Parameter cutoff;
    private NormalDistributionModel untruncated;
    private double sum;
    private boolean sumKnown;

    private boolean storedSumKnown;
    private double storedSum;
    private boolean untruncatedKnown;
    private boolean storedUntruncatedKnown;
    private NormalDistributionModel storedUntruncated;
    private Parameter data;

    public double logPdf(Parameter data) {
//        untruncatedKnown=false;
//        sumKnown=false;
        checkDistribution();
        if(sumKnown)
            return sum;
        else
        {
            sum=0;
        }

        if(cutoff!=null){
            if(data.getDimension()!=cutoff.getDimension()){
                throw new RuntimeException("Incorrect number of cutoffs");
            }
            for (int i = 0; i <data.getDimension() ; i++) {
                if (Math.sqrt(cutoff.getParameterValue(i)) - .05 > Math.abs(data.getParameterValue(i)) && data.getParameterValue(i)!=0){
                        return Double.NEGATIVE_INFINITY;                                                                          }
                    else if(data.getParameterValue(i)==0){
                        sum += 0;
                }
//                        sum+=-1000-Math.log(precision.getParameterValue(0));
                    else {
                        sum+=untruncated.logPdf(data.getParameterValue(i))
                                + getNormalizingConstant(i)
                                ;
                        }
            }
        }
        else{
            for (int i = 0; i <data.getDimension() ; i++) {
                sum+= untruncated.logPdf(data.getParameterValue(i)) + 2 * StrictMath.log(Math.abs(data.getParameterValue(i))) + StrictMath.log(precision.getParameterValue(0));
            }
        }
        sumKnown=true;
        return sum;

    }

    public double getNormalizingConstant(int i){
        return -Math.log(1 - (untruncated.cdf(Math.sqrt(cutoff.getParameterValue(i))) - untruncated.cdf(-Math.sqrt(cutoff.getParameterValue(i)))));
    }

    @Override
    public double logPdf(double[] x) {
        return 0;
    }

    public Parameter getCutoff(){return cutoff;}

    @Override
    public double[][] getScaleMatrix() {
        double[][] temp=new double[1][1];
        temp[0][0]=precision.getParameterValue(0);
        return temp;
//        return precision.getParameterAsMatrix();
    }

    @Override
    public double[] getMean() {
        return mean.getParameterValues();
    }

    @Override
    public String getType() {
        return "Moment Distribution Model";
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        sumKnown=false;
        if(variable==mean || variable==precision)
        {untruncatedKnown=false;}

    }

    @Override
    protected void storeState() {
        storedSumKnown=sumKnown;
        storedSum=sum;
        storedUntruncated=untruncated;
        storedUntruncatedKnown=untruncatedKnown;
    }

    @Override
    protected void restoreState() {
        sumKnown=storedSumKnown;
        sum=storedSum;
        untruncated=storedUntruncated;
        untruncatedKnown=storedUntruncatedKnown;

    }

    @Override
    protected void acceptState() {

    }

    private NormalDistributionModel createNewDistribution() {
        return new NormalDistributionModel(mean, precision, true);
//        return new NormalDistributionModel(new Parameter.Default(mean.getParameterValue(0)), new Parameter.Default(precision.getParameterValue(0)), true);
    }

    private void checkDistribution() {
        if (!untruncatedKnown) {
            untruncated = createNewDistribution();
            untruncatedKnown = true;
        }
    }

    @Override
    public double[] nextRandom() {
        return new double[0];
    }

    @Override
    public double logPdf(Object x) {
        if(x instanceof Parameter)
         return logPdf((Parameter) x);
        else
            return 0;
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        return logPdf(data);
    }

    @Override
    public void makeDirty() {
        sumKnown=false;
        untruncatedKnown=false;

    }

    // *****************************************************************
    // Interface DensityModel
    // *****************************************************************

    @Override
    public Variable<Double> getLocationVariable() {
        throw new UnsupportedOperationException("Not implemented");
    }

}
