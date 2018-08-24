/*
 * MultivariateNormalGibbsOperator.java
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

package dr.inference.operators;

import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.distribution.MultivariateNormalDistributionModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import dr.util.Attribute;

import java.util.List;

/**
@author Max Tolkoff
 */
public class MultivariateNormalGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private Matrix priorPrecision;
    private Vector priorMean;
    private MatrixParameter likelihoodPrecision;
    private Parameter likelihoodMean;
    private MultivariateDistributionLikelihood likelihood;
    private int dim;
    public static final String MVN_GIBBS="multivariateNormalGibbsOperator";


    public MultivariateNormalGibbsOperator(MultivariateDistributionLikelihood likelihood, MultivariateDistributionLikelihood prior, Double weight) throws IllegalDimension {

        
        MultivariateNormalDistribution tempPrior=(MultivariateNormalDistribution) prior.getDistribution();
        this.priorMean=new Vector(tempPrior.getMean());
        this.priorPrecision=new Matrix(tempPrior.getScaleMatrix());
         MultivariateNormalDistributionModel tempLikelihood=(MultivariateNormalDistributionModel) likelihood.getDistribution();
        this.likelihoodMean=tempLikelihood.getMeanParameter();
        this.likelihoodPrecision=tempLikelihood.getPrecisionMatrixParameter();
        this.likelihood=likelihood;
        this.dim=likelihoodMean.getValues().length;
//        if(dataTemp.contains(MatrixParameter.class))
//        {System.err.print("Well, at least you know it's there...\n");}
//        else{System.err.print("Nope, you screwed up\n");}

        setWeight(weight);
    }


    private void setParameterValue(Parameter set, double[] value){
        set.setDimension(value.length);
        for(int i=0; i<value.length; i++)
        {set.setParameterValueQuietly(i,value[i]);}
        set.fireParameterChangedEvent();
    }

    private double[] getMeanSum(){
        double[] answer=new double[dim];
        List<Attribute<double[]>> dataList = likelihood.getDataList();
        for(Attribute<double[]> d: dataList){
            for(int i=0; i<d.getAttributeValue().length; i++)
            {
                answer[i]+=d.getAttributeValue()[i];
            }
        }
/*
for(int i=0; i<dim; i++){
System.err.print(answer[i]);
System.err.print("\n");}
*/
        return answer;}

    private Matrix getPrecision() throws IllegalDimension {
        Matrix currentPrecision=new Matrix(likelihoodPrecision.getParameterAsMatrix());
        currentPrecision=currentPrecision.product(likelihood.getDataList().size());

/*
for(int i=0; i<currentPrecision.columns(); i++){
for(int j=0; j<currentPrecision.rows(); j++){
System.err.print(currentPrecision.toComponents()[i][j]);
System.err.print(" ");}
System.err.print("\n"); }
*/
        return priorPrecision.add(currentPrecision);
    }

    private Vector getMean() throws IllegalDimension {
        Vector meanSum=new Vector(getMeanSum());
        Matrix workingPrecision=new Matrix(likelihoodPrecision.getParameterAsMatrix());
        Vector meanPart=workingPrecision.product(meanSum);
        meanPart=meanPart.add(priorPrecision.product(priorMean));
        Matrix varPart=getPrecision().inverse();
        Vector answer=varPart.product(meanPart);
/*
for(int i=0; i<varPart.columns(); i++){
for(int j=0; j<varPart.rows(); j++){
System.err.print(varPart.toComponents()[i][j]);}
System.err.print("\n"); }
this.priorPrecision=new Matrix(prior.getDistribution().getScaleMatrix());
System.err.print(answer.toComponents()[0]);
System.err.print("\n");
System.out.print(answer.toComponents()[0]);
for(int i=0; i<answer.dimension(); i++){
System.err.print(answer.toComponents()[i]);}
System.err.print("\n");
*/

        return answer;
    }

//    private Vector getDraws() throws IllegalDimension{
//        double[] rUniform=new double[dim];
//        for(int i=0; i<dim; i++)
//        {rUniform[i]=}
//        Vector draws=new Vector(MultivariateNormalDistribution.);
//        return draws;
//    }

    @Override
    public String getOperatorName() {
        return MVN_GIBBS;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public double doOperation() {
        double[] draws=null;


//        try {
//            for(int i=0; i<getPrecision().columns(); i++){
//            for(int j=0; j<getPrecision().rows(); j++){
//            System.err.print(getPrecision().toComponents()[i][j]);
//                System.err.print(", ");}}
////            System.err.print(" ");}
////            System.err.print("\n"); }
//        } catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }



        try {
            draws=MultivariateNormalDistribution.nextMultivariateNormalPrecision(getMean().toComponents(), getPrecision().toComponents());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }

//        for(int i=0; i<dim; i++){
//            System.err.print(draws[i]);
//            System.err.print("\n");}
        setParameterValue(likelihoodMean, draws);

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

}
