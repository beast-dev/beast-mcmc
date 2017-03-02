/*
 * BitFlipOperator.java
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

import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * A generic operator that flips bits.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
public class BitFlipOperator extends SimpleMCMCOperator {

    public BitFlipOperator(Parameter parameter, double weight, boolean usesPriorOnSum, TreeModel treeModel) {

        this.parameter = parameter;
        this.usesPriorOnSum = usesPriorOnSum;

        if (treeModel != null) {
            indicators = new TreeParameterModel(treeModel, parameter, true);
            bitFlipHelper = new DriftBitFlipHelper();
            tree = treeModel;
        } else {
            bitFlipHelper = new BitFlipHelper();
        }

        setWeight(weight);
    }

    public BitFlipOperator(Parameter parameter, double weight, boolean usesPriorOnSum) {
        this(parameter, weight, usesPriorOnSum, null);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * Change the parameter and return the hastings ratio.
     * Flip (Switch a 0 to 1 or 1 to 0) for a random bit in a bit vector.
     * Return the hastings ratio which makes all subsets of vectors with the same number of 1 bits
     * equiprobable, unless usesPriorOnSum = false then all configurations are equiprobable
     */
    public final double doOperation() {
        final int dim = parameter.getDimension();
        //  final int dim = bitFlipHelper.getDim(parameter.getDimension());
        double sum = 0.0;
        // double sumNeg = 0.0;
        // double temp;

        if(usesPriorOnSum) {
            for (int i = 0; i < dim; i++) {
                //   if(parameter.getParameterValue(i)<0) {
                //     System.err.println(parameter.getParameterValue(i));
                // }
                sum += Math.abs(parameter.getParameterValue(i));
            }
// AR - removed a debugging printf
//            if (sum > 103) {
//                System.err.println("sum: " + sum);
//            }

            /*
            for (int i = 0; i < dim; i++) {
                temp = parameter.getParameterValue(i);
                if (temp < 0) {
                    sumNeg++;
                } else {
                   // sum += parameter.getParameterValue(i);
                    sum += temp;
                }
            }
            */
        }
        // make it so pos never corresponds to external nodes when used for drift diffusion model
        final int pos = MathUtils.nextInt(dim);

        final int value = (int) parameter.getParameterValue(pos);
        double logq = 0.0;
        if (value == 0) {
            logq = bitFlipHelper.flipZero(pos, dim, sum);
            //parameter.setParameterValue(pos, 1.0);
            //    if(sum==103){
            //        System.err.println("value: " + value);
            //       System.err.println("parameter.getParameterValue(pos): " + parameter.getParameterValue(pos));
            //    }

            //   if(usesPriorOnSum)
            //    logq = bitFlipHelper.getLogQFlipZero(dim, sum);
            // logq = -Math.log((dim - sum) / (sum + 1));

        } else if (value == 1) {
            // parameter.setParameterValue(pos, 0.0);
            logq = bitFlipHelper.flipOne(pos, dim, sum);
            //  if(usesPriorOnSum)
            //    logq = bitFlipHelper.getLogQFlipOne(dim, sum);
            // logq = -Math.log(sum / (dim - sum + 1));

        } else if (value == -1) {
            logq = bitFlipHelper.flipNegOne(pos, dim, sum);
            // logq = bitFlipHelper.getLogQFlipNegOne(dim,sum);
        } else {
            throw new RuntimeException("expected 1 or 0 or -1");
        }

        if (!usesPriorOnSum) {
            logq = 0;
        }

        // hastings ratio is designed to make move symmetric on sum of 1's
        return logq;
    }

    class BitFlipHelper {

        public BitFlipHelper() {
        }

        //  public int getDim (int paramDim){
        //      return paramDim;
        // }
        /*
        public double getLogQFlipZero(int dim, double sum, double sumNeg){
          //  System.err.println("I am NOT in DriftBitFlipHelper");
            return -Math.log((dim - sum) / (sum + 1));
        }

        public double getLogQFlipOne(int dim, double sum){
            return -Math.log(sum / (dim - sum + 1));
        }
        */
        public double flipOne(int pos, int dim, double sum) {
            parameter.setParameterValue(pos, 0.0);
            return -Math.log(sum / (dim - sum + 1));
        }

        public double flipZero(int pos, int dim, double sum) {
            parameter.setParameterValue(pos, 1.0);
            return -Math.log((dim - sum) / (sum + 1));
        }

        public double flipNegOne(int pos, int dim, double sum) {
            throw new RuntimeException("expected 1 or 0");
        }

    }

    class DriftBitFlipHelper extends BitFlipHelper {

        public DriftBitFlipHelper() {
        }

        //   public int getDim (int paramDim){
        //       return (int) (0.5*(paramDim+1)-1);
        //   }
        /*
        public double getLogQFlipZero(int dim, double sum){
         //   System.err.println("I am in DriftBitFlipHelper");
            if(isZeroPair){
                return Math.log((sum + 1)/(dim - 2*sum));
            }else{
                return 0;
            }
        }

        public double getLogQFlipOne(int dim, double sum){
            return Math.log((dim - 2*sum + 2)/sum);
        }

        public double getLogQFlipNegOne(int dim, double sum){
            return Math.log((dim - 2*sum + 2)/sum);
        }
        */

        public double flipOne(int pos, int dim, double sum) {
            // draw random number from [0,1]
            double rand = MathUtils.nextDouble();
            double internalDim = 0.5 * (dim + 1) - 1;

            if (rand < 0.5) {
                parameter.setParameterValue(pos, 0.0);
                return Math.log(2 * (internalDim + 1 - sum) / sum);
            } else {
                parameter.setParameterValue(pos, -1.0);
                return 0;
            }
        }

        public double flipNegOne(int pos, int dim, double sum) {
            double rand = MathUtils.nextDouble();
            double internalDim = 0.5 * (dim + 1) - 1;

            if (rand < 0.5) {
                parameter.setParameterValue(pos, 0.0);
                return Math.log(2 * (internalDim + 1 - sum) / sum);
            } else {
                parameter.setParameterValue(pos, 1.0);
                return 0;
            }
        }

        public double flipZero(int pos, int dim, double sum) {

            int nodeNum = indicators.getNodeNumberFromParameterIndex(pos);

            // if indicator corresponds to external node, keep default value of zero
            // would be unnecessary if pos never matches up with external nodes
            if (tree.isExternal(tree.getNode(nodeNum))) {
                //  if(nodeNum > 104) {
                //     System.err.println("external node num: " + nodeNum);
                // }
                parameter.setParameterValue(pos, 0.0);
                return 0;
            } else {

                double rand = MathUtils.nextDouble();
                double internalDim = 0.5 * (dim + 1) - 1;

                if (rand < 0.5) {
                    parameter.setParameterValue(pos, 1.0);
                    return Math.log((sum + 1) / (2 * (internalDim - sum)));
                } else {
                    parameter.setParameterValue(pos, -1.0);
                    return Math.log((sum + 1) / (2 * (internalDim - sum)));
                }
            }
        }
        /*
        public void flipZero(int pos){
            int nodeNumber = indicators.getNodeNumberFromParameterIndex(pos);
            NodeRef chosenNode = tree.getNode(nodeNumber);
            NodeRef parentNode = tree.getParent(chosenNode);
            NodeRef pairedNode = tree.getChild(parentNode, 0);
            int parentNumber = parentNode.getNumber();

            // check this
            if(pairedNode.equals(chosenNode)){
              //  System.err.println("I get HERE!!!!!");
                pairedNode = tree.getChild(parentNode, 1);
            }

                                    //    System.err.println("chosen value before: " + indicators.getNodeValue(tree,chosenNode));
                                    //   System.err.println("paired value after:" + indicators.getNodeValue(tree,pairedNode));


            if(indicators.getNodeValue(tree, pairedNode) == 0.0){
                // we have (0,0) pair
                isZeroPair = true;
                // flip (0,0) to (1,0)
                indicators.setNodeValue(tree, chosenNode, 1.0);
              //  parameter.setParameterValue(pos,1.0);
            }else if (indicators.getNodeValue(tree, pairedNode) == 1.0){
                // we have (0,1) pair
                isZeroPair = false;
                // flip (0,1) to (1,0)
                indicators.setNodeValue(tree, chosenNode, 1.0);
                indicators.setNodeValue(tree, pairedNode, 0.0);
             //   parameter.setParameterValue(pos,1.0);
             //   parameter.setParameterValue(indicators.getParameterIndexFromNodeNumber(pairedNode.getNumber()),0.0);
            }else {
                throw new RuntimeException("expected 1 or 0");
            }
        }

        protected boolean isZeroPair = false;
        */
    }


    // Interface MCMCOperator
    public final String getOperatorName() {
        return "bitFlip(" + parameter.getParameterName() + ")";
    }

    public final String getPerformanceSuggestion() {
        return "no performance suggestion";
    }

    public String toString() {
        return getOperatorName();
    }

    // Private instance variables

    private BitFlipHelper bitFlipHelper;
    private TreeModel tree;
    private TreeParameterModel indicators = null;
    private Parameter parameter = null;
//    private int bits;
    private boolean usesPriorOnSum = true;
}
