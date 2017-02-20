/*
 * ARGPartitioningOperator.java
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

package dr.evomodel.arg.operators;

import dr.evomodel.arg.ARGModel;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.logging.Logger;

public class ARGPartitioningOperator extends SimpleMCMCOperator {

    private final CompoundParameter partitioningParameters;
    private final ARGModel arg;

    public final static String OPERATOR_NAME = "argPartitionOperator";
    public static final String TOSS_SIZE = "tossSize";
    public static final String TOSS_ALL = "tossAll";

    private final boolean tossAll;
    private final boolean isRecombination;
    private final int tossSize;

    public ARGPartitioningOperator(ARGModel arg, int tossSize, int weight, boolean tossAll) {
        super.setWeight(weight);

        this.arg = arg;
        this.partitioningParameters = arg.getPartitioningParameters();
        this.tossSize = tossSize;
        this.isRecombination = arg.isRecombinationPartitionType();
        
        this.tossAll = tossAll;
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return partitioningParameters;
    }


    public final double doOperation() {
        double logq = 0;

        final int len = partitioningParameters.getParameterCount();

        if (len == 0) {
            return 0;
        }

        boolean[] updatePartition = new boolean[arg.getNumberOfPartitions()];

        if(tossAll){
            for(int i = 0 ; i < len; i++){
        		logq += doFlip(i,updatePartition);
            }	
        }else{
        	logq = doFlip(MathUtils.nextInt(len),updatePartition);
        }

        arg.fireModelChanged(new PartitionChangedEvent(partitioningParameters, updatePartition));
        return logq;
    }
    
    private double doFlip(int i, boolean[] updatePartition) {
    	if (isRecombination) {
            return doRecombination(partitioningParameters.getParameter(i),updatePartition);
        } 
        
    	return doReassortment(partitioningParameters.getParameter(i),updatePartition);
        
    }


    private double doRecombination(Parameter partition, boolean[] updatePartition) {

        assert checkValidRecombinationPartition(partition);

        int currentBreakLocation = 0;
        for (int i = 0, n = arg.getNumberOfPartitions(); i < n; i++) {
            if (partition.getParameterValue(i) == 1) {
                currentBreakLocation = i;
                break;
            }
        }

        assert currentBreakLocation > 0;

        if (MathUtils.nextBoolean()) {
            //Move break right 1
            partition.setParameterValueQuietly(currentBreakLocation, 0.0);
            updatePartition[currentBreakLocation] = true;
        } else {
            partition.setParameterValueQuietly(currentBreakLocation - 1, 1.0);
            updatePartition[currentBreakLocation-1] = true;
        }

        if (!checkValidRecombinationPartition(partition)) {
            return Double.NEGATIVE_INFINITY;
        }


        return 0;
    }

    public static boolean checkValidRecombinationPartition(Parameter partition) {
        int l = partition.getDimension();
//        if ((partition.getParameterValue(0) == 0 && partition.getParameterValue(l - 1) == 1))
//            return true;
//
//        return false;
        return (partition.getParameterValue(0) == 0 && partition.getParameterValue(l - 1) == 1);
    }


    private double doReassortment(Parameter partition, boolean[] updatePartition) {

        assert checkValidReassortmentPartition(partition);

        ArrayList<Integer> list = new ArrayList<Integer>(tossSize);

        while (list.size() < tossSize) {
            int a = MathUtils.nextInt(arg.getNumberOfPartitions() - 1) + 1;
            if (!list.contains(a)) {
                list.add(a);
            }
        }


        for (int a : list) {
            if (partition.getParameterValue(a) == 0) {
                partition.setParameterValueQuietly(a, 1);
            } else {
                partition.setParameterValueQuietly(a, 0);
            }
            updatePartition[a] = true;            
        }
        
        
        if (!checkValidReassortmentPartition(partition)) {
            return Double.NEGATIVE_INFINITY;
        }

        return 0;
    }

    public static boolean checkValidReassortmentPartition(Parameter partition) {
        if (partition.getParameterValue(0) != 0)
            return false;

        double[] a = partition.getParameterValues();

        double sum = 0;

        for (double b : a)
            sum += b;

//        if (sum == 0 || sum == a.length)
//            return false;
//
//        return true;
        return !(sum == 0 || sum == a.length);

    }

    @Override
    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    public class PartitionChangedEvent {
        Parameter partitioning;
        boolean[] updatePartition;

        public PartitionChangedEvent(Parameter partitioning, boolean[] updatePartition) {
            this.partitioning = partitioning;
            this.updatePartition = updatePartition;
        }

        public Parameter getParameter() {
            return partitioning;
        }

        public boolean[] getUpdatedPartitions() {
            return updatePartition;
        }
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public String[] getParserNames() {
            return new String[]{
                    OPERATOR_NAME,
                    "tossPartitioningOperator",
            };
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int weight = xo.getIntegerAttribute(WEIGHT);

            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            int tossSize = 1;
            if (xo.hasAttribute(TOSS_SIZE)) {
                tossSize = xo.getIntegerAttribute(TOSS_SIZE);

                if (tossSize <= 0 || tossSize >= arg.getNumberOfPartitions()) {
                    throw new XMLParseException("Toss size is incorrect");
                }
            }
            
            boolean tossAll = false;
            if(xo.hasAttribute(TOSS_ALL)){
            	tossAll = xo.getBooleanAttribute(TOSS_ALL);
            }

            Logger.getLogger("dr.evomodel").info("Creating ARGPartitionOperator: " + TOSS_SIZE + "=" + tossSize +
            		" " + TOSS_ALL + "=" + tossAll);


            return new ARGPartitioningOperator(arg, tossSize, weight, tossAll);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that picks a new partitioning uniformly at random.";
        }

        public Class getReturnType() {
            return ARGPartitioningOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule(WEIGHT),
                AttributeRule.newIntegerRule(TOSS_SIZE,true),
                AttributeRule.newBooleanRule(TOSS_ALL,true),
                new ElementRule(ARGModel.class)
        };
    };

    public String toString() {
        return "tossPartitioningOperator(" + partitioningParameters.getParameterName() + ")";
    }
}
