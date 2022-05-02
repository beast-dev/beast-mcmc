/*
 * JointOperator.java
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

import dr.math.MathUtils;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */
public class JointOperator extends AbstractAdaptableOperator {

    private final ArrayList<SimpleMCMCOperator> operatorList;
    private final ArrayList<Integer> operatorToOptimizeList;

    private int currentOptimizedOperator;

    public JointOperator(double weight, double targetAcceptanceProbability) {
        super(AdaptationMode.DEFAULT, targetAcceptanceProbability);

        operatorList = new ArrayList<SimpleMCMCOperator>();
        operatorToOptimizeList = new ArrayList<Integer>();

        setWeight(weight);
    }

    public void addOperator(SimpleMCMCOperator operation) {

        operatorList.add(operation);
        if (operation instanceof AdaptableMCMCOperator) {

            if (((AdaptableMCMCOperator) operation).getMode() == AdaptationMode.ADAPTATION_ON)

                operatorToOptimizeList.add(operatorList.size() - 1);

        }
    }

    public final double doOperation() {

        double logP = 0;

        for (SimpleMCMCOperator operation : operatorList) {

            logP += operation.doOperation();
            // todo After a failure, should not have to complete remaining operations, need to fake their operate();
        }

        return logP;
    }

//    private double old;

    @Override
    protected double getAdaptableParameterValue() {
        if (operatorToOptimizeList.size() > 0) {
            currentOptimizedOperator = operatorToOptimizeList.get(MathUtils.nextInt(operatorToOptimizeList.size()));
            return ((AdaptableMCMCOperator) operatorList.get(currentOptimizedOperator)).getAdaptableParameter();
        }
        throw new IllegalArgumentException();
    }

    public void setAdaptableParameterValue(double value) {
        if (operatorToOptimizeList.size() > 0) {
            ((AdaptableMCMCOperator) operatorList.get(currentOptimizedOperator)).setAdaptableParameter(value);
            return;
        }
        throw new IllegalArgumentException();
    }


    public int getNumberOfSubOperators() {
        return operatorList.size();
    }

    public double getRawParameter(int i) {
        if (i < 0 || i >= operatorList.size())
            throw new IllegalArgumentException();
        return ((AdaptableMCMCOperator) operatorList.get(i)).getRawParameter();
    }


    public double getRawParameter() {
        throw new RuntimeException("More than one raw parameter for a joint operator");
    }

    public AdaptationMode getMode() {
        if (operatorToOptimizeList.size() > 0)
            return AdaptationMode.ADAPTATION_ON;
        return AdaptationMode.ADAPTATION_OFF;
    }

    public MCMCOperator getSubOperator(int i) {
        return operatorList.get(i);
    }

    public AdaptationMode getSubOperatorMode(int i) {
        if (i < 0 || i >= operatorList.size())
            throw new IllegalArgumentException();
        if (operatorList.get(i) instanceof AdaptableMCMCOperator)
            return ((AdaptableMCMCOperator) operatorList.get(i)).getMode();
        return AdaptationMode.ADAPTATION_OFF;
    }

    public String getSubOperatorName(int i) {
        if (i < 0 || i >= operatorList.size())
            throw new IllegalArgumentException();
        return "Joint." + operatorList.get(i).getOperatorName();
    }

    public String getOperatorName() {
//        StringBuffer sb = new StringBuffer("Joint(\n");
//        for(SimpleMCMCOperator operation : operatorList)
//            sb.append("\t"+operation.getOperatorName()+"\n");
//        sb.append(") opt = "+optimizedOperator.getOperatorName());
//        return sb.toString();
        return "JointOperator";
    }

    public Element createOperatorElement(Document d) {
        throw new RuntimeException("not implemented");
    }

    public String getAdaptableParameterName() {
        return "";
    }

    public double getMinimumAcceptanceLevel() {
        double min = getTargetAcceptanceProbability() - 0.2;
        if (min < 0)
            min = 0.01;
        return min;
    }

    public double getMaximumAcceptanceLevel() {
        double max = getTargetAcceptanceProbability() + 0.2;
        if (max > 1)
            max = 0.9;
        return max;
    }

    public double getMinimumGoodAcceptanceLevel() {
        double min = getTargetAcceptanceProbability() - 0.1;
        if (min < 0)
            min = 0.01;
        return min;
    }

    public double getMaximumGoodAcceptanceLevel() {
        double max = getTargetAcceptanceProbability() + 0.2;
        if (max > 1)
            max = 0.9;
        return max;
    }

    public final String getPerformanceSuggestion() {

//		double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
//		double targetProb = getTargetAcceptanceProbability();
//		dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
//		double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
//		if (prob < getMinimumGoodAcceptanceLevel()) {
//			return "Try setting scaleFactor to about " + formatter.format(sf);
//		} else if (prob > getMaximumGoodAcceptanceLevel()) {
//			return "Try setting scaleFactor to about " + formatter.format(sf);
//		} else return "";
        return "";
    }


}
