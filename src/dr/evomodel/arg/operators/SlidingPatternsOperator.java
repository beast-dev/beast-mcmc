/*
 * SlidingPatternsOperator.java
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

import dr.evolution.alignment.SitePatterns;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.inference.operators.*;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 18, 2007
 * Time: 8:01:56 PM
 * To change this template use File | Settings | File Templates.
 */
public class SlidingPatternsOperator extends AbstractCoercableOperator {
//
//		SimpleMCMCOperator implements CoercableMCMCOperator {

    public static final String WINDOW_SIZE = "windowSize";
    public static final String OPERATOR_NAME = "slidingPatternsOperator";
    public static final String BREAK_POINTS = "breakPoints";

    public SlidingPatternsOperator(List<SitePatterns> list, Parameter breakPoints, int windowSize, int weight, CoercionMode mode) {
        super(mode);
        this.partitions = list;
        this.windowSize = windowSize;
//		this.weight = weight;
//		this.mode = mode;
        this.breakPoints = breakPoints;
//        System.err.println("Starting values: "+currentBreakPointsString() );
//        System.exit(-1) ;
    }


    public String currentBreakPointsString() {
        StringBuilder sb = new StringBuilder();
        sb.append("[");
        boolean first = true;
        for (double value : breakPoints.getParameterValues()) {
            int pt = (int) value;
            if (!first) {
                sb.append(",");
                first = false;
            }
            sb.append(pt);
        }
        sb.append("]");
        return sb.toString();

    }

    //                http://www.google.com/search?client=safari&rls=en&q=The+week+renewal+add+a+friend&ie=UTF-8&oe=UTF-8
    public void addNewSitePatterns(SitePatterns addPatterns) {
        // todo need to implement for a variable number of partitions
    }

    public String getOperatorName() {
        return OPERATOR_NAME;
    }

    public double doOperation() {

        // Select boundary to update, 0 => btw partition 0 and 1,  1 => btw partition 1 and 2, etc.
        int whichBoundary = MathUtils.nextInt(breakPoints.getDimension());

        int cBreakPt = (int) breakPoints.getParameterValue(whichBoundary);
        SitePatterns left = partitions.get(whichBoundary);
        SitePatterns right = partitions.get(whichBoundary + 1);
        int min = left.getFrom();
        int max = right.getTo();

        int pBreakPt = min;
        while (pBreakPt <= min || pBreakPt >= max) {        // cBreakPt + [windowSize-1, windowSize+1] (and not 0)
            if (MathUtils.nextBoolean())
                pBreakPt = cBreakPt + MathUtils.nextInt(windowSize) + 1;
            else
                pBreakPt = cBreakPt - MathUtils.nextInt(windowSize) - 1;
        }


        return 0;
    }

    public double getCoercableParameter() {
        return Math.log(windowSize);
    }

    public void setCoercableParameter(double value) {
        windowSize = (int) Math.exp(value);
    }

    public double getRawParameter() {
        return windowSize;
    }

//	public int getMode() {
//		return mode;
//	}

    public static boolean arePartitionsContiguous(List<SitePatterns> list) {
        int current = -1;
        int index = 0;
        for (SitePatterns patterns : list) {
            int start = patterns.getFrom();
            int end = patterns.getTo();
            /* System.err.println(start+" -> "+end+" : "+patterns.getSiteCount()+" "+patterns.getPatternCount());
                           int[] data = patterns.getSitePattern(0);
                           System.err.print("Data 0:");
                           for (int i : data)
                               System.err.print(" "+i);
                           System.err.println("");*/
//            if (current == -1)
//                current = end;
            if (current != -1 && start != (current + 1))
//                throw new NonContiguousPartitionsException("Partition #"+0+" does not start contiguously");
                return false;
            current = end;
        }

        return true;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

//	public int getWeight() {
//		return weight;
//	}

//	public void setWeight(int w) {
//		weight = w;
//	}

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();

        double ws = OperatorUtils.optimizeWindowSize(windowSize, totalAlignmentLength, prob, targetProb);

        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try decreasing windowSize to about " + ws;
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try increasing windowSize to about " + ws;
        } else return "";
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return OPERATOR_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			int mode = CoercableMCMCOperator.DEFAULT;
//			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
//				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
//					mode = CoercableMCMCOperator.COERCION_ON;
//				} else {
//					mode = CoercableMCMCOperator.COERCION_OFF;
//				}
//			}

            CoercionMode mode = CoercionMode.parseMode(xo);

            int weight = xo.getIntegerAttribute(WEIGHT);
            int windowSize = xo.getIntegerAttribute(WINDOW_SIZE);

            List<SitePatterns> list = new ArrayList<SitePatterns>();

            final int numChild = xo.getChildCount();
            for (int i = 0; i < numChild; i++) {
                Object obj = xo.getChild(i);
                if (obj instanceof SitePatterns) {
//                    System.err.println("Found: "+((SitePatterns)obj).getId()) ;
                    list.add((SitePatterns) obj);
                }
            }
            if (!arePartitionsContiguous(list))
                throw new XMLParseException("Only contiguous partitions are allowed");

            // Set current breakpoints
            int dim = list.size() - 1;
            XMLObject cxo = xo.getChild(BREAK_POINTS);
            Parameter breakPoints = new Parameter.Default(dim);

            ParameterParser.replaceParameter(cxo, breakPoints);

            breakPoints.setDimension(dim);
            for (int index = 0; index < dim; index++)
                breakPoints.setParameterValueQuietly(index,
                        list.get(index + 1).getFrom());

//            Parameter parameter = (Parameter)xo.getChild(Parameter.class);

            return new SlidingPatternsOperator(list, breakPoints, windowSize, weight, mode);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a sliding window operator on alignment sites.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule(WINDOW_SIZE),
                AttributeRule.newIntegerRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                new ElementRule(BREAK_POINTS, Parameter.class),
                new ElementRule(SitePatterns.class, 2, 100) // todo fix hard-coded 100
        };

    };


    private int totalAlignmentLength;
    private final List<SitePatterns> partitions;
    private int windowSize = 10;
    //	private int mode = CoercableMCMCOperator.DEFAULT;
    //	private int weight = 1;
    private final Parameter breakPoints;
}
