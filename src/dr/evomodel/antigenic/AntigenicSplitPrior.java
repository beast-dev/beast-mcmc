/*
 * AntigenicSplitPrior.java
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

package dr.evomodel.antigenic;

import dr.inference.model.*;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @author Marc Suchard
 * @version $Id$
 */
@Deprecated // drift prior incorporated into antigenicLikelihood
public class AntigenicSplitPrior extends AbstractModelLikelihood implements Citable {

    public final static String ANTIGENIC_SPLIT_PRIOR = "antigenicSplitPrior";

    public AntigenicSplitPrior(
            MatrixParameter locationsParameter,
            Parameter datesParameter,
            Parameter regressionSlopeParameter,
            Parameter regressionPrecisionParameter,
            Parameter splitTimeParameter,
            Parameter splitAngleParameter,
            Parameter splitAssignmentParameter,
            List<String> topBranchList,
            List<String> bottomBranchList
    ) {

        super(ANTIGENIC_SPLIT_PRIOR);

        this.locationsParameter = locationsParameter;
        addVariable(this.locationsParameter);

        this.datesParameter = datesParameter;
        addVariable(this.datesParameter);

        dimension = locationsParameter.getParameter(0).getDimension();
        count = locationsParameter.getParameterCount();

        earliestDate = datesParameter.getParameterValue(0);
        double latestDate = datesParameter.getParameterValue(0);
        for (int i=0; i<count; i++) {
            double date = datesParameter.getParameterValue(i);
            if (earliestDate > date) {
                earliestDate = date;
            }
            if (latestDate < date) {
                latestDate = date;
            }
        }
        double timeSpan = latestDate - earliestDate;

        this.regressionSlopeParameter = regressionSlopeParameter;
        addVariable(regressionSlopeParameter);
        regressionSlopeParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.regressionPrecisionParameter = regressionPrecisionParameter;
        addVariable(regressionPrecisionParameter);
        regressionPrecisionParameter.addBounds(new Parameter.DefaultBounds(Double.MAX_VALUE, 0.0, 1));

        this.splitTimeParameter = splitTimeParameter;
        addVariable(splitTimeParameter);
        splitTimeParameter.addBounds(new Parameter.DefaultBounds(50.0, 20.0, 1));

        this.splitAngleParameter = splitAngleParameter;
        addVariable(splitAngleParameter);
        splitAngleParameter.addBounds(new Parameter.DefaultBounds(0.5*Math.PI, 0.01, 1));

        this.splitAssignmentParameter = splitAssignmentParameter;
        addVariable(splitAssignmentParameter);
        splitAssignmentParameter.addBounds(new Parameter.DefaultBounds(1.0, 0.0, 1));

        String[] labelArray = new String[count];
        splitAssignmentParameter.setDimension(count);

        topBranchIndices = new int[topBranchList.size()];
        bottomBranchIndices = new int[bottomBranchList.size()];

        for (int i = 0; i < count; i++) {
            String iName = datesParameter.getDimensionName(i);
            labelArray[i] = iName;
            splitAssignmentParameter.setParameterValueQuietly(i, 0.0);

            // top branch is 1
            for (int j=0; j < topBranchList.size(); j++) {
                String jName = topBranchList.get(j);
                if (jName.equals(iName)) {
                    topBranchIndices[j] = i;
                    splitAssignmentParameter.setParameterValueQuietly(i, 1.0);
                }
            }

           // bottom branch is 0
           for (int j=0; j < bottomBranchList.size(); j++) {
                String jName = bottomBranchList.get(j);
                if (jName.equals(iName)) {
                    bottomBranchIndices[j] = i;
                    splitAssignmentParameter.setParameterValueQuietly(i, 0.0);
                }
            }

        }

        splitAssignmentParameter.setDimensionNames(labelArray);

        likelihoodKnown = false;

    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (variable == locationsParameter || variable == datesParameter
            || variable == regressionSlopeParameter || variable == regressionPrecisionParameter
            || variable == splitTimeParameter || variable == splitAngleParameter
            || variable == splitAssignmentParameter) {
            likelihoodKnown = false;
        }
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = false;
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
        }
        return logLikelihood;
    }

    private double computeLogLikelihood() {

        double logLikelihood = 0.0;

        double precision = regressionPrecisionParameter.getParameterValue(0);
        logLikelihood += (0.5 * Math.log(precision) * count) - (0.5 * precision * sumOfSquaredResiduals());
        logLikelihood += assignmentLikelihood();

        likelihoodKnown = true;
        return logLikelihood;

    }

    // set likelihood to negative infinity is an assignment is broken
    protected double assignmentLikelihood() {

        double logLikelihood = 0.0;

        // these must be 1
        for (int index : topBranchIndices) {
            int assignment = (int) splitAssignmentParameter.getParameterValue(index);
            if (assignment == 0) {
            //    logLikelihood = Double.NEGATIVE_INFINITY;
                logLikelihood -= 1000;
            }
        }

        // these must be 0
        for (int index : bottomBranchIndices) {
            int assignment = (int) splitAssignmentParameter.getParameterValue(index);
            if (assignment == 1) {
            //    logLikelihood = Double.NEGATIVE_INFINITY;
                logLikelihood -= 1000;
            }
        }

        return logLikelihood;

    }

    // go through each location and compute sum of squared residuals from regression line
    protected double sumOfSquaredResiduals() {

        double ssr = 0.0;

        for (int i=0; i < count; i++) {

            Parameter loc = locationsParameter.getParameter(i);

            double x = loc.getParameterValue(0);
            double y = expectedAG1(i);
            ssr += (x - y) * (x - y);

            if (dimension > 1) {
                x = loc.getParameterValue(1);
                y = expectedAG2(i);
                ssr += (x - y) * (x - y);
            }

            for (int j=2; j < dimension; j++) {
                x = loc.getParameterValue(j);
                ssr += x*x;
            }

        }

        return ssr;
    }

    // given a location index, calculate the expected AG1 value
    protected double expectedAG1(int index) {

        double date = datesParameter.getParameterValue(index);

        double ag1 = 0;
        double time = date - earliestDate;
        double splitTime = splitTimeParameter.getParameterValue(0);
        double splitAngle = splitAngleParameter.getParameterValue(0);
        double beta = regressionSlopeParameter.getParameterValue(0);

        if (time <= splitTime) {
            ag1 = beta * time;
        }
        else {
            ag1 = (beta * splitTime) + (beta * (time-splitTime) * Math.cos(splitAngle));
        }

        return ag1;

    }

    // given a location index, calculate the expected AG2 value of top branch
    protected double expectedAG2(int index) {

        double date = datesParameter.getParameterValue(index);
        int assignment = (int) splitAssignmentParameter.getParameterValue(index);

        double ag2 = 0;
        double time = date - earliestDate;
        double splitTime = splitTimeParameter.getParameterValue(0);
        double splitAngle = splitAngleParameter.getParameterValue(0);
        double beta = regressionSlopeParameter.getParameterValue(0);

        if (time <= splitTime) {
            ag2 = 0;
        }
        else {
            ag2 = beta * (time-splitTime) * Math.sin(splitAngle);
        }

        if (assignment == 1) {
            ag2 = -1*ag2;
        }

        return ag2;

    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
    }

    private final int dimension;
    private final int count;
    private final Parameter datesParameter;
    private final MatrixParameter locationsParameter;
    private final Parameter regressionSlopeParameter;
    private final Parameter regressionPrecisionParameter;
    private final Parameter splitTimeParameter;
    private final Parameter splitAngleParameter;
    private final Parameter splitAssignmentParameter;
    private final int[] topBranchIndices;
    private final int[] bottomBranchIndices;

    private double earliestDate;
    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
    private boolean likelihoodKnown = false;

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String LOCATIONS = "locations";
        public final static String DATES = "dates";
        public final static String REGRESSIONSLOPE = "regressionSlope";
        public final static String REGRESSIONPRECISION = "regressionPrecision";
        public final static String SPLITTIME = "splitTime";
        public final static String SPLITANGLE = "splitAngle";
        public final static String SPLITASSIGNMENT = "splitAssignment";
        public static final String TOPBRANCH = "topBranch";
        public static final String BOTTOMBRANCH = "bottomBranch";

        public String getParserName() {
            return ANTIGENIC_SPLIT_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild(LOCATIONS);
            Parameter datesParameter = (Parameter) xo.getElementFirstChild(DATES);
            Parameter regressionSlopeParameter = (Parameter) xo.getElementFirstChild(REGRESSIONSLOPE);
            Parameter regressionPrecisionParameter = (Parameter) xo.getElementFirstChild(REGRESSIONPRECISION);
            Parameter splitTimeParameter = (Parameter) xo.getElementFirstChild(SPLITTIME);
            Parameter splitAngleParameter = (Parameter) xo.getElementFirstChild(SPLITANGLE);
            Parameter splitAssignmentParameter = (Parameter) xo.getElementFirstChild(SPLITASSIGNMENT);  List<String> virusLocationStatisticList = null;

            List<String> topBranchList = null;
            String[] topBranch = xo.getStringArrayAttribute(TOPBRANCH);
            if (topBranch != null) {
                topBranchList = Arrays.asList(topBranch);
            }

            List<String> bottomBranchList = null;
            String[] bottomBranch = xo.getStringArrayAttribute(BOTTOMBRANCH);
            if (bottomBranch != null) {
                bottomBranchList = Arrays.asList(bottomBranch);
            }

            AntigenicSplitPrior AGDP = new AntigenicSplitPrior(
                locationsParameter,
                datesParameter,
                regressionSlopeParameter,
                regressionPrecisionParameter,
                splitTimeParameter,
                splitAngleParameter,
                splitAssignmentParameter,
                topBranchList,
                bottomBranchList);

//            Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Utils.getCitationString(AGL));

            return AGDP;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a vector of coordinates in some multidimensional 'antigenic' space based on an expected relationship with time.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(LOCATIONS, MatrixParameter.class),
                new ElementRule(DATES, Parameter.class),
                new ElementRule(REGRESSIONSLOPE, Parameter.class),
                new ElementRule(REGRESSIONPRECISION, Parameter.class),
                new ElementRule(SPLITTIME, Parameter.class),
                new ElementRule(SPLITANGLE, Parameter.class),
                new ElementRule(SPLITASSIGNMENT, Parameter.class),
                AttributeRule.newStringArrayRule(TOPBRANCH, true, "A list of virus names to assign to the top branch."),
                AttributeRule.newStringArrayRule(BOTTOMBRANCH, true, "A list of virus names to assign to the bottom branch.")
        };

        public Class getReturnType() {
            return ContinuousAntigenicTraitLikelihood.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian Antigenic Cartography framework";
    }

    public List<Citation> getCitations() {
        return Arrays.asList(CommonCitations.BEDFORD_2015_INTEGRATING);
    }
}