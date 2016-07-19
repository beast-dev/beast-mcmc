/*
 * LineageModelLikelihood.java
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

package dr.oldevomodel.lineage;

import dr.inference.model.*;
import dr.xml.*;

/**
 * Package: LineageModelLikelihood
 * Description:
 * <p/>
 * <p/>
 * Created by
 *
 * @author Alexander V. Alekseyenko (alexander.alekseyenko@gmail.com)
 *         Date: 10/14/13
 *         Time: 12:13 PM
 */
public class LineageModelLikelihood extends AbstractModelLikelihood {

    protected int numSamples;
    protected int numSNPs;
    protected int numLineages;
    protected double normalization;
    protected MatrixParameter mixtureMatrix;
    protected LineageSitePatterns patterns;
    protected Parameter errorRate;
    protected MatrixParameter refData;
    protected MatrixParameter nonData;

    public LineageModelLikelihood(LineageSitePatterns patterns, MatrixParameter mixtureMatrix, Parameter errorRate,
                                  MatrixParameter refData, MatrixParameter nonData)
    {
        super(LINEAGE_MODEL);

        numSamples = mixtureMatrix.getRowDimension();
        numSNPs = patterns.getSiteCount();
        numLineages = mixtureMatrix.getColumnDimension();

        this.mixtureMatrix = mixtureMatrix;
        this.patterns = patterns;
        this.errorRate = errorRate;
        this.refData = refData;
        this.nonData = nonData;

        normalization = 0.0;
        for (int j = 0; j < numSNPs; j++){
            for (int i =0; i < numSamples; i++){
                normalization += dr.math.Binomial.logChoose(Math.round(refData.getParameterValue(i,j) +
                        nonData.getParameterValue(i,j)), Math.round(nonData.getParameterValue(i,j)));
            }
        }

        addVariable(errorRate);
        addModel(patterns);
        addVariable(mixtureMatrix);
    }
    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if(model == mixtureMatrix || model == patterns || model == errorRate)
            makeDirty();
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
    }

    @Override
    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    protected double calculateLogLikelihood() {
        double logLike=normalization, p;

        int i, j, k;

        for (j = 0; j < numSNPs; j++)
        {
            for (i =0; i < numSamples; i++)
            {
                p = 0;
                for (k =0; k < numLineages; k++)
                {
                    p += mixtureMatrix.getParameterValue(k,i)*(1-patterns.getState(k, j));
                }
                p = p - 2*errorRate.getParameterValue(0)*p +errorRate.getParameterValue(0);
                logLike += refData.getParameterValue(i,j)*Math.log(p) + nonData.getParameterValue(i,j)*Math.log(1-p);
            }
        }
        return logLike;
    }


    public void makeDirty() {
        likelihoodKnown = false;
    }
    protected boolean likelihoodKnown = false;
    protected double logLikelihood = 0;
    private double storedLogLikelihood;
    private boolean storedLikelihoodKnown = false;

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

//    public Element createElement(Document d) {
//        throw new RuntimeException("Not implemented yet!");
//    }

    public static final String LINEAGE_MODEL = "LINEAGE_MODEL";
    public static final String LINEAGE_MODEL_PARSER = "lineageModel";
    public static final String MIXTURE = "mixture";
    public static final String REFERENCE = "ref";
    public static final String NON_REFERENCE = "non";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LINEAGE_MODEL_PARSER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            MatrixParameter nonData=null, refData=null, mixtureMatrix=null;
            Parameter errorRate=null;
            LineageSitePatterns patterns=null;

            for(int i=0; i<xo.getChildCount(); ++i){
                if(xo.getChild(i) instanceof Parameter)
                    errorRate = (Parameter) xo.getChild(i);
                else if(xo.getChild(i) instanceof LineageSitePatterns)
                    patterns = (LineageSitePatterns) xo.getChild(i);
            }
            mixtureMatrix = (MatrixParameter)xo.getElementFirstChild(MIXTURE);
            refData = (MatrixParameter)xo.getElementFirstChild(REFERENCE);
            nonData = (MatrixParameter)xo.getElementFirstChild(NON_REFERENCE);

            if(errorRate==null){
                throw new XMLParseException("An element of class Parameter corresponding to error rate needs to be provided.");
            }
            if(patterns == null){
                throw new XMLParseException("Lineage model-compatible site patterns need to be provided.");
            }

            if(nonData.getColumnDimension() != refData.getColumnDimension() || nonData.getRowDimension() != refData.getRowDimension() || mixtureMatrix.getRowDimension() != refData.getColumnDimension()) {
                System.err.println("REF " + refData.getRowDimension() + " x " + refData.getColumnDimension() + "\n");
                System.err.println("NON " + nonData.getRowDimension() + " x " + nonData.getColumnDimension() + "\n");
                System.err.println("MIXTURE " + mixtureMatrix.getRowDimension() + " x " + mixtureMatrix.getColumnDimension() + "\n");
                throw new XMLParseException("Some dimensions do not match, check your input data.");
            }

            return new LineageModelLikelihood(patterns, mixtureMatrix, errorRate, refData, nonData);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A matrix parameter constructed from its component parameters.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MIXTURE, new XMLSyntaxRule[]{
                                    new ElementRule(MatrixParameter.class, false)
                            }, false),
                new ElementRule(REFERENCE, new XMLSyntaxRule[]{
                                                    new ElementRule(MatrixParameter.class, false)
                                            }, false),
                new ElementRule(NON_REFERENCE, new XMLSyntaxRule[]{
                                                    new ElementRule(MatrixParameter.class, false)
                                            }, false),
                new ElementRule(LineageSitePatterns.class),
                new ElementRule(Parameter.class)
        };

        public Class getReturnType() {
            return LineageModelLikelihood.class;
        }
    };
}
