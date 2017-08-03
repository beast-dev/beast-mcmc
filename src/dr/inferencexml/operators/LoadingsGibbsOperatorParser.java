/*
 * LoadingsGibbsOperatorParser.java
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

package dr.inferencexml.operators;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MomentDistributionModel;
import dr.inference.model.LatentFactorModel;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.operators.LoadingsGibbsOperator;
import dr.inference.operators.LoadingsGibbsTruncatedOperator;
import dr.xml.*;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/23/14
 * Time: 1:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class LoadingsGibbsOperatorParser extends AbstractXMLObjectParser {
    public static final String LOADINGS_GIBBS_OPERATOR = "loadingsGibbsOperator";
    public static final String WEIGHT = "weight";
    private final String RANDOM_SCAN = "randomScan";
    private final String WORKING_PRIOR = "workingPrior";
    private final String CUTOFF_PRIOR = "cutoffPrior";
    private final String MULTI_THREADED = "multiThreaded";
    private final String NUM_THREADS = "numThreads";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        String weightTemp = (String) xo.getAttribute(WEIGHT);
        double weight = Double.parseDouble(weightTemp);
        LatentFactorModel LFM = (LatentFactorModel) xo.getChild(LatentFactorModel.class);
        DistributionLikelihood prior = (DistributionLikelihood) xo.getChild(DistributionLikelihood.class);
        MomentDistributionModel prior2 = (MomentDistributionModel) xo.getChild(MomentDistributionModel.class);
        DistributionLikelihood cutoffPrior = null;
        if(xo.hasChildNamed(CUTOFF_PRIOR)){
            cutoffPrior = (DistributionLikelihood) xo.getChild(CUTOFF_PRIOR).getChild(DistributionLikelihood.class);
        }
        boolean randomScan = xo.getAttribute(RANDOM_SCAN, true);
        int numThreads = xo.getAttribute(NUM_THREADS, 4);
        MatrixParameterInterface loadings=null;
        if(xo.getChild(MatrixParameterInterface.class)!=null){
            loadings=(MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
        }
        DistributionLikelihood WorkingPrior = null;
        if(xo.getChild(WORKING_PRIOR) != null){
            WorkingPrior = (DistributionLikelihood) xo.getChild(WORKING_PRIOR).getChild(DistributionLikelihood.class);
        }
        boolean multiThreaded = xo.getAttribute(MULTI_THREADED, false);

        if(prior!=null)
        return new LoadingsGibbsOperator(LFM, prior, weight, randomScan, WorkingPrior, multiThreaded, numThreads);  //To change body of implemented methods use File | Settings | File Templates.
        else
            return new LoadingsGibbsTruncatedOperator(LFM, prior2, weight, randomScan, loadings, cutoffPrior);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LatentFactorModel.class),

            new XORRule(
            new ElementRule(DistributionLikelihood.class),
                    new AndRule(
            new ElementRule(MomentDistributionModel.class),

                    new ElementRule(CUTOFF_PRIOR, new XMLSyntaxRule[]{new ElementRule(DistributionLikelihood.class)}))
            ),
//            new ElementRule(CompoundParameter.class),
            AttributeRule.newDoubleRule(WEIGHT),
            AttributeRule.newBooleanRule(MULTI_THREADED, true),
            AttributeRule.newIntegerRule(NUM_THREADS, true),
            new ElementRule(WORKING_PRIOR, new XMLSyntaxRule[]{
                    new ElementRule(DistributionLikelihood.class)
            }, true),
    };

    @Override
    public String getParserDescription() {
        return "Gibbs sampler for the loadings matrix of a latent factor model";  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public Class getReturnType() {
        return LoadingsGibbsOperator.class;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getParserName() {
        return LOADINGS_GIBBS_OPERATOR;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
