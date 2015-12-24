/*
 * CompoundGaussianProcessParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.AbstractDistributionLikelihood;
import dr.inference.distribution.CachedDistributionLikelihood;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Variable;
import dr.math.distributions.CompoundGaussianProcess;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.util.Attribute;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class CompoundGaussianProcessParser extends AbstractXMLObjectParser {

    public static final String NAME = "compoundGaussianProcess";

    public String getParserName() {
        return NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        List<GaussianProcessRandomGenerator> gpList = new ArrayList<GaussianProcessRandomGenerator>();
        List<Likelihood> likelihoodList = new ArrayList<Likelihood>();
        List<Integer> copyList = new ArrayList<Integer>();

        for (int i = 0; i < xo.getChildCount(); ++i) {
            Object obj = xo.getChild(i);
            GaussianProcessRandomGenerator gp = null;
            Likelihood likelihood = null;
            int copies = -1;
            if (obj instanceof DistributionLikelihood) {
                DistributionLikelihood dl = (DistributionLikelihood) obj;
                if (!(dl.getDistribution() instanceof GaussianProcessRandomGenerator)) {
                    throw new XMLParseException("Not a Gaussian process");
                }
                likelihood = dl;
                gp = (GaussianProcessRandomGenerator) dl.getDistribution();
                copies = 0;
                for (Attribute<double[]> datum : dl.getDataList()) {
//                    Double draw = (Double) gp.nextRandom();
//                    System.err.println("DL: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + "1");
                    copies += datum.getAttributeValue().length;
                }
            } else if (obj instanceof MultivariateDistributionLikelihood) {
                MultivariateDistributionLikelihood mdl = (MultivariateDistributionLikelihood) obj;
                if (!(mdl.getDistribution() instanceof GaussianProcessRandomGenerator)) {
                    throw new XMLParseException("Not a Gaussian process");
                }
                likelihood = mdl;
                gp = (GaussianProcessRandomGenerator) mdl.getDistribution();
                copies = 0;
                double[] draw = (double[]) gp.nextRandom();
                for (Attribute<double[]> datum : mdl.getDataList()) {
//                    System.err.println("ML: " + datum.getAttributeName() + " " + datum.getAttributeValue().length + " " + draw.length);
                    copies += datum.getAttributeValue().length / draw.length;
                }
            } else if (obj instanceof GaussianProcessRandomGenerator) {
                gp = (GaussianProcessRandomGenerator) obj;
                likelihood = gp.getLikelihood();
                copies = 1;
            } else {
                throw new XMLParseException("Not a Gaussian process");
            }
            gpList.add(gp);
            likelihoodList.add(likelihood);
            copyList.add(copies);
        }

//        System.exit(-1);
        return new CompoundGaussianProcess(gpList, likelihoodList, copyList);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(GaussianProcessRandomGenerator.class, 1, Integer.MAX_VALUE),
    };

    public String getParserDescription() {
        return "Returned a Gaussian process formed from an ordered list of independent Gaussian processes";
    }

    public Class getReturnType() {
        return GaussianProcessRandomGenerator.class;
    }
}
