/*
 * CompoundLikelihoodParser.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inferencexml.model;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 */
public class CompoundLikelihoodParser extends AbstractXMLObjectParser {
    public static final String COMPOUND_LIKELIHOOD = "compoundLikelihood";
    public static final String THREADS = "threads";
    @Deprecated public static final String POSTERIOR = "posterior"; // replaced with JOINT
    public static final String JOINT = "joint";
    public static final String PRIOR = "prior";
    public static final String LIKELIHOOD = "likelihood";
    public static final String PSEUDO_PRIOR = "pseudoPrior";
    public static final String REFERENCE_PRIOR = "referencePrior";
    public static final String WORKING_PRIOR = "workingPrior";

    public String getParserName() {
        return COMPOUND_LIKELIHOOD;
    }

    public String[] getParserNames() {
        return new String[]{getParserName(), POSTERIOR, JOINT, PRIOR, LIKELIHOOD, PSEUDO_PRIOR, REFERENCE_PRIOR, WORKING_PRIOR};
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        // the default is -1 threads (automatic thread pool size) but an XML attribute can override it
        int threads = xo.getAttribute(THREADS, -1);

        // both the XML attribute and a system property can override it
        if (System.getProperty("thread.count") != null) {

            threads = Integer.parseInt(System.getProperty("thread.count"));
            if (threads < -1 || threads > 1000) {
                // put an upper limit here - may be unnecessary?
                threads = -1;
            }
        }
//        }

        List<Likelihood> likelihoods = new ArrayList<Likelihood>();
        for (int i = 0; i < xo.getChildCount(); i++) {
            final Object child = xo.getChild(i);
            if (child instanceof Likelihood) {

                if (likelihoods.contains(child)) {
                    throw new XMLParseException("The likelihood element, '" + ((Likelihood) child).getId() + "', is already present in the likelihood or prior density.");
                }
                likelihoods.add((Likelihood) child);
                
//            } else if (child instanceof BeagleBranchLikelihoods){
//                
//            	//TODO
//            	likelihoods.addAll( ((BeagleBranchLikelihoods)child).getBranchLikelihoods());
                
            } else {

                throw new XMLParseException("An element (" + child + ") which is not a likelihood has been added to a "
                        + COMPOUND_LIKELIHOOD + " element");
            }
        }

        CompoundLikelihood compoundLikelihood;

        if (xo.getName().equalsIgnoreCase(LIKELIHOOD)) {
            compoundLikelihood = new CompoundLikelihood(threads, likelihoods);
            switch (threads) {
                case -1:
                    Logger.getLogger("dr.evomodel").info("\nLikelihood computation is using an auto sizing thread pool.");
                    break;
                case 0:
                    Logger.getLogger("dr.evomodel").info("\nLikelihood computation is using a single thread.");
                    break;
                default:
                    Logger.getLogger("dr.evomodel").info("\nLikelihood computation is using a pool of " + threads + " threads.");
                    break;
            }
        } else {
            compoundLikelihood = new CompoundLikelihood(likelihoods);
        }


//		TODO
//        System.err.println("CompundLikelihood consists of " + compoundLikelihood.getLikelihoodCount() + " likelihood element(s)");
        
        return compoundLikelihood;
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A likelihood function which is simply the product of its component likelihood functions.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(THREADS, true),
            new ElementRule(Likelihood.class, -1, Integer.MAX_VALUE)
    };

    public Class getReturnType() {
        return CompoundLikelihood.class;
    }
}
