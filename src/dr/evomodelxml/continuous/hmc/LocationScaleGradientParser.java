/*
 * NodeHeightGradientParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.HyperParameterBranchRateGradient;
import dr.evomodel.treedatalikelihood.discrete.LocationGradient;
import dr.evomodel.treedatalikelihood.discrete.ScaleGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.Parameter;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

public class LocationScaleGradientParser extends AbstractXMLObjectParser {

    private static final String NAME = "locationScaleGradient";
    private static final String LOCATION = "location";
    private static final String SCALE = "scale";

    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String USE_HESSIAN = "useHessian";

    public String getParserName(){ return NAME; }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        boolean useHessian = xo.getAttribute(USE_HESSIAN, false);

        final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
//        final Parameter locationScaleParameter = (Parameter) xo.getChild(Parameter.class);
        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

        if (delegate instanceof  ContinuousDataLikelihoodDelegate) {
            throw new XMLParseException("Not yet implemented! ");
        } else if (delegate instanceof  BeagleDataLikelihoodDelegate) {
            BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) delegate;


            if (branchRateModel instanceof DefaultBranchRateModel || branchRateModel instanceof ArbitraryBranchRates) {
                if (xo.hasChildNamed(LOCATION)) {
                    Object locationObject = xo.getElementFirstChild(LOCATION);
                    BranchSpecificFixedEffects location;

                    if (locationObject instanceof Parameter) {
                        location = new BranchSpecificFixedEffects.None((Parameter) xo.getElementFirstChild(LOCATION));
                    } else if (locationObject instanceof BranchSpecificFixedEffects) {
                        location = (BranchSpecificFixedEffects) locationObject;
                    } else {
                        throw new XMLParseException("Poorly formed");
                    }

                    return new LocationGradient(traitName, treeDataLikelihood, beagleData, location, useHessian);

                } else if (xo.hasChildNamed(SCALE)) {

                    Parameter scale = (Parameter) xo.getElementFirstChild(SCALE);
                    return new ScaleGradient(traitName, treeDataLikelihood, beagleData, scale, useHessian);

                } else {
                    throw new XMLParseException("Poorly formed");
                }
            } else {
                throw new XMLParseException("Only implemented for an arbitrary rates model");
            }
        } else {
            throw new XMLParseException("Unknown likelihood delegate type");
        }
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),

            new ElementRule(TreeDataLikelihood.class),

            new XORRule(
                new ElementRule(LOCATION, new XMLSyntaxRule[]{
                        new XORRule(
                                new ElementRule(BranchSpecificFixedEffects.class),
                                new ElementRule(Parameter.class)
                        )
                }),
                new ElementRule(SCALE, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                })
            ),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return HyperParameterBranchRateGradient.class;
    }
}
