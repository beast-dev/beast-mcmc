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

import dr.evomodel.branchmodel.ArbitrarySubstitutionParameterBranchModel;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.BranchSpecificFixedEffects;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.*;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.SumDerivative;
import dr.inference.model.BranchParameter;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

import static dr.evomodelxml.continuous.hmc.BranchRateGradientParser.checkBranchRateModels;
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

        final Object child = xo.getChild(TreeDataLikelihood.class);

        if (child != null) {
            return parseTreeDataLikelihood(xo, (TreeDataLikelihood) child, traitName, useHessian);
        } else {

            CompoundLikelihood compoundLikelihood = (CompoundLikelihood) xo.getChild(CompoundLikelihood.class);
            List<GradientWrtParameterProvider> providers = new ArrayList<>();

            for (Likelihood likelihood : compoundLikelihood.getLikelihoods()) {
                if (!(likelihood instanceof TreeDataLikelihood)) {
                    throw new XMLParseException("Unknown likelihood type");
                }

                GradientWrtParameterProvider provider = parseTreeDataLikelihood(xo, (TreeDataLikelihood) likelihood,
                        traitName, useHessian);

                providers.add(provider);
            }

            checkBranchRateModels(providers);

            return new SumDerivative(providers);
        }
    }

    private GradientWrtParameterProvider parseTreeDataLikelihood(XMLObject xo, TreeDataLikelihood treeDataLikelihood,
                                                                 String traitName,
                                                                 boolean useHessian) throws XMLParseException {


        BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();

        if (delegate instanceof  ContinuousDataLikelihoodDelegate) {
            throw new XMLParseException("Not yet implemented! ");
        } else if (delegate instanceof  BeagleDataLikelihoodDelegate) {
            BeagleDataLikelihoodDelegate beagleData = (BeagleDataLikelihoodDelegate) delegate;

            BranchModel branchModel = beagleData.getBranchModel();

            if (branchRateModel instanceof DefaultBranchRateModel || branchRateModel instanceof ArbitraryBranchRates) {
                if (xo.hasChildNamed(LOCATION)) {

                    BranchSpecificFixedEffects location = parseLocation(xo);

                    return new LocationGradient(traitName, treeDataLikelihood, beagleData, location, useHessian);

                } else if (xo.hasChildNamed(SCALE)) {

                    Parameter scale = (Parameter) xo.getElementFirstChild(SCALE);
                    return new ScaleGradient(traitName, treeDataLikelihood, beagleData, scale, useHessian);

                } else {
                    throw new XMLParseException("Poorly formed");
                }
            } else if (branchModel instanceof ArbitrarySubstitutionParameterBranchModel){

                BranchParameter branchParameter = (BranchParameter) xo.getChild(BranchParameter.class);

                if (xo.hasChildNamed(LOCATION)) {

                    BranchSpecificFixedEffects location = parseLocation(xo);

                    return new BranchSubstitutionParameterLocationGradient(traitName, treeDataLikelihood, beagleData, branchParameter,
                            useHessian, location);

                } else if (xo.hasChildNamed(SCALE)) {

                    Parameter scale = (Parameter) xo.getElementFirstChild(SCALE);
                    return new BranchSubstitutionParameterScaleGradient(traitName, treeDataLikelihood, beagleData, branchParameter, scale, useHessian);

                } else {
                    throw new XMLParseException("Not yet implemented.");
                }

            } else {
                throw new XMLParseException("Only implemented for an arbitrary rates model");
            }
        } else {
            throw new XMLParseException("Unknown likelihood delegate type");
        }
    }

    private BranchSpecificFixedEffects parseLocation(XMLObject xo) throws XMLParseException {
        Object locationObject = xo.getElementFirstChild(LOCATION);
        BranchSpecificFixedEffects location;

        if (locationObject instanceof Parameter) {
            location = new BranchSpecificFixedEffects.None((Parameter) xo.getElementFirstChild(LOCATION));
        } else if (locationObject instanceof BranchSpecificFixedEffects) {
            location = (BranchSpecificFixedEffects) locationObject;
        } else if (locationObject instanceof ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal) {
            location = ((ArbitraryBranchRates.BranchRateTransform.LocationScaleLogNormal) locationObject).getLocationObject();
        } else {
            throw new XMLParseException("Poorly formed");
        }
        return location;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(TRAIT_NAME),

            new XORRule(
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(CompoundLikelihood.class)
            ),

            new XORRule(
                new ElementRule(LOCATION, new XMLSyntaxRule[]{
                        new XORRule(
                                new XMLSyntaxRule[]{
                                        new ElementRule(BranchSpecificFixedEffects.class),
                                        new ElementRule(Parameter.class),
                                        new ElementRule(ArbitraryBranchRates.BranchRateTransform.class)}
                        )
                }),
                new ElementRule(SCALE, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class),
                })
            ),

            new ElementRule(Parameter.class, true)
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
