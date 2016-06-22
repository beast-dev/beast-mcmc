/*
 * IntervalLatentLiabilityLikelihood.java
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

package dr.evomodel.continuous;

import dr.evolution.alignment.PatternList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.math.distributions.Distribution;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.xml.*;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.logging.Logger;


/**
 * A class to model multivariate unit-interval data as realizations from a latent (liability) multivariate Brownian diffusion
 *
 * @author Marc A. Suchard
 * @author Philippe Lemey
 * @version $Id$
 */

public class IntervalLatentLiabilityLikelihood extends AbstractModelLikelihood implements LatentTruncation, Citable, SoftThresholdLikelihood {

    public final static String LATENT_LIABILITY_LIKELIHOOD = "intervalLatentLiabilityLikelihood";

    public IntervalLatentLiabilityLikelihood(TreeModel treeModel, CompoundParameter tipTraitParameter) {
        super(LATENT_LIABILITY_LIKELIHOOD);
        this.treeModel = treeModel;
        this.patternList = null;
        this.tipTraitParameter = tipTraitParameter;

        addVariable(tipTraitParameter);

        setTipDataValuesForAllNodes();

//        System.err.println("Name: " + tipTraitParameter.getId());
//        System.exit(-1);

        for (int i = 0; i < tipTraitParameter.getParameterCount(); ++i) {
            Parameter p = tipTraitParameter.getParameter(i);
            p.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, p.getDimension()));
        }

//        tipTraitParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY, tipTraitParameter.getDimension()));

        StringBuilder sb = new StringBuilder();
        sb.append("Constructing a unit interval latent liability likelihood model:\n");
//        sb.append("\tBinary patterns: ").append(patternList.getId()).append("\n");
        sb.append("\tPlease cite:\n").append(Utils.getCitationString(this));
        Logger.getLogger("dr.evomodel.continuous").info(sb.toString());
    }

    private void setTipDataValuesForAllNodes() {
        System.err.println(tipTraitParameter.getParameterCount());
        System.err.println(tipTraitParameter.getDimension());
//        System.exit(-1);
        if (tipData == null) {
//            tipData = new int[treeModel.getExternalNodeCount()][patternList.getPatternCount()];
            tipData = new long[tipTraitParameter.getDimension()];
        }


//        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
//            NodeRef node = treeModel.getExternalNode(i);
//            String id = treeModel.getTaxonId(i);
//            int index = patternList.getTaxonIndex(id);
//            setTipDataValuesForNode(node, index);
//        }

//        for (int tip = 0; tip < treeModel.getExternalNodeCount(); ++tip) {
//            System.err.println("Tip #" + tip);
//        Parameter oneTipTraitParameter = tipTraitParameter.getParameter(tip);
//        int[] data = tipData[tip];
        for (int index = 0; index < tipData.length; ++index) {
//            int datum = data[index];
//            double trait = oneTipTraitParameter.getParameterValue(index);
//            valid = Math.round(trait) == datum;
            tipData[index] = Math.round(tipTraitParameter.getParameterValue(index));
//            System.err.print(" " + tipData[index]);

//        }
        }
//        System.exit(-1);
    }

//    private void setTipDataValuesForNode(NodeRef node, int indexFromPatternList) {
//        // Set tip data values
//        int index = node.getNumber();
//        if (index != indexFromPatternList) {
//            throw new RuntimeException("Need to figure out the indexing");
//        }
//
//        for (int datum = 0; datum < patternList.getPatternCount(); ++datum) {
//            tipData[index][datum] = patternList.getPattern(datum)[index]  == 1;
//            if (DEBUG) {
//                Parameter oneTipTraitParameter = tipTraitParameter.getParameter(index);
//                System.err.println("Data = " + tipData[index][datum] + " : " + oneTipTraitParameter.getParameterValue(datum));
//            }
//        }
//    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = true;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        likelihoodKnown = false;
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = computeLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public String toString() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    protected double computeLogLikelihood() {
        boolean valid = true;
//        for (int tip = 0; tip < tipData.length && valid; ++tip) {
//            valid = validTraitForTip(tip);
//        }
        for (int index = 0; index < tipData.length && valid; ++index) {
            double raw = tipTraitParameter.getParameterValue(index);
            long round = Math.round(raw);
            valid = round == tipData[index];
//            System.err.println(tipData[index] + " " + round + " " + raw);
            // TODO Handle missing values
        }
//        System.err.println("valid = " + valid);

        // check
        boolean valid2 = true;
        for (int tip = 0; tip < treeModel.getExternalNodeCount() && valid2; ++tip) {
            if (!validTraitForTip(tip)) {
                valid2 = false;
            }
        }

//        System.err.println(valid + " " + valid2);
        if (valid != valid2) {
            throw new RuntimeException("Error in computing validity of tips values");
        }

        if (valid) {
            return 0.0;
        } else {
//            System.exit(-1);
            return Double.NEGATIVE_INFINITY;
        }
    }

    public boolean validTraitForTip(int tip) {
        boolean valid = true;
        Parameter oneTipTraitParameter = tipTraitParameter.getParameter(tip);
        final int offset = oneTipTraitParameter.getDimension() * tip;
        for (int index = 0; index < oneTipTraitParameter.getDimension() && valid; ++index) {
            double raw = oneTipTraitParameter.getParameterValue(index);
            long round = Math.round(raw);
            valid = round == tipData[index + offset];
        }
        return valid;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String TIP_TRAIT = "tipTrait";

        public String getParserName() {
            return LATENT_LIABILITY_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AbstractMultivariateTraitLikelihood traitLikelihood = (AbstractMultivariateTraitLikelihood)
                    xo.getChild(AbstractMultivariateTraitLikelihood.class);
//            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);

            int numTaxa = treeModel.getTaxonCount();
            int numData = traitLikelihood.getNumData();
            int dimTrait = traitLikelihood.getDimTrait();

            if (tipTraitParameter.getDimension() != numTaxa * numData * dimTrait) {
                throw new XMLParseException("Tip trait parameter is wrong dimension in latent liability model");
            }

//            if (!(patternList.getDataType() instanceof TwoStates)) {
//                throw new XMLParseException("Latent liability model currently only works for binary data");
//            }

//            if (patternList.getPatternCount() != numData * dimTrait) {
//                throw new XMLParseException("Binary data is wrong dimension in latent liability model");
//            }

            return new IntervalLatentLiabilityLikelihood(treeModel, tipTraitParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a latent liability model on multivariate-binary trait data";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(AbstractMultivariateTraitLikelihood.class, "The model for the latent random variables"),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree"),
//                new ElementRule(PatternList.class, "The binary tip data"),
                new ElementRule(TreeModel.class, "The tree model"),
        };

        public Class getReturnType() {
            return IntervalLatentLiabilityLikelihood.class;
        }
    };

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Intervaled latent liability model";
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(CommonCitations.CYBIS_2015_ASSESSING);
        return citations;
    }

    public double getNormalizationConstant(Distribution working) {
        return normalizationDelegate.getNormalizationConstant(working); // delegate to abstract Delegate
    }

    private final LatentTruncation.Delegate normalizationDelegate = new Delegate() {

        protected double computeNormalizationConstant(Distribution working) {
            double constant = 0.0;
            for (long datum : tipData) {
                constant += Math.log(working.cdf(datum + 0.5) - working.cdf(datum - 0.5));
            }
            return -constant; // Note minus sign
//            return 16.30411;
        }
    };

    public void setPathParameter(double beta){
        pathParameter=beta;
    }

    @Override
    public double getLikelihoodCorrection() {
        return 0;
    }

    private TreeModel treeModel;
    private PatternList patternList;
    private CompoundParameter tipTraitParameter;

    private long[] tipData;

    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private static final boolean DEBUG = false;

    private double pathParameter=1;

}