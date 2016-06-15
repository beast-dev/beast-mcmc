/*
 * OrderedLatentLiabilityLikelihood.java
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
import dr.evolution.tree.NodeRef;
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
 * A class to model multivariate ordered data as realizations from a latent (liability) multivariate Brownian diffusion
 *
 * @author Gabrila Cybis
 * @author Marc A. Suchard
 * @author Joe Felsenstein
 * @version $Id$
 */

public class OrderedLatentLiabilityLikelihood extends AbstractModelLikelihood implements LatentTruncation, Citable, SoftThresholdLikelihood {

    public final static String ORDERED_LATENT_LIABILITY_LIKELIHOOD = "orderedLatentLiabilityLikelihood";

    public OrderedLatentLiabilityLikelihood(TreeModel treeModel, PatternList patternList, CompoundParameter tipTraitParameter, CompoundParameter thresholdParameter, Parameter numClasses, boolean isUnordered) {
        super(ORDERED_LATENT_LIABILITY_LIKELIHOOD);
        this.treeModel = treeModel;
        this.patternList = patternList;
        this.tipTraitParameter = tipTraitParameter;
        this.thresholdParameter = thresholdParameter;
        this.numClasses = numClasses;
        this.isUnordered = isUnordered;

        addVariable(tipTraitParameter);
        addVariable(thresholdParameter);


        for (int i=0; i<tipTraitParameter.getParameterCount();i++){
            tipTraitParameter.getParameter(i).addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY,Double.NEGATIVE_INFINITY, tipTraitParameter.getParameter(0).getDimension()));
        }



        setTipDataValuesForAllNodes();

        StringBuilder sb = new StringBuilder();
        sb.append("Constructing a latent liability likelihood model:\n");
        sb.append("\tBinary patterns: ").append(patternList.getId()).append("\n");
        sb.append("\tPlease cite:\n").append(Citable.Utils.getCitationString(this));
        Logger.getLogger("dr.evomodel.continous").info(sb.toString());
    }

    private void setTipDataValuesForAllNodes() {
        if (tipData == null) {
            tipData = new int[treeModel.getExternalNodeCount()][patternList.getPatternCount()];
        }

        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef node = treeModel.getExternalNode(i);
            String id = treeModel.getTaxonId(i);
            int index = patternList.getTaxonIndex(id);
            setTipDataValuesForNode(node, index);

            System.err.println("\t For node: " + i + " with ID " + id + " you get taxon " + index + " with ID " + patternList.getTaxonId(index));


        }
    }

    private void setTipDataValuesForNode(NodeRef node, int index) {
        // Set tip data values
        int Nindex = node.getNumber();
        //   if (index != indexFromPatternList) {
        //       throw new RuntimeException("Need to figure out the indexing");
        //  }

        for (int datum = 0; datum < patternList.getPatternCount(); ++datum) {
            tipData[Nindex][datum] = (int) patternList.getPattern(datum)[index];
            if (DEBUG) {
                Parameter oneTipTraitParameter = tipTraitParameter.getParameter(Nindex);
                System.err.println("Data = " + tipData[Nindex][datum] + " : " + oneTipTraitParameter.getParameterValue(datum));


            }
        }
    }

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

    @Override
    public void setPathParameter(double beta) {
        pathParameter=beta;
    }

    @Override
    public double getLikelihoodCorrection() {
        boolean valid=true;
        for (int tip = 0; tip < tipData.length && valid; ++tip) {
            valid = validTraitForTip(tip);

        }
        if(valid)
        {
            return 0;
        }
        else{
            return 1/(1-pathParameter);
        }
    }

    public String toString() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    protected double computeLogLikelihood() {
        boolean valid = true;

        for (int tip = 0; tip < tipData.length && valid; ++tip) {
            valid = validTraitForTip(tip);


        }



        if (valid) {
            return 0.0;
        } else {
            if(pathParameter==1)
                return Double.NEGATIVE_INFINITY;
            else{
                return Math.log(1-pathParameter);
            }
        }
    }

    public boolean validTraitForTip(int tip) {
        boolean valid = true;
        Parameter oneTipTraitParameter = tipTraitParameter.getParameter(tip);
        int[] data = tipData[tip];


       if (!isUnordered){




        int threshNum = 0;

        for (int index = 0; index < data.length && valid; ++index) {

            int datum = data[index];
            double trait = oneTipTraitParameter.getParameterValue(index);
            int dim = (int) numClasses.getParameterValue(index);




            if (dim == 1.0) {
                valid = true;
            } else if (dim == 2.0){
                if (trait == 0) {
                    valid = true;
                }
                else if(datum>1){
                    valid=true;
                } else {
                    boolean positive = trait > 0.0;
                    if (positive) {
                        valid = (datum == 1.0);
                    } else {
                        valid = (datum == 0.0);
                    }
                }
            } else {
                if (datum == 0) {
                    valid = trait <= 0.0;
                } else if (datum == 1) {
                    valid = (trait >= 0 &&
                            trait <= thresholdParameter.getParameter(threshNum).getParameterValue(0));
                } else if (datum == (dim - 1)) {
                    valid = trait >= thresholdParameter.getParameter(threshNum).getParameterValue(dim - 3);
                } else if (datum > (dim-1)) {
                    valid=true;
                }else {
                    valid = (trait >= thresholdParameter.getParameter(threshNum).getParameterValue(datum - 2) &&
                            trait <= thresholdParameter.getParameter(threshNum).getParameterValue(datum - 1));
                }

                threshNum++;
            }
        }
       }else{
           int LLpointer = 0;

           for (int index = 0; index < data.length && valid; ++index) {

               int datum = data[index];
               int dim = (int) numClasses.getParameterValue(index);


               if (dim == 1.0) {
                   valid = true;
                   LLpointer++;
               } else if (dim == 2.0) {

                   double trait = oneTipTraitParameter.getParameterValue(LLpointer);


                   if (trait == 0) {
                       valid = true;
                   } else {


                       boolean positive = trait > 0.0;
                       if (positive) {
                           valid = (datum == 1.0);
                       } else {
                           valid = (datum == 0.0);
                       }
                   }
                   LLpointer++;
               } else {
                   double[] trait = new double[dim];
                   trait[0] = 0.0;
                   for (int l=1; l<dim; l++){
                       trait[l]=oneTipTraitParameter.getParameterValue(LLpointer+l-1);
                   }

                   valid=isMax(trait, datum);


                   LLpointer += dim-1;


               }
           }
       }

        return valid;
    }



  /*
    public boolean unorderedValidTraitForTip(int tip) {
        boolean valid = true;
        Parameter oneTipTraitParameter = tipTraitParameter.getParameter(tip);
        int[] data = tipData[tip];
        int LLpointer = 0;

        for (int index = 0; index < data.length && valid; ++index) {

            int datum = data[index];
            int dim = (int) numClasses.getParameterValue(index);


            if (dim == 1.0) {
                valid = true;
                LLpointer++;
            } else if (dim == 2.0) {

                double trait = oneTipTraitParameter.getParameterValue(LLpointer);


                if (trait == 0) {
                    valid = true;
                } else {


                    boolean positive = trait > 0.0;
                    if (positive) {
                        valid = (datum == 1.0);
                    } else {
                        valid = (datum == 0.0);
                    }
                }
                LLpointer++;
            } else {
                double[] trait = new double[dim];
                trait[0] = 0.0;
                for (int l=1; l<dim; l++){
                    trait[l]=oneTipTraitParameter.getParameterValue(LLpointer+l-1);
                }

                valid=isMax(trait, datum);


                LLpointer += dim-1;


            }
        }
        return valid;
    }
  */

    private boolean isMax(double[] trait, int datum) {

        boolean isMax = true;


        for (int j = 0; j < trait.length && isMax; j++) {
            isMax = (trait[datum] >= trait[j]);
        }


        return isMax;
    }






    public double getNormalizationConstant(Distribution working) {
        return normalizationDelegate.getNormalizationConstant(working); // delegate to abstract Delegate
    }

    private final LatentTruncation.Delegate normalizationDelegate = new Delegate() {

        protected double computeNormalizationConstant(Distribution working) {
            double constant = 0.0;
            // TODO
            return constant;
        }
    };
    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public final static String TIP_TRAIT = "tipTrait";
        public final static String THRESHOLD_PARAMETER = "threshold";
        public final static String NUM_CLASSES = "numClasses";
        public static final String IS_UNORDERED = "isUnordered";
        public final static String N_DATA="NData";
        public final static String N_TRAITS="NTraits";


        public String getParserName() {
            return ORDERED_LATENT_LIABILITY_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            int numData;
            int dimTrait;
            if(xo.hasAttribute(N_DATA)&&xo.hasAttribute(N_TRAITS)){
                String numDataTemp=(String) xo.getAttribute(N_DATA);
                String dimTraitTemp=(String) xo.getAttribute(N_TRAITS);
                numData= Integer.parseInt(numDataTemp);
                dimTrait= Integer.parseInt(dimTraitTemp);
            }
            else {
                AbstractMultivariateTraitLikelihood traitLikelihood = (AbstractMultivariateTraitLikelihood)
                        xo.getChild(AbstractMultivariateTraitLikelihood.class);
                numData = traitLikelihood.getNumData();
                dimTrait = traitLikelihood.getDimTrait();
            }
            PatternList patternList = (PatternList) xo.getChild(PatternList.class);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
            CompoundParameter thresholdParameter = (CompoundParameter) xo.getElementFirstChild(THRESHOLD_PARAMETER);
            Parameter numClasses = (Parameter) xo.getElementFirstChild(NUM_CLASSES);
            boolean isUnorderd = xo.getAttribute(IS_UNORDERED,false);


            int numTaxa = treeModel.getTaxonCount();




                if (tipTraitParameter.getDimension() != numTaxa * numData * dimTrait) {
                    throw new XMLParseException("Tip trait parameter is wrong dimension in latent liability model");
                }

            if (!isUnorderd){

                if (patternList.getPatternCount() != numData * dimTrait) {
                    throw new XMLParseException("Data is wrong dimension in latent liability model");
                }
            }



            return new OrderedLatentLiabilityLikelihood(treeModel, patternList, tipTraitParameter, thresholdParameter, numClasses,isUnorderd);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Provides the likelihood of a latent liability model on multivariate ordered trait data";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new OrRule(
                new ElementRule(AbstractMultivariateTraitLikelihood.class, "The model for the latent random variables"),
                new AndRule(AttributeRule.newIntegerRule(N_DATA), AttributeRule.newIntegerRule(N_TRAITS))
                        ),
                new ElementRule(TIP_TRAIT, CompoundParameter.class, "The parameter of tip locations from the tree"),
                new ElementRule(THRESHOLD_PARAMETER, CompoundParameter.class, "The parameter with nonzero thershold values"),
                new ElementRule(NUM_CLASSES, Parameter.class, "Number of multinomial classes in each dimention"),
                new ElementRule(PatternList.class, "The binary tip data"),
                new ElementRule(TreeModel.class, "The tree model"),
                AttributeRule.newBooleanRule(IS_UNORDERED,true),
        };

        public Class getReturnType() {
            return OrderedLatentLiabilityLikelihood.class;
        }
    };

    @Override
    public String getCategory() {
        return "Trait Model";
    }

    @Override
    public String getDescription() {
        return "Latent Liability model";
    }

    @Override
    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(CommonCitations.SUCHARD_2012_LATENT);
        return citations;
    }

    public TreeModel treeModel;
    private PatternList patternList;
    public CompoundParameter tipTraitParameter;
    private CompoundParameter thresholdParameter;
    public Parameter numClasses;
    private Parameter containsMissing;

    private boolean isUnordered = false;
    private int[][] tipData;


    private boolean likelihoodKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    private static final boolean DEBUG = false;

    private double pathParameter=1;
}