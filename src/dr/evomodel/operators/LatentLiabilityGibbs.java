


/*
 * LatentLiabilityGibbs.java
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


/*
 * Gabriela Cybis
 *
 *   Gibbs operator for latent variable in latent liability model
 */

package dr.evomodel.operators;

import dr.evolution.tree.MutableTreeModel;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.LatentTruncation;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.math.matrixAlgebra.Vector;
import dr.xml.*;
import dr.math.matrixAlgebra.Matrix;


import java.util.logging.Logger;

//import dr.evomodel.continuous.MultinomialLatentLiabilityLikelihood;
//import dr.inference.operators.GibbsOperator;


public class LatentLiabilityGibbs extends SimpleMCMCOperator {
//	public class LatentLiabilityGibbs  extends SimpleMCMCOperator implements GibbsOperator {

    public static final String LATENT_LIABILITY_GIBBS_OPERATOR = "latentLiabilityGibbsOperator";
    public static final String TREE_MODEL = "treeModel";


    private final LatentTruncation latentLiability;
//	    private final MultinomialLatentLiabilityLikelihood latentLiability;

    private final FullyConjugateMultivariateTraitLikelihood traitModel;
    private final CompoundParameter tipTraitParameter;


    protected double[] rootPriorMean;
    protected double rootPriorSampleSize;

    private final MatrixParameter precisionParam;
    private final MutableTreeModel treeModel;
    private final int dim;

    public double[][] postMeans;
    public double[][] preMeans;
    public double[] preP;
    public double[] postP;
    private Parameter mask;
    private boolean hasMask=false;
    //private int estourou=0;
    //private double ultimoAccept=0;
    private   int numFixed = 0;
    private int numUpdate = 0;
    private int[] doUpdate;
    private int[] dontUpdate;



    public LatentLiabilityGibbs(
            FullyConjugateMultivariateTraitLikelihood traitModel,
            LatentTruncation LatentLiability, CompoundParameter tipTraitParameter,  Parameter mask,
            double weight) {
        super();


        /// Set up all required variables


        this.latentLiability = LatentLiability;
        this.traitModel = traitModel;
        this.tipTraitParameter = tipTraitParameter;


        this.rootPriorMean = traitModel.getPriorMean();
        this.rootPriorSampleSize = traitModel.getPriorSampleSize();

        this.precisionParam = (MatrixParameter) traitModel.getDiffusionModel().getPrecisionParameter();
        this.treeModel = traitModel.getTreeModel();
        dim = precisionParam.getRowDimension(); // assumed to be square

        this.mask=mask;
        if(mask!=null) {
            hasMask=true;
          //  printInformation(mask);
        }


        postMeans = new double[treeModel.getNodeCount()][dim];
        preMeans = new double[treeModel.getNodeCount()][dim];
        preP = new double[treeModel.getNodeCount()];
        postP = new double[treeModel.getNodeCount()];




        dontUpdate=new int[dim];
        doUpdate=new int[dim];

        if (hasMask) {
             for (int i = 0; i < dim; i++) {
                if (mask.getParameterValue(i) == 0.0) {
                    dontUpdate[numFixed] = i;
                    numFixed++;
                } else {
                    doUpdate[numUpdate] = i;
                    numUpdate++;
                }

            }
        }

        setWeight(weight);
    }


    public int getStepCount() {
        return 1;
    }


    private void printInformation(MatrixParameter par) {
        StringBuffer sb = new StringBuffer("\n \n parameter \n");
        for (int j = 0; j < dim; j++) {
            sb.append(par.getParameterValue(0, j));
        }

        Logger.getLogger("dr.evomodel").info(sb.toString());
    }




    private void printInformation(double[] par) {
        StringBuffer sb = new StringBuffer("\n \n double vector \n");
        for (int j = 0; j < treeModel.getNodeCount(); j++) {
            sb.append(par[j]);
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }





    private void printInformation(double[][] par) {
        StringBuffer sb = new StringBuffer("\n \n double matrix \n");
        for (int i = 0; i < 1; i++) {
            for (int j = 0; j < treeModel.getNodeCount(); j++) {

                sb.append(par[j][i]);
            }
        }
        Logger.getLogger("dr.evomodel").info(sb.toString());
    }


    private void printInformation(double par) {
        StringBuffer sb = new StringBuffer("\n \n double \n");

        sb.append(par);

        Logger.getLogger("dr.evomodel").info(sb.toString());
    }


    private void printInformation(double par, String lala) {
        StringBuffer sb = new StringBuffer("\n");

        sb.append(lala);
        sb.append("\t\t");
        sb.append(par);

        Logger.getLogger("dr.evomodel").info(sb.toString());
    }


    public double doOperation() {


//        doPostOrderTraversal(treeModel.getRoot());
//        doPreOrderTraversal(treeModel.getRoot());
//printInformation(postP);
//printInformation(preP);
//printInformation(postMeans);
//printInformation(preMeans);


        final int pos = MathUtils.nextInt(treeModel.getExternalNodeCount());
//	 final int pos = 1;

        NodeRef node = treeModel.getExternalNode(pos);


        double logq = sampleNode2(node);


        tipTraitParameter.fireParameterChangedEvent();

        return logq;
    }


    //Fill out partial mean and precision values in post order
    public void doPostOrderTraversal(NodeRef node) {          // TODO This is already computed IntegratedMultivariateTraitLikelihood

        final int thisNumber = node.getNumber();

        if (treeModel.isExternal(node)) {

            // writes trait values and precision values for tips
            double[] traitValue = getNodeTrait(node);

            for (int j = 0; j < dim; j++) {
                postMeans[thisNumber][j] = traitValue[j];
            }

            postP[thisNumber] = 1 / traitModel.getRescaledBranchLengthForPrecision(node);


            return;

        }


        final NodeRef childNode0 = treeModel.getChild(node, 0);
        final NodeRef childNode1 = treeModel.getChild(node, 1);

        doPostOrderTraversal(childNode0);
        doPostOrderTraversal(childNode1);


        if (!treeModel.isRoot(node)) {

            final int childNumber0 = childNode0.getNumber();
            final int childNumber1 = childNode1.getNumber();


            // precision values
            final double precision0 = postP[childNumber0];
            final double precision1 = postP[childNumber1];
            final double thisPrecision = 1 / traitModel.getRescaledBranchLengthForPrecision(node);
            double tp = precision0 + precision1;
            postP[thisNumber] = tp * thisPrecision / (tp + thisPrecision);


            //mean values
            for (int j = 0; j < dim; j++) {
                postMeans[thisNumber][j] = (precision0 * postMeans[childNumber0][j] + precision1 * postMeans[childNumber1][j]) / (precision0 + precision1);
            }
        }

    }


    public double[] getNodeTrait(NodeRef node) {
        int index = node.getNumber();
        double[] traitValue = tipTraitParameter.getParameter(index).getParameterValues();
        return traitValue;
    }


    public double getNodeTrait(NodeRef node, int entry) {
        int index = node.getNumber();
        double traitValue = tipTraitParameter.getParameter(index).getParameterValue(entry);
        return traitValue;
    }


    public void setNodeTrait(NodeRef node, double[] traitValue) {
        int index = node.getNumber();
        for (int i = 0; i < dim; i++) {

            tipTraitParameter.getParameter(index).setParameterValue(i, traitValue[i]);
        }
        traitModel.getTraitParameter().getParameter(index).fireParameterChangedEvent();
    }


    public void setNodeTrait(NodeRef node, int entry, double traitValue) {
        int index = node.getNumber();

        tipTraitParameter.getParameter(index).setParameterValue(entry, traitValue);

        //	 traitModel.getTraitParameter().getParameter(index).fireParameterChangedEvent();
    }


    //Fill out partial mean and precision values in pre order


    public void doPreOrderTraversal(NodeRef node) {     // TODO This should be computed IntegratedMultivariateTraitLikelihood


        final int thisNumber = node.getNumber();


        if (treeModel.isRoot(node)) {
            preP[thisNumber] = rootPriorSampleSize;
            for (int j = 0; j < dim; j++) {
                preMeans[thisNumber][j] = rootPriorMean[j];
            }


        } else {

            final NodeRef parentNode = treeModel.getParent(node);
            final NodeRef sibNode = getSisterNode(node);

            final int parentNumber = parentNode.getNumber();
            final int sibNumber = sibNode.getNumber();



	/*

			  if (treeModel.isRoot(parentNode)){
				  //partial precisions
				    final double precisionParent = rootPriorSampleSize;
			        final double precisionSib = postP[sibNumber];
			        final double thisPrecision=1/treeModel.getBranchLength(node);
			        double tp= precisionParent + precisionSib;
			        preP[thisNumber]= tp*thisPrecision/(tp+thisPrecision);

			        //partial means

			        for (int j =0; j<dim;j++){
			        	preMeans[thisNumber][j] = (precisionParent*preMeans[parentNumber][j] + precisionSib*rootPriorMean[j])/(precisionParent+precisionSib);
			        }

			  }else{
	*/
            //partial precisions
            final double precisionParent = preP[parentNumber];
            final double precisionSib = postP[sibNumber];
            final double thisPrecision = 1 / traitModel.getRescaledBranchLengthForPrecision(node);
            double tp = precisionParent + precisionSib;
            preP[thisNumber] = tp * thisPrecision / (tp + thisPrecision);

            //partial means

            for (int j = 0; j < dim; j++) {
                preMeans[thisNumber][j] = (precisionParent * preMeans[parentNumber][j] + precisionSib * postMeans[sibNumber][j]) / (precisionParent + precisionSib);
            }
        }

        if (treeModel.isExternal(node)) {
            return;
        } else {
            doPreOrderTraversal(treeModel.getChild(node, 0));
            doPreOrderTraversal(treeModel.getChild(node, 1));

        }

    }


    public NodeRef getSisterNode(NodeRef node) {
        NodeRef sib0 = treeModel.getChild(treeModel.getParent(node), 0);
        NodeRef sib1 = treeModel.getChild(treeModel.getParent(node), 1);


        if (sib0 == node) {
            return sib1;
        } else return sib0;

    }


    public double sampleNode(NodeRef node) {

        final int thisNumber = node.getNumber();
        double[] traitValue = getNodeTrait(node);


        double[] mean = new double[dim];
        for (int i = 0; i < dim; i++) {
            mean[i] = preMeans[thisNumber][i];
        }

        double p = preP[thisNumber];

        double[][] thisP = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {

                thisP[i][j] = p * precisionParam.getParameterValue(i, j);

            }
        }





 	/*	  Sample it all traits together as one multivariate normal
            *
               MultivariateNormalDistribution distribution  = new MultivariateNormalDistribution(mean, thisP);


 		double[] oldValue = getNodeTrait(node);
 		double[] value = distribution.nextMultivariateNormal();

 		setNodeTrait(node,value);



		double pOld = distribution.logPdf(oldValue[]);

		double pNew = distribution.logPdf(value[]);


 	//	printInformation(oldValue[0]);

 	//	printInformation(value[0]);

  */
        // double[] newTraitValue = getNodeTrait(node);
        //double pNew = distribution.logPdf(newTraitValue);


        ///////////  Individually gibbs sample each entry in the vector


        //  for(int entry=0;entry<dim; entry++){


        int entry = MathUtils.nextInt(dim);

        double thisMean = getConditionalMean(entry, thisP, traitValue, mean);

        double SD = Math.sqrt(1 / thisP[entry][entry]);


        double oldValue = getNodeTrait(node, entry);


        double value = MathUtils.nextGaussian();
        value *= SD;
        value += thisMean;


        //	 printInformation(oldValue);
        //    printInformation(value);


        NormalDistribution distribution = new NormalDistribution(thisMean, SD);


        double pOld = distribution.logPdf(oldValue);


        double pNew = distribution.logPdf(value);


        setNodeTrait(node, entry, value);


        double logq = pOld - pNew;


        traitModel.getTraitParameter().getParameter(thisNumber).fireParameterChangedEvent();


        return logq;


    }



    public double sampleNode2(NodeRef node) {

        final int thisNumber = node.getNumber();



//        double[] mean = new double[dim];
//        for (int i = 0; i < dim; i++) {
//            mean[i] = preMeans[thisNumber][i];
//        }
//
//        double p = preP[thisNumber];
//
//        double[][] thisP = new double[dim][dim];
//
//        for (int i = 0; i < dim; i++) {
//            for (int j = 0; j < dim; j++) {
//
//                thisP[i][j] = p * precisionParam.getParameterValue(i, j);
//
//            }
//        }

        double[] mean=traitModel.getConditionalMean(thisNumber);
        double[][] thisP=traitModel.getConditionalPrecision(thisNumber);


        double[] oldValue = getNodeTrait(node);
        double[] value = oldValue;

        int attempt = 0;
        boolean validTip = false;

        if(hasMask) {

            double[] oldCompVal =new double[numUpdate];
            double[] newValue =new double[numUpdate];
            double[][] condP=new double[numUpdate][numUpdate];

            for (int i=0; i<numUpdate;i++){
                oldCompVal[i]=  value[doUpdate[i]];
                for (int j=0;j<numUpdate;j++){
                    condP[i][j]= thisP[doUpdate[i]][doUpdate[j]];
                }
            }





            double[] condMean = getComponentConditionalMean(thisP, oldValue, mean, condP);



            // Start sampling


            MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(condMean, condP);




            //estourou++;
            while (!validTip & attempt < 10000) {
               newValue = distribution.nextMultivariateNormal();

                 for(int i=0; i<numUpdate; i++){
                   value[doUpdate[i]]=newValue[i];
                 }

                setNodeTrait(node, value);

                if (latentLiability.validTraitForTip(thisNumber)) {
                    validTip = true;
               // estourou--;
                }
                attempt++;
            }           // TODO Failure rate should be stored somewhere and polled later for diagnostics

      /*      if (Math.floorMod(getCount(),10000)==1) {
              printInformation(estourou/10000.0 , "Could not sample truncated normal");
              estourou=0;
              printInformation(((double)getAcceptCount()-ultimoAccept)/10000.0, "Accept probability");
              ultimoAccept=(double)getAcceptCount();

          }

        */


            double pOld = distribution.logPdf(oldCompVal);

            double pNew = distribution.logPdf(newValue);

            double logq= pOld -pNew;


            traitModel.getTraitParameter().getParameter(thisNumber).fireParameterChangedEvent();




            return logq;

              }else {


            //	  Sample it all traits together as one multivariate normal

           MultivariateNormalDistribution distribution = new MultivariateNormalDistribution(mean, thisP);


            while (!validTip & attempt < 10000) {

                value = distribution.nextMultivariateNormal();

                setNodeTrait(node, value);

                if (latentLiability.validTraitForTip(thisNumber)) {
                    validTip = true;
                }
                attempt++;
            }


            double pOld = distribution.logPdf(oldValue);

            double pNew = distribution.logPdf(value);


            double logq = pOld - pNew;


            traitModel.getTraitParameter().getParameter(thisNumber).fireParameterChangedEvent();


            return logq;

        }





    }


   private double[] getComponentConditionalMean(double[][] thisP, double[] oldValue, double[] mean, double[][] condP){

       double[] condMean =new double[numUpdate];

       double[][]  H =new   double[numUpdate][numFixed];
       Matrix  prod  =new Matrix(numUpdate,numFixed);
       Vector dif =new Vector(numUpdate);
       double[] contMeans = new double[numFixed] ;


       for (int i=0; i<numUpdate;i++){
           for (int j=0;j<numFixed;j++){
               H[i][j]= thisP[doUpdate[i]][dontUpdate[j]];
           }
       }

       for (int i=0; i<numFixed;i++){
           contMeans[i]= oldValue[dontUpdate[i]]-mean[dontUpdate[i]];
       }


       Matrix invK= new SymmetricMatrix(condP).inverse() ;
       Matrix HH= new Matrix(H);

       try {
           prod = invK.product(HH);
           dif= prod.product(new Vector(contMeans)) ;

       } catch (IllegalDimension illegalDimension) {
           illegalDimension.printStackTrace();
       }

       for(int i=0; i<numUpdate;i++){
           condMean[i]= mean[doUpdate[i]] - dif.component(i);

       }


       return condMean;
   }






    private double getConditionalMean(int entry, double[][] thisP, double[] traitValue, double[] mean) {
        double sumProd = 0;
        for (int i = 0; i < dim; i++) {
            if (i != entry) sumProd += thisP[entry][i] * (traitValue[i] - mean[i]);
        }

        double condMean = mean[entry] - sumProd / thisP[entry][entry];
        return condMean;
    }


    public String getPerformanceSuggestion() {
        return null;
    }


    public String getOperatorName() {
        return LATENT_LIABILITY_GIBBS_OPERATOR;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public final static String MASK = "mask";


        public String getParserName() {
            return LATENT_LIABILITY_GIBBS_OPERATOR;
        }


        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            if (xo.getChildCount() < 3) {
                throw new XMLParseException(
                        "Element with id = '" + xo.getName() + "' should contain:\n" +
                                "\t 1 conjugate multivariateTraitLikelihood, 1 latentLiabilityLikelihood and one parameter \n"
                );
            }

            double weight = xo.getDoubleAttribute(WEIGHT);

            FullyConjugateMultivariateTraitLikelihood traitModel = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
            LatentTruncation LLModel = (LatentTruncation) xo.getChild(LatentTruncation.class);
            //       MultinomialLatentLiabilityLikelihood LLModel = (MultinomialLatentLiabilityLikelihood) xo.getChild(MultinomialLatentLiabilityLikelihood.class);

            CompoundParameter tipTraitParameter = (CompoundParameter) xo.getChild(CompoundParameter.class);


                Parameter mask =null;



            if (xo.hasChildNamed(MASK)) {
                mask = (Parameter) xo.getElementFirstChild(MASK);
            }






            return new LatentLiabilityGibbs(traitModel, LLModel, tipTraitParameter, mask, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a gibbs sampler on tip latent trais for latent liability model.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(WEIGHT),
                new ElementRule(FullyConjugateMultivariateTraitLikelihood.class, "The model for the latent random variables"),
                new ElementRule(LatentTruncation.class, "The model that links latent and observed variables"),
                new ElementRule(MASK, dr.inference.model.Parameter.class, "Mask: 1 for latent variables that should be sampled",true),
//	                new ElementRule(MultinomialLatentLiabilityLikelihood.class, "The model that links latent and observed variables"),
                new ElementRule(CompoundParameter.class, "The parameter of tip locations from the tree")


        };
    };
}



