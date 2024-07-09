/*
 * AnnotateMuTreeTrait.java
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

package dr.evomodel.antigenic.phyloclustering.statistics;


import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


public class AnnotateMuTreeTrait implements TreeTraitProvider {
	
    protected Helper treeTraits = new Helper();
    private TreeModel treeModel;
    public static final String ANNOTATEMUPARAMETERTREETRAIT = "annotateMuTreeTrait";
    public final static String muParameterStr = "mu";
    //public final static String serumDriftStr = "serumDrift";  //but maybe serum drift is off
    public final static String indicatorsStr = "indicators";
    public final static String hasDriftStr = "hasDrift";

    public final static String MU1_SCALE_PARAMETER = "mu1Scale";
    public final static String MU2_SCALE_PARAMETER = "mu2Scale";
    
    public final static String MUMEAN_PARAMETER = "muMean";
    
    
    private MatrixParameter muParameter;
    private Parameter indicators;
   // private Parameter serumDrift;
    private Parameter mu1ScaleParameter;
    private Parameter mu2ScaleParameter;
    private Parameter muMeanParameter;
    private boolean hasDrift;

    
    
	
	// public AnnotateMuTreeTrait(TreeModel treeModel_in, MatrixParameter muTreeNode_in, Parameter indicators_in, Parameter serumDrift_in, boolean driftTreeClusterMu){
    public AnnotateMuTreeTrait(TreeModel treeModel_in, MatrixParameter muTreeNode_in, Parameter indicators_in, boolean driftTreeClusterMu, Parameter mu1Scale, Parameter mu2Scale, Parameter muMean) {	
    
		 this.treeModel = treeModel_in;
		 this.muParameter = muTreeNode_in;
		 this.indicators = indicators_in;
		 //this.serumDrift = serumDrift_in;
		 this.mu1ScaleParameter = mu1Scale;
		 this.mu2ScaleParameter = mu2Scale;
		 this.muMeanParameter = muMean;
		 this.hasDrift = driftTreeClusterMu;
		 
	        treeTraits.addTrait(new TreeTrait.IA() {
	        	
	            public String getTraitName() {
	            	//System.out.println("print label");
	               // return tag;
	            	return "mu";
	            }

	            public String getTraitString(Tree tree, NodeRef node) {
	            	
	            	if(tree != treeModel){
	            		System.out.println("Something is wrong. Why is tree not equal to treeModel?");
	            		System.exit(0);
	            	}
	            	
	            	String outputStr = "{0,0}";
	            	if( (int)indicators.getParameterValue(node.getNumber() )  ==1 ){
	            		double firstCoord = muParameter.getParameter(node.getNumber()).getParameterValue(0);
	            		double secondCoord =  muParameter.getParameter(node.getNumber()).getParameterValue(1);
	            	
	            	
	            		//if(hasDrift && serumDrift != null ){
	            		//	//multiply the first coordinate by the drift term
	            		//	firstCoord = firstCoord * serumDrift.getParameterValue(0);
	            		//}
	               		if(hasDrift && mu1ScaleParameter != null && muMeanParameter != null ){
	            			//multiply the first coordinate by the drift term
	            			firstCoord =  firstCoord * mu1ScaleParameter.getParameterValue(0);
	            		}
	               		if(hasDrift && mu2ScaleParameter != null ){
	            			//multiply the first coordinate by the drift term
	            			secondCoord = secondCoord * mu2ScaleParameter.getParameterValue(0);
	            		}
	            		
	            		outputStr = "{"+ firstCoord + "," + secondCoord +"}";
	            	}
	            	return outputStr;
	            }
	            
	            
	            public Intent getIntent() {
	            	//System.out.println("getIntent");
	                return Intent.NODE;
	            }

	            public Class getTraitClass() {
	            	System.out.println("getTraitClass ran. Not expected. Quit now");
	            	System.exit(0);
	                return int[].class;
	            }

	            
	            public int[] getTrait(Tree tree, NodeRef node) {
	              //  return getStatesForNode(tree, node);
	            	System.out.println("getTrait ran. Not expected. Quit now");
	            	System.exit(0);
	            	//int x[] = new int[10];
	            	return null;
	            }


	        });

		 
	 }
	 
	 

	    public TreeTrait[] getTreeTraits() {
	    	//System.out.println("hihi");
	        return treeTraits.getTreeTraits();
	    }
	
	
	    public TreeTrait getTreeTrait(String key) {
	    	System.out.println("not expected to run getTreeTrait. Quit now");
	    	System.exit(0);
	        return treeTraits.getTreeTrait(key);
	    }
	 

	    

	    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {


	        public String getParserName() {
	            return ANNOTATEMUPARAMETERTREETRAIT;
	        }

	        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
	            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

	            XMLObject cxo = xo.getChild(muParameterStr);
                MatrixParameter muParameter = (MatrixParameter) cxo.getChild(MatrixParameter.class);
           
                 cxo = xo.getChild(indicatorsStr);
                Parameter indicators = (Parameter) cxo.getChild(Parameter.class);
                
              //   cxo = xo.getChild(serumDriftStr);
               // Parameter serumDrift = (Parameter) cxo.getChild(Parameter.class);
           
                boolean hasDrift = xo.getAttribute(hasDriftStr, false);

                Parameter mu1Scale = null;
                if (xo.hasChildNamed(MU1_SCALE_PARAMETER)) {
                	mu1Scale = (Parameter) xo.getElementFirstChild(MU1_SCALE_PARAMETER);
                }
                
                Parameter mu2Scale = null;
                if (xo.hasChildNamed(MU2_SCALE_PARAMETER)) {
                	mu2Scale = (Parameter) xo.getElementFirstChild(MU2_SCALE_PARAMETER);
                }
                
                Parameter muMean = null;
                if(xo.hasChildNamed(MUMEAN_PARAMETER)){
                	muMean = (Parameter) xo.getElementFirstChild(MUMEAN_PARAMETER);
                }
	        	// return new AnnotateMuTreeTrait( treeModel, clusterLabelsTreeNode, indicators, serumDrift,driftTreeClusterMu);
                return new AnnotateMuTreeTrait( treeModel, muParameter, indicators, hasDrift, mu1Scale, mu2Scale, muMean);

	        }

	        //************************************************************************
	        // AbstractXMLObjectParser implementation
	        //************************************************************************

	        public String getParserDescription() {
	            return "Integrate mu coordinates into the tree.";
	        }

	        public Class getReturnType() {
	            return AnnotateMuTreeTrait.class;
	        }

	        public XMLSyntaxRule[] getSyntaxRules() {
	            return rules;
	        }

	        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
	            new ElementRule(TreeModel.class),
                new ElementRule(muParameterStr, MatrixParameter.class),
                new ElementRule(indicatorsStr,Parameter.class),
           //     new ElementRule(serumDriftStr, Parameter.class),
                new ElementRule(MU1_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the first dimension of mu"),
                new ElementRule(MU2_SCALE_PARAMETER, Parameter.class, "Optional parameter for scaling the second dimension of mu"),
                new ElementRule(MUMEAN_PARAMETER, Parameter.class),
                AttributeRule.newBooleanRule(hasDriftStr, true, "whether to multiple the mu by the drift term"),

	        };
	    };
	
	
}
