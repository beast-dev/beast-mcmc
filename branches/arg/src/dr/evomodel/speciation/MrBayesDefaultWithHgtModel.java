/*
 * YuleModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.ARGModel;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;


/**
 * This class contains methods that describe a Yule speciation with Poisson HGT event model.
 *
 * @author Marc Suchard
 */

public class MrBayesDefaultWithHgtModel extends SpeciationModel {

	public static final String MRBAYES_DEFAULT_MODEL = "mrbayesDefaultModel";
	public static String BIRTH_RATE = "birthRate";
	public static final String HGT_RATE = "hgtRate";


	public MrBayesDefaultWithHgtModel(Parameter birthRateParameter, Parameter hgtRateParameter, int units) {

		super(MRBAYES_DEFAULT_MODEL, units);

		this.birthRateParameter = birthRateParameter;
		this.hgtRateParameter = hgtRateParameter;
		addParameter(birthRateParameter);
		addParameter(hgtRateParameter);
		birthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
		hgtRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

	}

	public double getBirthRate() {
		return birthRateParameter.getParameterValue(0);
	}

	public double getHgtRate() {
		return hgtRateParameter.getParameterValue(0);
	}

	public void setBirthRate(double birthRate) {

		birthRateParameter.setParameterValue(0, birthRate);
	}

	//
	// functions that define a speciation model
	//
	public double logTreeProbability(int taxonCount) {
		if (true)
			throw new RuntimeException("Why was 'logTreeProbability' called?");
		return 0.0;
	}

	private static double logGamma(int i) {
		return logFactorial(i - 1);
	}

	private static double logFactorial(int i) {
		double result = 0;
		if (i > 1) {
			for (int j = 2; j <= i; j++)
				result += Math.log(j);
		}
		return result;
	}

	public double logReassortmentProbability(ARGModel arg) {

		int count = arg.getReassortmentNodeCount();
//        System.err.println(count);
		double lambda = getHgtRate();

		if (lambda == 0) {
//            System.err.println("yo");
			if (count > 1)
				return Double.NEGATIVE_INFINITY;
			else
				return 0.0;
		}

		/*     double logFactorial = 0.0;
		  if (count > 1) {
			  for(int j=2; j<=count; j++)
				  logFactorial += Math.log(j);
		  }

		  System.err.printf("logGamma(%d) = %5.4f\n",count,logGamma(count));
  */
//	    System.err.println(lambda);
		
		return -lambda + count * Math.log(lambda) - logFactorial(count);
				
	}


	//
	// functions that define a speciation model
	//
	public double logNodeProbability(Tree tree, NodeRef node) {
		if (tree.getRoot() != node)
			//  tree.set
			return 0.0;

//        ARGModel arg = (ARGModel)tree;

//        if (arg.isReassortment(node))
//            return 0.0;

//        double nodeHeight = tree.getNodeHeight(node);
		double rootHeight = tree.getNodeHeight(node);

		double lambda = getBirthRate();
//        System.err.println(getBirthRate());
		//return Math.log((lambda * Math.exp(-lambda * nodeHeight)) / (1 - Math.exp(-lambda * rootHeight)));
		int randomHeights = tree.getInternalNodeCount() - 1; // do no count root
//	    System.err.println("rhc = "+randomHeights);
//        return -lambda * rootHeight - logGamma(randomHeights);
		return 0;
	}

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

	public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
		throw new RuntimeException("createElement not implemented");
	}

	/**
	 * Parses an element from an DOM document into a SpeciationModel. Recognises
	 * YuleModel.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MRBAYES_DEFAULT_MODEL;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			//int units = XMLParser.Utils.getUnitsAttr(xo);
			int units = XMLParser.Utils.getUnitsAttr(xo);

			XMLObject cxo = (XMLObject) xo.getChild(BIRTH_RATE);
			Parameter brParameter = (Parameter) cxo.getChild(Parameter.class);
			cxo = (XMLObject) xo.getChild(HGT_RATE);
			Parameter hgtParameter = (Parameter) cxo.getChild(Parameter.class);

			return new MrBayesDefaultWithHgtModel(brParameter, hgtParameter, units);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A speciation model of a Yule with Horizontal Gene Transfer process.";
		}

		public Class getReturnType() {
			return YuleModel.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(BIRTH_RATE,
						new XMLSyntaxRule[]{new ElementRule(Parameter.class), new ElementRule(Parameter.class)}),
				XMLUnits.SYNTAX_RULES[0]
		};
	};


	//Protected stuff
	private Parameter birthRateParameter;
	private Parameter hgtRateParameter;
}
