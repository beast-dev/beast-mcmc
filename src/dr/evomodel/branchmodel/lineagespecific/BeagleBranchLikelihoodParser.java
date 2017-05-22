/*
 * BeagleBranchLikelihoodParser.java
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

package dr.evomodel.branchmodel.lineagespecific;

import java.util.ArrayList;
import java.util.List;

import dr.evomodel.tree.TreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BeagleBranchLikelihoodParser extends AbstractXMLObjectParser {

	public static final String BEAGLE_BRANCH_LIKELIHOODS = "beagleBranchLikelihood";

	public static final String UNIQUE_LIKELIHOODS = "uniqueLikelihoods";

	@Override
	public String getParserName() {
		return BEAGLE_BRANCH_LIKELIHOODS;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

		// if(xo.hasChildNamed(TreeModel.TREE_MODEL)) {
		//
		// treeModel = (TreeModel) xo.getChild(TreeModel.class);
		// }

		Parameter zParameter = (Parameter) xo
				.getElementFirstChild(DirichletProcessPriorParser.CATEGORIES);

		List<Likelihood> likelihoods = new ArrayList<Likelihood>();

		XMLObject cxo = (XMLObject) xo.getChild(UNIQUE_LIKELIHOODS);
		for (int i = 0; i < cxo.getChildCount(); i++) {

			Likelihood likelihood = (Likelihood) cxo.getChild(i);
			likelihoods.add(likelihood);
		}

		return null;
//		new BeagleBranchLikelihood(
////				treeModel, likelihoods, zParameter
//				);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public String getParserDescription() {
		return BEAGLE_BRANCH_LIKELIHOODS;
	}

	@Override
	public Class getReturnType() {
		return BeagleBranchLikelihood.class;
	}

}
