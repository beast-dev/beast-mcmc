/*
 * RandomBranchModelParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.branchmodel;

import java.util.logging.Logger;

import dr.evomodel.branchmodel.RandomBranchModel;
import dr.evomodel.substmodel.codon.GY94CodonModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 * 
 */

public class RandomBranchModelParser extends AbstractXMLObjectParser {

	 public static final String BASE_MODEL = "baseSubstitutionModel";
	 public static final String SEED = "seed";
	 
	 public static final String RATE = "rate";
	 
	@Override
	public String getParserName() {
		return RandomBranchModel.RANDOM_BRANCH_MODEL;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Logger.getLogger("dr.evomodel").info("\nUsing random assignment branch model.");
        
        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        GY94CodonModel baseSubstitutionModel = (GY94CodonModel) xo .getElementFirstChild(BASE_MODEL);
        
        
        long seed = -1;
        boolean hasSeed = false;
        if (xo.hasAttribute(SEED)) {
        	seed = xo.getLongIntegerAttribute(SEED);
            hasSeed = true;
        }
        
        
        double rate = 1;
        if (xo.hasAttribute(RATE)) {
        	rate = xo.getDoubleAttribute(RATE);
        }
        
        
		return new RandomBranchModel(treeModel, baseSubstitutionModel, rate, hasSeed, seed);
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {

		return new XMLSyntaxRule[] {

		AttributeRule.newDoubleRule(RATE, true, "Rate of the exponentially distributed random component."),		//
		new ElementRule(TreeModel.class, false), //
        new ElementRule(BASE_MODEL,
                new XMLSyntaxRule[]{ new ElementRule(SubstitutionModel.class, 1, 1) })

		};
	}// END: XMLSyntaxRule

	@Override
	public String getParserDescription() {
		return RandomBranchModel.RANDOM_BRANCH_MODEL;
	}

	@Override
	public Class getReturnType() {
		return RandomBranchModel.class;
	}


}//END: class
