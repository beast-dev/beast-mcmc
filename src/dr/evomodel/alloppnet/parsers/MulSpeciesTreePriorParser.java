/*
 * MulSpeciesTreePriorParser.java
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

package dr.evomodel.alloppnet.parsers;


import dr.evomodel.alloppnet.speciation.MulSpeciesTreeModel;
import dr.evomodel.alloppnet.speciation.MulSpeciesTreePrior;
import dr.evomodel.speciation.SpeciationModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class MulSpeciesTreePriorParser extends AbstractXMLObjectParser {
	public static final String MUL_SPECIES_TREE_PRIOR = "mulSpeciesTreePrior";
	public static final String MODEL = "model";
	public static final String MUL_SPECIES_TREE = "mulTree";


	public String getParserName() {
		return MUL_SPECIES_TREE_PRIOR;
	}


	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final XMLObject mxo = xo.getChild(MODEL);
		final SpeciationModel sppm = (SpeciationModel) mxo.getChild(SpeciationModel.class);
		final XMLObject mulsptxo = xo.getChild(MUL_SPECIES_TREE);
		final MulSpeciesTreeModel mulspt = (MulSpeciesTreeModel) mulsptxo.getChild(MulSpeciesTreeModel.class);
		return new MulSpeciesTreePrior(sppm, mulspt);	
	}

	private  XMLSyntaxRule[] modelRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(SpeciationModel.class)
        };
    }

    private  XMLSyntaxRule[] mulsptRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(MulSpeciesTreeModel.class)
        };
    }
	
	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[]{
				new ElementRule(MODEL, modelRules()),
				new ElementRule(MUL_SPECIES_TREE, mulsptRules()),
				
		};
	}	
	


	@Override
	public String getParserDescription() {
		return "Prior for a multiply-labelled species tree for allopolyploids.";
	}

	@Override
	public Class getReturnType() {
		return MulSpeciesTreePrior.class;
	}

}
