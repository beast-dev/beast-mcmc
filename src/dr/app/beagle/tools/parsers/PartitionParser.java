/*
 * PartitionParser.java
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

package dr.app.beagle.tools.parsers;

import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.app.beagle.tools.Partition;
import dr.evolution.sequence.Sequence;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.tree.TreeModel;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;
import dr.xml.XORRule;

/**
 * @author Filip Bielejec
 * @version $Id$
 */
public class PartitionParser extends AbstractXMLObjectParser {

	public static final String PARTITION = "partition";
    public static final String FROM = "from";
    public static final String TO = "to";
    public static final String EVERY = "every";
	
	@Override
	public String getParserName() {
		return PARTITION;
	}

	@Override
	public String getParserDescription() {
		return "Partition element";
	}

	@Override
	public Class<Partition> getReturnType() {
		return Partition.class;
	}
	
	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
        int from = 0;
        int to = -1;
        int every = xo.getAttribute(EVERY, 1);
		
		if (xo.hasAttribute(FROM)) {
			from = xo.getIntegerAttribute(FROM) - 1;
			
			if (from < 0) {
				throw new XMLParseException(
						"Illegal 'from' attribute in patterns element");
			}
			
		}// END: from check

		if (xo.hasAttribute(TO)) {
			
			to = xo.getIntegerAttribute(TO) - 1;
			if (to < 0 || to < from) {
				throw new XMLParseException(
						"Illegal 'to' attribute in patterns element");
			}
			
		}// END: to check

		if (every <= 0) {
			throw new XMLParseException(
					"Illegal 'every' attribute in patterns element");
		}// END: every check

		TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);
		
		GammaSiteRateModel siteModel = (GammaSiteRateModel) xo.getChild(GammaSiteRateModel.class);
		FrequencyModel freqModel = (FrequencyModel) xo.getChild(FrequencyModel.class);
		Sequence rootSequence = (Sequence) xo.getChild(Sequence.class);
		
		BranchRateModel rateModel = (BranchRateModel) xo.getChild(BranchRateModel.class);
		if (rateModel == null) {
			rateModel = new DefaultBranchRateModel();
		} 
		
		BranchModel branchModel = (BranchModel) xo.getChild(BranchModel.class);
		if (branchModel == null) {

			SubstitutionModel substitutionModel = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
			branchModel = new HomogeneousBranchModel(substitutionModel);
		
		}
		
		Partition partition = new Partition(tree, branchModel, siteModel, rateModel, freqModel, from, to, every);

		if (rootSequence != null) {
				partition.setRootSequence(rootSequence);
		}// END: ancestralSequence check

		return partition;
	}//END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		
		return new XMLSyntaxRule[] {
				AttributeRule.newIntegerRule(FROM, true,
								"The site position to start at, default is 1 (the first position)"), //
				AttributeRule.newIntegerRule(TO, true,
						"The site position to finish at, must be greater than <b>"
								+ FROM
								+ "</b>, default is number of sites"), //
				AttributeRule.newIntegerRule(
								EVERY,
								true,
								"Determines how many sites are selected. A value of 3 will select every third site starting from <b>"
										+ FROM
										+ "</b>, default is 1 (every site)"), //
				new ElementRule(TreeModel.class), //
				new XORRule(new ElementRule(BranchModel.class),
						new ElementRule(SubstitutionModel.class), false), //
				new ElementRule(GammaSiteRateModel.class), //
				new ElementRule(BranchRateModel.class, true), //
				new ElementRule(FrequencyModel.class), //
				new ElementRule(Sequence.class, true) //
		};
		
	}// END: getSyntaxRules
	
}// END: class
