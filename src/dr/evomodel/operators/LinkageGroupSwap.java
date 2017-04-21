/*
 * LinkageGroupSwap.java
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

package dr.evomodel.operators;

import dr.evolution.alignment.Alignment;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.HiddenLinkageModel;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * @author Aaron Darling
 */
public class LinkageGroupSwap extends SimpleMCMCOperator {

	HiddenLinkageModel hlm;	
	int groupCount;
	int columnCount;

	public LinkageGroupSwap(HiddenLinkageModel hlm, double weight){
		this.hlm = hlm;
		groupCount = hlm.getLinkageGroupCount();
		columnCount = hlm.getData().getAlignment().getSiteCount();
        setWeight(weight);
	}

	@Override
	public double doOperation() {
		if(MathUtils.nextBoolean()){
			// swap all taxa in one group to a new group
			int A = MathUtils.nextInt(groupCount);
			int B = MathUtils.nextInt(groupCount);
			HashSet<Taxon> aTaxa = new HashSet<Taxon>(hlm.getGroup(A));
			HashSet<Taxon> bTaxa = new HashSet<Taxon>(hlm.getGroup(B));
			for(Taxon taxon : aTaxa){
				hlm.moveReadGroup(taxon, A, B);
			}
			for(Taxon taxon : bTaxa){
				hlm.moveReadGroup(taxon, B, A);
			}
		}else{
			// pick two linkage groups uniformly A, B
			// pick an alignment column uniformly X
			// move all reads in group A that are defined in column X to group B, and
			// vice versa from B to A.
			// pick a read uniformly at random, add it to a linkage group uniformly at random
			ArrayList<Taxon> aTaxa = new ArrayList<Taxon>();
			ArrayList<Taxon> bTaxa = new ArrayList<Taxon>();
			int A=0, B=0;
	
			// iterate until we actually find something to move
			// this could be done more efficiently, e.g. by limiting
			// choice of X and Y to groups that have something in col X.
			while(aTaxa.size()==0&&bTaxa.size()==0)
			{
				int X = MathUtils.nextInt(columnCount);
		
				A = MathUtils.nextInt(groupCount);
				B = MathUtils.nextInt(groupCount);
				if(A==B)
					continue;	// nothing to do.
		
				// find all reads intersecting column X
				Alignment aln = hlm.getData().getAlignment();
				for(int i=0; i<aln.getTaxonCount(); i++){
					int state = aln.getPatternState(i, X);
					if(state==hlm.getDataType().getGapState() || state==hlm.getDataType().getUnknownState())
						continue;	// seq undefined in this column
					if(hlm.getGroup(A).contains(aln.getTaxon(i))){
						aTaxa.add(aln.getTaxon(i));
					}
					if(hlm.getGroup(B).contains(aln.getTaxon(i))){
						bTaxa.add(aln.getTaxon(i));
					}
				}
			}
	
			// move taxa from A to B and from B to A
			for(Taxon taxon : aTaxa){
				hlm.moveReadGroup(taxon, A, B);
			}
			for(Taxon taxon : bTaxa){
				hlm.moveReadGroup(taxon, B, A);
			}
		}
		return 0;
	}

	public String getOperatorName() {
		return "linkageGroupSwap";
	}

	public String getPerformanceSuggestion() {
		return "Ask Aaron Darling to write a better operator";
	}

    public static XMLObjectParser LINKAGE_GROUP_SWAP_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return "linkageGroupSwap";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        	HiddenLinkageModel hlm = (HiddenLinkageModel)xo.getChild(HiddenLinkageModel.class);
            double weight = xo.getDoubleAttribute("weight");

            return new LinkageGroupSwap(hlm, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents an operator that swaps taxa among two linkage groups. ";
        }

        public Class getReturnType() {
            return LinkageGroupSwap.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(HiddenLinkageModel.class)
        };

    };
}
