/*
 * AlloppSpeciesBindingsParser.java
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

import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.tree.TreeModel;
import dr.xml.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Parser for list of possibly-allopolyploid species containing individuals
 * in turn containing taxons (which are diploid genomes).
 *
 * @author Graham Jones
 *         Date: 18/04/2011
 */




/*
 *
 * Parses something like this:
 
<alloppspecies id="alloppspecies">
  <sp id="Alpha" ploidylevel = 2>
    <individual id = "1">  
      <taxon idref="1_Alpha" /> 
    </individual>
    <individual id = "2">  
      <taxon idref="2_Alpha" /> 
    </individual>
    <individual id = "3">  
      <taxon idref="3_Alpha" /> 
    </individual>
  </sp>
  <sp id="Beta" ploidylevel = 4>
    <individual id = "4">  
      <taxon idref="4_Beta_A" /> 
      <taxon idref="4_Beta_B" />
    </individual>
    <individual id = "5">  
      <taxon idref="5_Beta_A" /> 
      <taxon idref="5_Beta_B" />
    </individual>
  </sp>
   ...
   (more species, then genetrees)
</alloppspecies>

 */


// I have adapted code from SpeciesBindingsParser.
// I changed 'ploidy' to 'popfactor' to reduce confusion with the ploidy level
// of a species which is independent of which gene is considered.
// Use of popfactors is untested. Maybe chloroplast data would use it.


public class AlloppSpeciesBindingsParser extends AbstractXMLObjectParser {
    public static final String ALLOPPSPECIES = "alloppspecies";
    public static final String GENE_TREES = "geneTrees";
    public static final String GTREE = "gtree";
    public static final String POPFACTOR = "popfactor";
    public static final String MIN_GENENODE_HEIGHT = "minGeneNodeHeight";


    public String getParserName() {
	    return ALLOPPSPECIES;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		
		List<AlloppSpeciesBindings.ApSpInfo> apsp = new ArrayList<AlloppSpeciesBindings.ApSpInfo>();
        for (int k = 0; k < xo.getChildCount(); ++k) {
            final Object child = xo.getChild(k);
            if (child instanceof AlloppSpeciesBindings.ApSpInfo) {
                apsp.add((AlloppSpeciesBindings.ApSpInfo) child);
            }
        }
        final double mingenenodeheight = xo.getDoubleAttribute(MIN_GENENODE_HEIGHT);
        final XMLObject xogt = xo.getChild(GENE_TREES);
        final int nTrees = xogt.getChildCount();
        final TreeModel[] trees = new TreeModel[nTrees];
        double[] popFactors = new double[nTrees];

        for (int nt = 0; nt < trees.length; ++nt) {
            Object child = xogt.getChild(nt);
            if (!(child instanceof TreeModel)) {
                assert child instanceof XMLObject;
                popFactors[nt] = ((XMLObject) child).getDoubleAttribute(POPFACTOR);
                child = ((XMLObject) child).getChild(TreeModel.class);

            } else {
                popFactors[nt] = -1;
            }
            trees[nt] = (TreeModel) child;
        }

        try {
            return new AlloppSpeciesBindings(apsp.toArray(new AlloppSpeciesBindings.ApSpInfo[apsp.size()]),
            		                         trees, mingenenodeheight, popFactors);
        } catch (Error e) {
            throw new XMLParseException(e.getMessage());
        }
	}

	
	// I have adapted code from SpeciesBindingsParser 
	// I changed 'Ploidy' to 'PopFactors' to reduce confusion -
	// the only use I can think of for popfactors in an AlloppNetwork is chloroplast data 
    ElementRule treeWithPopFactors = new ElementRule(GTREE,
            new XMLSyntaxRule[]{AttributeRule.newDoubleRule(POPFACTOR),
                    new ElementRule(TreeModel.class)}, 0, Integer.MAX_VALUE);

    
	@Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
        		AttributeRule.newDoubleRule(MIN_GENENODE_HEIGHT),
                new ElementRule(AlloppSpeciesBindings.ApSpInfo.class, 2, Integer.MAX_VALUE),
                new ElementRule(GENE_TREES,
                        new XMLSyntaxRule[]{
                                new ElementRule(TreeModel.class, 0, Integer.MAX_VALUE),
                                treeWithPopFactors
                        }),
        };
    }	

	
	@Override
	public String getParserDescription() {
        return "Binds taxa to gene trees with information about possibly allopolyploid species.";
	}

	@Override
	public Class getReturnType() {
		return AlloppSpeciesBindings.class;
	}

}
