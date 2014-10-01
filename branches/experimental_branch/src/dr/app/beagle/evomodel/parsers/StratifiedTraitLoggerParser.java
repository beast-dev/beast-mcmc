/*
 * StratifiedTraitLoggerParser.java
 *
 * Copyright (C) 2002-2012 Alexei Drummond, Andrew Rambaut & Marc A. Suchard
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

package dr.app.beagle.evomodel.parsers;

import dr.app.beagle.evomodel.utilities.TreeTraitLogger;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.xml.*;
import dr.app.beagle.evomodel.substmodel.CodonPartitionedRobustCounting;
import dr.app.beagle.evomodel.substmodel.StratifiedTraitOutputFormat;
import dr.evomodel.tree.TreeModel;

/**
 * @author Marc A. Suchard
 */

public class StratifiedTraitLoggerParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "stratifiedTraitLogger";
    public static final String TRAIT_NAME = "traitName";
    public static final String PARTITION = "partition";
    //public static final String LOG_FORMAT = "format";

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);
        TreeTraitProvider traitProvider =
                (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

        String traitName = (String) xo.getAttribute(TRAIT_NAME);

        TreeTrait trait = traitProvider.getTreeTrait(traitName);
        if (trait == null || trait.getIntent() != TreeTrait.Intent.WHOLE_TREE) {
            StringBuilder sb = new StringBuilder("Unable to find whole tree trait '" + traitName + "' in '" + xo.getId() + "\n");
            sb.append("\tPossible traits:");
            for (TreeTrait existingTrait : traitProvider.getTreeTraits()) {
                if (existingTrait.getIntent() == TreeTrait.Intent.WHOLE_TREE) {
                    sb.append(" " +existingTrait.getTraitName());
                }
            }
            sb.append("\n");
            throw new XMLParseException(sb.toString());
        }

        boolean partition = xo.getAttribute(PARTITION, false);

        if (trait.getTraitClass() == Double.class || trait.getTraitClass() == Integer.class || !partition) {
            return new TreeTraitLogger( treeModel, new TreeTrait[] { trait });
        }
                                                    
        int length;
        Object obj = trait.getTrait(treeModel, treeModel.getNode(0));

        if (obj instanceof double[]) {
           length = ((double[])obj).length;
        } else if (obj instanceof int[]) {
            length = ((int[])obj).length;
        } else {
            throw new XMLParseException("Unknown trait type for partitioning");
        }
       
        TreeTrait[] partitionedTraits = new TreeTrait[length];
        for (int i = 0; i < length; i++) {
            partitionedTraits[i] = new TreeTrait.PickEntryD(trait, i);
        }
        return new TreeTraitLogger( treeModel, partitionedTraits );        
    }

  private static final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
          AttributeRule.newStringRule(TRAIT_NAME),
          AttributeRule.newBooleanRule(PARTITION, true),
//          AttributeRule.newStringRule(LOG_FORMAT, true),
          new ElementRule(TreeModel.class),
          new ElementRule(TreeTraitProvider.class),
    };

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    public String getParserDescription() {
        return "A parser to stratify traits to be logged";
    }

    public Class getReturnType() {
        return TreeTraitLogger.class;
    }

    public String getParserName() {
        return PARSER_NAME;
    }
}
