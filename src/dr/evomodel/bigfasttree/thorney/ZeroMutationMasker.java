/*
 * ZeroMutationMasker.java
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

package dr.evomodel.bigfasttree.thorney;

import dr.evomodel.treedatalikelihood.discrete.MaskProvider;
import dr.inference.model.Parameter;

import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 * @author JT McCrone
 */
public class ZeroMutationMasker implements  MaskProvider {

    private final ConstrainedTreeModel tree;
    private final Parameter mask;
    private final MutationBranchMap branchLengthProvider;


    public ZeroMutationMasker(ConstrainedTreeModel tree,
                              Parameter mask,
                              MutationBranchMap branchLengthProvider
                            ) {
        this.tree = tree;
        this.mask = mask;
        this.branchLengthProvider = branchLengthProvider;
        updateMask();
    }
    @Override
    public Parameter getMask() {
        return mask;
    }

    @Override
    public void updateMask() {
        for (int i = 0; i < mask.getDimension(); i++) {
            final double currentMaskValue = branchLengthProvider.getMutations(tree.getInternalNode(i)).getMutationCount() >0?1.0:0.0;
            mask.setParameterValueQuietly(i, currentMaskValue);
        }
        mask.fireParameterChangedEvent();
    }

// **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String MASKER = "ZeroMutationMasker";
        final static String MASK = "mask";

        public String getParserName() {
            return MASKER;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter mask = (Parameter) xo.getChild(MASK).getChild(Parameter.class);
            ConstrainedTreeModel treeModel = (ConstrainedTreeModel) xo.getChild(ConstrainedTreeModel.class);
            MutationBranchMap branchLengthProvider = (MutationBranchMap) xo.getChild(MutationBranchMap.class);
            if (mask.getDimension() == 1) {
                mask.setDimension(treeModel.getInternalNodeCount());
            }
            return new ZeroMutationMasker(treeModel, mask,branchLengthProvider);
        }

        public String getParserDescription() {
            return "A utility to craft mask for filtering dimensions in where nodes do not have any mutations";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MASK, Parameter.class),
                new ElementRule(ConstrainedTreeModel.class),
                new ElementRule(MutationBranchMap.class),

        };

        public Class getReturnType() {
            return ZeroMutationMasker.class;
        }
    };


}
