/*
 * DriftedTraitsLogger.java
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

package dr.evomodel.antigenic;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.inference.model.*;
import dr.evomodel.continuous.AbstractMultivariateTraitLikelihood;
import dr.xml.*;

/**
 * @author Trevor Bedford
 * @author Andrew Rambaut
 */
public class DriftedTraitsLogger implements TreeTraitProvider {

	public static final String DRIFTED_TRAITS_LOGGER = "driftedTraits";

	private final AbstractMultivariateTraitLikelihood multivariateTraits;
    private final Parameter locationDriftParameter;
    private TreeTrait[] treeTraits = null;
    private double maxHeight = -1.0;

	public DriftedTraitsLogger(AbstractMultivariateTraitLikelihood multivariateTraits, Parameter locationDriftParameter) {
		this.multivariateTraits = multivariateTraits;
        this.locationDriftParameter = locationDriftParameter;
	}

    @Override
    public TreeTrait[] getTreeTraits() {
        if (treeTraits == null) {
            treeTraits = new TreeTrait[] {
                new TreeTrait.DA() {
                    public String getTraitName() {
                        return multivariateTraits.getTraitName();
                    }

                    public Intent getIntent() {
                        return Intent.NODE;
                    }

                    public Class getTraitClass() {
                        return Double.class;
                    }

                    public double[] getTrait(Tree tree, NodeRef node) {
                        double t[] = multivariateTraits.getTraitForNode(tree, node, multivariateTraits.getTraitName());
                        computeMaxHeight(tree);

                        // drift first dimension
                        double nodeHeight = tree.getNodeHeight(node);
                        double offset = locationDriftParameter.getParameterValue(0) * (maxHeight-nodeHeight);
                        t[0] = t[0] + offset;

                        return t;
                    }
                }
            };
        }
        return treeTraits;
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        TreeTrait[] tts = getTreeTraits();
        for (TreeTrait tt : tts) {
            if (tt.getTraitName().equals(key)) {
                return tt;
            }
        }
        return null;
    }

    private void computeMaxHeight(Tree tree) {
        if (maxHeight < 0) {
            int m = tree.getTaxonCount();
            for (int i = 0; i < m; i++) {
                Taxon taxon = tree.getTaxon(i);
                double height = taxon.getHeight();
                if (height > maxHeight) {
                    maxHeight = height;
                }
            }
        }
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return DRIFTED_TRAITS_LOGGER;
		}

		/**
		 * @return an object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            AbstractMultivariateTraitLikelihood multivariateTraits = (AbstractMultivariateTraitLikelihood) xo.getChild(AbstractMultivariateTraitLikelihood.class);
            Parameter locationDrift = (Parameter) xo.getChild(Parameter.class);

		    return new DriftedTraitsLogger(multivariateTraits, locationDrift);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				new ElementRule(AbstractMultivariateTraitLikelihood.class, "The tree trait provider which is to be drifted."),
                new ElementRule(Parameter.class, "The parameter specifying location drift rate.")
		};

		public String getParserDescription() {
			return null;
		}

		public String getExample() {
			return null;
		}

		public Class getReturnType() {
			return TreeTraitProvider.class;
		}
	};

}
