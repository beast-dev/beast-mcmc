/*
 * DriftedLocationLogger.java
 *
 * Copyright (C) 2002-2013 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.continuous.IntegratedMultivariateTraitLikelihood;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 */
public class DriftedLocationLogger implements TreeTraitProvider {

	public static final String DRIFTED_TRAIT_LOGGER = "driftedTraits";

	private final TreeTraitProvider treeTraitProvider;

	public DriftedLocationLogger(TreeTraitProvider treeTraitProvider) {
		this.treeTraitProvider = treeTraitProvider;
	}

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraitProvider.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraitProvider.getTreeTrait(key);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return DRIFTED_TRAIT_LOGGER;
		}

		/**
		 * @return an object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeTraitProvider treeTraitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

		    return new DriftedLocationLogger(treeTraitProvider);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				new ElementRule(TreeTraitProvider.class, "The tree trait provider which is to be drifted")
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
