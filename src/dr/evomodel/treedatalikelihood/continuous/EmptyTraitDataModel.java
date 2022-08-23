/*
 * ContinuousTraitData.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.xml.*;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class EmptyTraitDataModel implements ContinuousTraitPartialsProvider {

    private final String name;
    private final int dimTrait;
    private final CompoundParameter traitParameter;
    private final PrecisionType precisionType;

    public EmptyTraitDataModel(String name, CompoundParameter traitParameter,
                               final int dimTrait, PrecisionType precisionType) {
        this.name = name;
        this.traitParameter = traitParameter;
        this.dimTrait = dimTrait;
        this.precisionType = precisionType;
    }

    @Override
    public boolean bufferTips() {
        return true;
    } // TODO maybe should be false

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getTraitDimension() {
        return dimTrait;
    }

    @Override
    public String getTipTraitName() {
        return null;
    }

    @Override
    public void setTipTraitName(String name) {
        // do nothing
    }

    @Override
    public PrecisionType getPrecisionType() {
        return precisionType;
    }

    @Override
    public CompoundParameter getParameter() {
        return traitParameter;
    }

    @Override
    public String getModelName() {
        return name;
    }

    @Override
    public boolean usesMissingIndices() {
        return false;
    }

    @Override
    public ContinuousTraitPartialsProvider[] getChildModels() {
        return new ContinuousTraitPartialsProvider[0];
    }

    @Override
    public List<Integer> getMissingIndices() {
        return null;
    }

    @Override
    public boolean[] getDataMissingIndicators() {
        return null;
    }

    @Override
    public double[] getTipPartial(int taxonIndex, boolean fullyObserved) {
        return new double[precisionType.getPartialsDimension(dimTrait)];
    }
//
//    private static final String EMPTY_TRAIT_MODEL = "emptyTraitModel";
//
//    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
//        @Override
//        public Object parseXMLObject(XMLObject xo) throws XMLParseException{
//            MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(TreeModel.class);
//            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();
//
//            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
//                    utilities.parseTraitsFromTaxonAttributes(xo, TreeTraitParserUtilities.DEFAULT_TRAIT_NAME,
//                            treeModel, true);
//
//            String traitName = returnValue.traitName;
//            MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel)
//                    xo.getChild(MultivariateDiffusionModel.class);
//
//            PrecisionType precision = PrecisionType.SCALAR;
//
//            return new EmptyTraitDataModel(
//                    traitName,
//                    diffusionModel.getPrecisionParameter().getRowDimension(),
//                    precision
//
//            );
//        }
//
//        @Override
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return new XMLSyntaxRule[0];
//        }
//
//        @Override
//        public String getParserDescription() {
//            return null;
//        }
//
//        @Override
//        public Class getReturnType() {
//            return null;
//        }
//
//        @Override
//        public String getParserName() {
//            return null;
//        }
//    };
}

