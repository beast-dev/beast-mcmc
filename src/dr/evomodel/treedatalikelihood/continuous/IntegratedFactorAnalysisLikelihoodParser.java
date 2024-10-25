/*
 * IntegratedFactorAnalysisLikelihoodParser.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.app.beauti.components.BeautiModelIDProvider;
import dr.app.beauti.components.BeautiParameterIDProvider;
import dr.evolution.tree.MutableTreeModel;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.util.TaskPool;
import dr.xml.*;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.STANDARDIZE;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TARGET_SD;

public class IntegratedFactorAnalysisLikelihoodParser extends AbstractXMLObjectParser implements BeautiModelIDProvider {

    public static final String INTEGRATED_FACTOR_Model = "integratedFactorModel";
    public static final String LOADINGS = "loadings";
    public static final String PRECISION = "precision";
    private static final String NUGGET = "nugget";
    private static final String CACHE_PRECISION = "cachePrecision";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MutableTreeModel treeModel = (MutableTreeModel) xo.getChild(MutableTreeModel.class);
        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();

        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo,
                        treeModel, true);
        CompoundParameter traitParameter = returnValue.traitParameter;
        boolean[] missingIndicators = returnValue.getMissingIndicators();

        MatrixParameterInterface loadings = (MatrixParameterInterface) xo.getElementFirstChild(LOADINGS);
        Parameter traitPrecision = (Parameter) xo.getElementFirstChild(PRECISION);

        double nugget = xo.getAttribute(NUGGET, 0.0);

        TaskPool taskPool = (TaskPool) xo.getChild(TaskPool.class);

        IntegratedFactorAnalysisLikelihood.CacheProvider cacheProvider;
        boolean useCache = xo.getAttribute(CACHE_PRECISION, false);
        if (useCache) {
            cacheProvider = IntegratedFactorAnalysisLikelihood.CacheProvider.USE_CACHE;
        } else {
            cacheProvider = IntegratedFactorAnalysisLikelihood.CacheProvider.NO_CACHE;
        }


        return new IntegratedFactorAnalysisLikelihood(xo.getId(), traitParameter, missingIndicators,
                loadings, traitPrecision, nugget, taskPool, cacheProvider);
    }

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(LOADINGS, new XMLSyntaxRule[]{
                    new ElementRule(MatrixParameterInterface.class),
            }),
            new ElementRule(PRECISION, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class),
            }),
            // Tree trait parser
            new ElementRule(MutableTreeModel.class),
            AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
            new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }),
            new ElementRule(TreeTraitParserUtilities.MISSING, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
            AttributeRule.newDoubleRule(NUGGET, true),
            AttributeRule.newBooleanRule(STANDARDIZE, true),
            new ElementRule(TaskPool.class, true),
            AttributeRule.newDoubleRule(TARGET_SD, true),
            AttributeRule.newBooleanRule(CACHE_PRECISION, true)

    };

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return IntegratedFactorAnalysisLikelihood.class;
    }

    @Override
    public String getParserName() {
        return INTEGRATED_FACTOR_Model;
    }

    @Override
    public String getParserTag() {
        return getParserName();
    }

    @Override
    public String getId(String modelName) {
        return modelName + ".integratedFactorModel";
    }

    private static final String[] ALLOWABLE_PARAMETERS = new String[]{PRECISION, LOADINGS};

    public BeautiParameterIDProvider getBeautiParameterIDProvider(String parameterKey) {
        for (int i = 0; i < ALLOWABLE_PARAMETERS.length; i++) {
            if (parameterKey.equalsIgnoreCase(ALLOWABLE_PARAMETERS[i])) {
                if (parameterKey.equalsIgnoreCase(PRECISION)) {
                    parameterKey = "residualPrecision";
                }
                return new BeautiParameterIDProvider(parameterKey);
            }
        }
        throw new IllegalArgumentException("Unrecognized parameter key '" + parameterKey + "'");
    }
}
