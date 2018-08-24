/*
 * CategoryOperator.java
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
import dr.math.MathUtils;
import dr.oldevomodel.sitemodel.CategorySampleModel;
import dr.inference.model.Parameter;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;


/**
 * An operator on categories of sites
 *
 * @author Roald Forsberg
 */
public class CategoryOperator extends SimpleMCMCOperator {

    public static final String CATEGORY_OPERATOR = "categoryOperator";

    // dimension of categoryParameter should be set beforehand
    public CategoryOperator(CategorySampleModel siteModel, int siteCount,
                            Parameter categoryParameter, double weight) {
        this.categoryParameter = categoryParameter;
        setWeight(weight);
        this.siteModel = siteModel;
        this.categoryCount = siteModel.getCategoryCount();
        this.siteCount = siteCount;
    }


    /**
     * Alter the category of one site
     */
    public final double doOperation() {

        int randomSite = (int) (MathUtils.nextDouble() * siteCount);

        int currentCategory = (int) categoryParameter.getParameterValue(randomSite);

        siteModel.subtractSitesInCategoryCount(currentCategory);

        int[] temp = new int[categoryCount - 1];

        int count = 0;

        for (int i = 0; i < categoryCount; i++) {
            if (i != currentCategory) {
                temp[count] = i;
                count++;
            }
        }

        int newCategory = temp[(int) (MathUtils.nextDouble() * temp.length)];

        categoryParameter.setParameterValue(randomSite, newCategory);
        siteModel.addSitesInCategoryCount(newCategory);

        return 0.0;

    }


    // Interface MCMCOperator
    public final String getOperatorName() {
        return CATEGORY_OPERATOR;
    }


    /**
     * Create the Operator part of this model parameter!
     */
    public Element createOperatorElement(Document d) {
        throw new RuntimeException("Not implemented!");
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CATEGORY_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Parameter catParam = (Parameter) xo.getChild(Parameter.class);
            CategorySampleModel siteModel = (CategorySampleModel) xo.getChild(CategorySampleModel.class);
            Alignment alignment = (Alignment) xo.getChild(Alignment.class);

            double weight = xo.getDoubleAttribute(WEIGHT);

            return new CategoryOperator(siteModel, alignment.getSiteCount(),
                    catParam, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator on categories of sites.";
        }

        public Class getReturnType() {
            return CategoryOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule("weight"),
                new ElementRule(Parameter.class),
                new ElementRule(CategorySampleModel.class),
                new ElementRule(Alignment.class)
        };

    };

    public String toString() {
        return getOperatorName();
    }

    public String getPerformanceSuggestion() {
        return "";
    }

    // Private instance variables
    private Parameter categoryParameter;

    private CategorySampleModel siteModel;

    private int categoryCount;

    private int siteCount;

}

