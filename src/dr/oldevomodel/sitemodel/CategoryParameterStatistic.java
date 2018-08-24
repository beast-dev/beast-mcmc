/*
 * CategoryParameterStatistic.java
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

package dr.oldevomodel.sitemodel;

import dr.inference.model.BooleanStatistic;
import dr.inference.model.Statistic;
import dr.xml.*;


/**
 * Performs test of whether a category parameter has a given minimum number of sites in each category
 *
 * @author Roald Forsberg
 */
public class CategoryParameterStatistic extends Statistic.Abstract implements BooleanStatistic {

    private static String MINIMUM_NUMBER = "minimumNumber";

    public CategoryParameterStatistic(String name,
                                      SampleStateAndCategoryModel siteModel,
                                      int minimumNumber) {

        super(name);
        this.minimumNumber = minimumNumber;
        this.siteModel = siteModel;
        this.categoryCount = siteModel.getCategoryCount();
    }

    public int getDimension() {
        return 1;
    }

    /**
     * @return boolean result of test.
     */
    public double getStatisticValue(int dim) {
        return getBoolean(dim) ? 1.0 : 0.0;
    }

    /**
     * @return boolean result of test.
     */
    public boolean getBoolean(int dim) {
        for (int i = 0; i < categoryCount; i++) {
            if (siteModel.getSitesInCategoryCount(i) < minimumNumber)
                return false;
        }

        return true;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return "categoryParameterStatistic";
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            SampleStateAndCategoryModel siteModel = null;
            String name = xo.getAttribute("name", null);

            int minimum = xo.getAttribute(MINIMUM_NUMBER, 0);

            for (int i = 0; i < xo.getChildCount(); i++) {
                if (xo.getChild(i) instanceof SampleStateAndCategoryModel) {
                    siteModel = (SampleStateAndCategoryModel) xo.getChild(i);
                }
            }

            if (siteModel == null)
                throw new XMLParseException(getParserName() + " must contain a SampleStateAndCategoryModel.");
            if (minimum < 1) throw new XMLParseException(getParserName() + " minimum number must be greater than 0.");
            return new CategoryParameterStatistic(name, siteModel, minimum);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns true if the minimum number of sites in a category are present";
        }

        public Class getReturnType() {
            return CategoryParameterStatistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(NAME, "A name for this statistic for the purposes of logging"),
                AttributeRule.newIntegerRule(MINIMUM_NUMBER),
                new ElementRule(SampleStateAndCategoryModel.class)
        };

    };

    private int minimumNumber;
    private int categoryCount;
    private SampleStateAndCategoryModel siteModel;

}
