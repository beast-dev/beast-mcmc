/*
 * RootTraitStatistic.java
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

package dr.evomodel.continuous;

import dr.inference.model.Statistic;
import dr.xml.*;

/**
 *
 * Gabriela Cybis
 *
 *  */
public class RootTraitStatistic extends Statistic.Abstract{


    private int dimension;
    public static final String SAMPLED_ROOT_TRAITS = "sampledRootTraits";
    private IntegratedMultivariateTraitLikelihood likelihood;

    public RootTraitStatistic( IntegratedMultivariateTraitLikelihood likelihood, String id){
        this.likelihood=likelihood;
        setId(id);

        this.dimension = likelihood.dimTrait;



    }




    public int getDimension() {
        return dimension;
    }
    public String getDimensionName(int dim){
        return "root." + likelihood.getTraitName() + dim;
    }

    public double getStatisticValue(int dim) {

    return    likelihood.getRootNodeTrait()[dim];
    }








public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        IntegratedMultivariateTraitLikelihood likelihood = (IntegratedMultivariateTraitLikelihood)
                xo.getChild(IntegratedMultivariateTraitLikelihood.class);

        String id = xo.getId();

        return new RootTraitStatistic(likelihood, id);
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[] {
                new ElementRule(IntegratedMultivariateTraitLikelihood.class),

        };
    }

    public String getParserDescription() {
        return null;
    }

    public Class getReturnType() {
        return RootTraitStatistic.class;
    }

    public String getParserName() {
        return SAMPLED_ROOT_TRAITS;
    }
};

}