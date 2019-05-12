/*
 * PiecewisePopulationModelParser.java
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

package dr.evomodelxml.coalescent;

import com.sun.javafx.tools.packager.Param;
import dr.evolution.util.Units;
import dr.evomodel.coalescent.PiecewisePopulationSizeModel;
import dr.evomodel.coalescent.PopulationSizeModel;
import dr.evomodel.coalescent.demographicmodels.PiecewisePopulationModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Andrew Rambaut
 */
public class PiecewisePopulationSizeModelParser extends AbstractXMLObjectParser {

    public static final String PIECEWISE_POPULATION_SIZE = "piecewisePopulationSize";
    public static final String EPOCHS = "epochs";
    public static final String POPULATION_SIZE = "populationSize";
    public static final String EPOCH_DURATIONS = "epochDurations";

    public static final String WIDTHS = "widths";
    public static final String LINEAR = "linear";

    public String getParserName() {
        return PIECEWISE_POPULATION_SIZE;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        Parameter epochDurations = (Parameter)xo.getElementFirstChild(EPOCH_DURATIONS);
        Parameter populationSize = (Parameter)xo.getElementFirstChild(POPULATION_SIZE);

        XMLObject cxo = xo.getChild(EPOCHS);
        List<PopulationSizeModel> populationSizeModels = cxo.getAllChildren(PopulationSizeModel.class);

        return new PiecewisePopulationSizeModel(PIECEWISE_POPULATION_SIZE, populationSize,
                populationSizeModels,
                epochDurations, units);
    }


    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "This element represents a piecewise population size model";
    }

    public Class getReturnType() {
        return PiecewisePopulationSizeModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(EPOCHS,
                    new XMLSyntaxRule[]{new ElementRule(PopulationSizeModel.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(POPULATION_SIZE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(EPOCH_DURATIONS,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}
