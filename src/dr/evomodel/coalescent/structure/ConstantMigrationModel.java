/*
 * ConstantMigrationModel.java
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

package dr.evomodel.coalescent.structure;

import dr.evolution.colouring.ColourChangeMatrix;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

/**
 * A wrapper for ConstantPopulation.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: ConstantMigrationModel.java,v 1.4 2006/09/08 14:28:07 rambaut Exp $
 */
public class ConstantMigrationModel extends MigrationModel {
    //
    // Public stuff
    //

    public static String CONSTANT_MIGRATION_MODEL = "constantMigrationModel";
    public static String MIGRATION_RATES = "migrationRates";

    /**
     * Construct demographic model with default settings
     */
    public ConstantMigrationModel(int demeCount, Parameter migrationParameter) {

        this(CONSTANT_MIGRATION_MODEL, demeCount, migrationParameter);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstantMigrationModel(String name, int demeCount, Parameter migrationParameter) {

        super(name);

        this.demeCount = demeCount;
        this.migrationParameter = migrationParameter;
        addVariable(migrationParameter);
        migrationParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, migrationParameter.getDimension()));
    }

    // general functions

    public ColourChangeMatrix getMigrationMatrix() {
        if (colourChangeMatrix == null) {
            colourChangeMatrix = new ColourChangeMatrix(migrationParameter.getParameterValues(), demeCount);
        }
        return colourChangeMatrix;
    }

    public double[] getMigrationRates(double time) {
        return migrationParameter.getParameterValues();
    }


    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        colourChangeMatrix = null;
    }

    protected void restoreState() {
        colourChangeMatrix = null;
    }

    /**
     * Parses an element from an DOM document into a ConstantPopulation.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONSTANT_MIGRATION_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int demeCount = 2;

            XMLObject cxo = xo.getChild(MIGRATION_RATES);
            Parameter migrationParameter = (Parameter) cxo.getChild(Parameter.class);


            return new ConstantMigrationModel(demeCount, migrationParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A migration model representing constant migration rates through time.";
        }

        public Class getReturnType() {
            return ConstantMigrationModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MIGRATION_RATES,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
        };
    };

    //
    // protected stuff
    //

    private int demeCount;
    private Parameter migrationParameter;
    private ColourChangeMatrix colourChangeMatrix = null;
}
