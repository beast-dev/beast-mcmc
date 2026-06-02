/*
 * AntigenicLikelihoodParser.java
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

package dr.evomodelxml.antigenic;

import dr.evomodel.antigenic.NewAntigenicLikelihood;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.DataTable;
import dr.xml.*;

import java.io.FileReader;
import java.io.IOException;
import java.util.logging.Logger;

import static dr.evomodel.antigenic.NewAntigenicLikelihood.ANTIGENIC_LIKELIHOOD;

public class AntigenicLikelihoodParser extends AbstractXMLObjectParser {

    public static final String FILE_NAME = "fileName";
    public static final String TIP_TRAIT = "tipTrait";
    public static final String LOCATION_SAMPLING = "locationSampling";
    //        public static final String VIRUS_LOCATIONS = "virusLocations";
    public static final String SERUM_LOCATIONS = "serumLocations";
    public static final String MDS_DIMENSION = "mdsDimension";
    public static final String MERGE_SERUM_ISOLATES = "mergeSerumIsolates";
    public static final String DRIFT_INITIAL_LOCATIONS = "driftInitialLocations";
    public static final String INTERVAL_WIDTH = "intervalWidth";
    public static final String MDS_PRECISION = "mdsPrecision";
    public static final String LOCATION_DRIFT = "locationDrift";
    public static final String VIRUS_DRIFT = "virusDrift";
    public static final String SERUM_DRIFT = "serumDrift";
    public static final String VIRUS_AVIDITIES = "virusAvidities";
    public static final String SERUM_POTENCIES = "serumPotencies";
    public static final String SERUM_BREADTHS = "serumBreadths";
    public static final String VIRUS_OFFSETS = "virusOffsets";
    public static final String SERUM_OFFSETS = "serumOffsets";

    private static final String START_DIMENSION = "tipStartDim";

    public String getParserName() {
        return ANTIGENIC_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String fileName = xo.getStringAttribute(FILE_NAME);

        DataTable<String[]> assayTable;
        try {
            assayTable = DataTable.Text.parse(new FileReader(fileName), true, false);
        } catch (IOException e) {
            throw new XMLParseException("Unable to read assay data from file: " + e.getMessage());
        }
        System.out.println("Loaded HI table file: " + fileName);

        boolean mergeSerumIsolates = xo.getAttribute(MERGE_SERUM_ISOLATES, false);

        int mdsDimension = xo.getIntegerAttribute(MDS_DIMENSION);
        double intervalWidth = 0.0;
        if (xo.hasAttribute(INTERVAL_WIDTH)) {
            intervalWidth = xo.getDoubleAttribute(INTERVAL_WIDTH);
        }

        double driftInitialLocations = 0.0;
        if (xo.hasAttribute(DRIFT_INITIAL_LOCATIONS)) {
            driftInitialLocations = xo.getDoubleAttribute(DRIFT_INITIAL_LOCATIONS);
        }

        CompoundParameter tipTraitParameter = null;
        if (xo.hasChildNamed(TIP_TRAIT)) {
            tipTraitParameter = (CompoundParameter) xo.getElementFirstChild(TIP_TRAIT);
        }

        MatrixParameter samplingParameter = null;
        if (xo.hasChildNamed(LOCATION_SAMPLING)) {
            samplingParameter = (MatrixParameter) xo.getElementFirstChild(LOCATION_SAMPLING);
            // TOD Remove
        }

        MatrixParameterInterface serumLocationsParameter = null;
        if (xo.hasChildNamed(SERUM_LOCATIONS)) {
            serumLocationsParameter = (MatrixParameterInterface) xo.getElementFirstChild(SERUM_LOCATIONS);
        }

        Parameter mdsPrecision = (Parameter) xo.getElementFirstChild(MDS_PRECISION);

        Parameter locationDrift = null;
        if (xo.hasChildNamed(LOCATION_DRIFT)) {
            locationDrift = (Parameter) xo.getElementFirstChild(LOCATION_DRIFT);
        }

        Parameter virusDrift = null;
        if (xo.hasChildNamed(VIRUS_DRIFT)) {
            virusDrift = (Parameter) xo.getElementFirstChild(VIRUS_DRIFT);
        }

        Parameter serumDrift = null;
        if (xo.hasChildNamed(SERUM_DRIFT)) {
            serumDrift = (Parameter) xo.getElementFirstChild(SERUM_DRIFT);
        }

        Parameter virusOffsetsParameter = null;
        if (xo.hasChildNamed(VIRUS_OFFSETS)) {
            virusOffsetsParameter = (Parameter) xo.getElementFirstChild(VIRUS_OFFSETS);
        }

        Parameter serumOffsetsParameter = null;
        if (xo.hasChildNamed(SERUM_OFFSETS)) {
            serumOffsetsParameter = (Parameter) xo.getElementFirstChild(SERUM_OFFSETS);
        }

        Parameter serumPotenciesParameter = null;
        if (xo.hasChildNamed(SERUM_POTENCIES)) {
            serumPotenciesParameter = (Parameter) xo.getElementFirstChild(SERUM_POTENCIES);
        }

        Parameter serumBreadthsParameter = null;
        if (xo.hasChildNamed(SERUM_BREADTHS)) {
            serumBreadthsParameter = (Parameter) xo.getElementFirstChild(SERUM_BREADTHS);
        }

        Parameter virusAviditiesParameter = null;
        if (xo.hasChildNamed(VIRUS_AVIDITIES)) {
            virusAviditiesParameter = (Parameter) xo.getElementFirstChild(VIRUS_AVIDITIES);
        }

        int startDim = xo.getAttribute(START_DIMENSION, 1) - 1;

        NewAntigenicLikelihood AGL = new NewAntigenicLikelihood(
                mdsDimension,
                mdsPrecision,
                locationDrift,
                virusDrift,
                serumDrift,
                samplingParameter,
                serumLocationsParameter,
                tipTraitParameter,
                virusOffsetsParameter,
                serumOffsetsParameter,
                serumPotenciesParameter,
                serumBreadthsParameter,
                virusAviditiesParameter,
                assayTable,
                mergeSerumIsolates,
                intervalWidth,
                driftInitialLocations,
                startDim);

        Logger.getLogger("dr.evomodel").info("Using EvolutionaryCartography model. Please cite:\n" + Citable.Utils.getCitationString(AGL));

        return AGL;
    }

//************************************************************************
// AbstractXMLObjectParser implementation
//************************************************************************

    public String getParserDescription() {
        return "Provides the likelihood of immunological assay data such as Hemagglutinin inhibition (HI) given vectors of coordinates" +
                "for viruses and sera/antisera in some multidimensional 'antigenic' space.";
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newStringRule(FILE_NAME, false, "The name of the file containing the assay table"),
            AttributeRule.newIntegerRule(MDS_DIMENSION, false, "The dimension of the space for MDS"),
            AttributeRule.newBooleanRule(MERGE_SERUM_ISOLATES, true, "Should multiple serum isolates from the same strain have their locations merged (defaults to false)"),
            AttributeRule.newDoubleRule(INTERVAL_WIDTH, true, "The width of the titre interval in log 2 space"),
            AttributeRule.newDoubleRule(DRIFT_INITIAL_LOCATIONS, true, "The degree to drift initial virus and serum locations, defaults to 0.0"),
            new ElementRule(TIP_TRAIT, CompoundParameter.class, "Optional parameter of tip locations from the tree", true),
//                new ElementRule(VIRUS_LOCATIONS, MatrixParameter.class, "Parameter of locations of all virus"),
            new ElementRule(SERUM_LOCATIONS, MatrixParameterInterface.class, "Parameter of locations of all sera"),
            new ElementRule(VIRUS_OFFSETS, Parameter.class, "Optional parameter for virus dates to be stored", true),
            new ElementRule(SERUM_OFFSETS, Parameter.class, "Optional parameter for serum dates to be stored", true),
            new ElementRule(SERUM_POTENCIES, Parameter.class, "Optional parameter for serum potencies", true),
            new ElementRule(SERUM_BREADTHS, Parameter.class, "Optional parameter for serum breadths", true),
            new ElementRule(VIRUS_AVIDITIES, Parameter.class, "Optional parameter for virus avidities", true),
            new ElementRule(MDS_PRECISION, Parameter.class, "Parameter for precision of MDS embedding"),
            new ElementRule(LOCATION_DRIFT, Parameter.class, "Optional parameter for drifting locations with time", true),
            new ElementRule(VIRUS_DRIFT, Parameter.class, "Optional parameter for drifting only virus locations, overrides locationDrift", true),
            new ElementRule(SERUM_DRIFT, Parameter.class, "Optional parameter for drifting only serum locations, overrides locationDrift", true),
            AttributeRule.newIntegerRule(START_DIMENSION, true),
    };

    public Class getReturnType() {
        return NewAntigenicLikelihood.class;
    }
}

