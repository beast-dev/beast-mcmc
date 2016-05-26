/*
 * CartogramDiffusionModel.java
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

import dr.inference.model.Parameter;
import dr.xml.*;
import dr.geo.cartogram.CartogramMapping;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.IOException;
import java.util.logging.Logger;

/**
 * @author Marc A. Suchard
 */

public class CartogramDiffusionModel extends MultivariateDiffusionModel {

    public static final String DIFFUSION_PROCESS = "cartogramDiffusionModel";
    public static final String FILENAME = "cartogramFileName";
    public static final String DENSITY = "density";

    public static final String BOUNDING_BOX = "boundingBox";
    public static final String MIN_X = "minX";
    public static final String MAX_X = "maxX";
    public static final String MIN_Y = "minY";
    public static final String MAX_Y = "maxY";

    public static final String XSIZE = "xGridSize";
    public static final String YSIZE = "yGridSize";

    public CartogramDiffusionModel(String name, Parameter precision) {
        super();
        this.precision = precision;
        addVariable(precision);
        setId(name);
        Logger.getLogger("dr.evomodel.continuous").info(
            "Constructing cartogram diffusion model '"+ getId() +"': \n" +
            "\tIf you use this model, please reference: Lemey, Drummond and Suchard (in preparation)\n" +
            "\tPrecision: " + precision.getId()
        );
    }

    public void addMapping(CartogramMapping mapping) {
        this.mapping = mapping;
        Logger.getLogger("dr.evomodel.continuous").info(
            "\tMapping  : " + mapping.toString() + "\n"
        );
    }

    protected CartogramMapping getMapping() { return mapping; }

    protected double calculateLogDensity(double[] start, double[] stop, double time) {

        Point2D realStart = new Point2D.Double(start[0],start[1]);
        Point2D realStop  = new Point2D.Double(stop[0], stop[1]);

        final CartogramMapping mapping = getMapping();

        Point2D mappedStart = mapping.map(realStart);
        Point2D mappedStop  = mapping.map(realStop);

        if (mappedStart == null || mappedStop == null) // points outside of cartogram bounding box
            return Double.NEGATIVE_INFINITY;

        final double factor = mapping.getAverageDensity(); // Weighted by average density of different cartograms
        final double distance = mappedStop.distance(mappedStart); // Euclidean distance in mapped space
        final double inverseVariance = precision.getParameterValue(0) / time / Math.pow(factor,0.25);
        // I believe this is a 2D (not 1D) Normal diffusion approx; hence the precision is squared
        // in the normalization constant
        return -LOG2PI + Math.log(inverseVariance) - 0.00 * Math.log(factor) - 0.5*(distance * distance * inverseVariance);
    }

    protected void calculatePrecisionInfo() {
    }


    public static Rectangle2D parseRectangle2D(XMLObject xo) throws XMLParseException {

        XMLObject cxo = (XMLObject) xo.getChild(BOUNDING_BOX);
        double minX = cxo.getAttribute(MIN_X, 0.0);
        double maxX = cxo.getAttribute(MAX_X, 0.0);
        double minY = cxo.getAttribute(MIN_Y, 0.0);
        double maxY = cxo.getAttribute(MAX_Y, 0.0);

        if (maxX - minX <= 0 || maxY - minY <= 0)
            throw new XMLParseException("Bounding box must contain volume");

        return new Rectangle2D.Double(minX, minY, maxX - minX, maxY - minY);
    }

    public static CartogramMapping parseCartogramMapping(XMLObject xo, Rectangle2D boundingBox) throws XMLParseException {
        int xGridSize = xo.getAttribute(XSIZE, 0);
        int yGridSize = xo.getAttribute(YSIZE, 0);

        if (xGridSize <= 1 || yGridSize <= 1)
            throw new XMLParseException("Strictly positive grid sizes required");

        CartogramMapping mapping = new CartogramMapping(xGridSize, yGridSize, boundingBox);

        String fileName = xo.getAttribute(FILENAME, "NONE");

        Logger.getLogger("dr.evomodel.continuous").info(
                "Loading cartogram file: " + fileName + "\n");

        if (xo.hasAttribute(FILENAME)) {

            try {
                mapping.readCartogramOutput(fileName);
            } catch (IOException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        double density = xo.getAttribute(DENSITY,1.0);
        mapping.setAverageDensity(density);
                
        return mapping;
    }

    protected static ElementRule boundingBoxRules = new ElementRule(BOUNDING_BOX, new XMLSyntaxRule[]{
                    AttributeRule.newDoubleRule(MIN_X),
                    AttributeRule.newDoubleRule(MAX_X),
                    AttributeRule.newDoubleRule(MIN_Y),
                    AttributeRule.newDoubleRule(MAX_Y),
                 });

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

         public String getParserName() {
             return DIFFUSION_PROCESS;
         }

         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             Rectangle2D boundingBox = parseRectangle2D(xo);

             CartogramMapping mapping = parseCartogramMapping(xo, boundingBox);

             Parameter diffusionParam = (Parameter) xo.getChild(Parameter.class);

             CartogramDiffusionModel model = new CartogramDiffusionModel(xo.getId(), diffusionParam);
             model.addMapping(mapping);

             return model;
         }

         //************************************************************************
         // AbstractXMLObjectParser implementation
         //************************************************************************

         public String getParserDescription() {
             return "Describes a bivariate diffusion process using cartogram distances.";
         }

         public XMLSyntaxRule[] getSyntaxRules() {
             return rules;
         }

         private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                 new ElementRule(Parameter.class),
                 AttributeRule.newStringRule(FILENAME,true),
                 AttributeRule.newIntegerRule(XSIZE),
                 AttributeRule.newIntegerRule(YSIZE),
                 boundingBoxRules,
         };

         public Class getReturnType() {
             return MultivariateDiffusionModel.class;
         }
     };

    private Parameter precision;
    private CartogramMapping mapping;

}