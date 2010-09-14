/*
 * GammaParameterPrior.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

//package dr.inference.prior;
//
//import dr.inference.model.Parameter;
//import dr.math.distributions.GammaDistribution;
//import org.w3c.dom.Document;
//import org.w3c.dom.Element;

/**
 * This prior provides a gamma distributed expectation on a single parameter
 *
 * @author Alexei Drummond
 * @version $Id: GammaParameterPrior.java,v 1.4 2005/05/24 20:26:00 rambaut Exp $
 */
//public class GammaParameterPrior extends AbstractParameterPrior {
//
//    /**
//     * the gamma distribution
//     */
//    GammaDistribution gammaDistribution;
//
//    public GammaParameterPrior(Parameter parameter, double gammaMean, double gammaShape) {
//        this(parameter, -1, gammaMean, gammaShape);
//    }
//
//    public GammaParameterPrior(Parameter parameter, int dimension, double gammaMean, double gammaShape) {
//        gammaDistribution = new GammaDistribution(gammaShape, gammaMean / gammaShape);
//        setParameter(parameter);
//        setDimension(dimension);
//    }
//
//    public double getLogPriorComponent(double value) {
//        return gammaDistribution.logPdf(value);
//       // return Math.log(gammaDistribution.pdf(value));
//    }
//
//    public Element createElement(Document d) {
//        Element e = d.createElement("gammaPrior");
//        e.setAttribute("mean", gammaDistribution.mean() + "");
//        e.setAttribute("shape", gammaDistribution.getShape() + "");
//        return e;
//    }
//
//    public final double getMean() {
//        return gammaDistribution.mean();
//    }
//
//    public final double getShape() {
//        return gammaDistribution.getShape();
//    }
//
//
//    public String toString() {
//        return "Gamma(" + formatter.format(getMean()).trim() + ", " + formatter.format(getShape()).trim() + ")";
//    }
//
//    public String toHTML() {
//        return "<font color=\"#FF00FF\">" + toString() + "</font>";
//    }
//}
