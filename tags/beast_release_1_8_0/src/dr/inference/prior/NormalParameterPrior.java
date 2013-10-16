/*
 * NormalParameterPrior.java
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
//import dr.math.distributions.NormalDistribution;
//import org.w3c.dom.Document;
//import org.w3c.dom.Element;

/**
 * This prior provides a normally distributed expectation on a single parameter
 *
 * @author Alexei Drummond
 * @version $Id: NormalParameterPrior.java,v 1.4 2005/05/24 20:26:00 rambaut Exp $
 */
//public class NormalParameterPrior extends AbstractParameterPrior {
//
//    /**
//     * the normal distribution
//     */
//    NormalDistribution normalDistribution;
//
//    public NormalParameterPrior(Parameter parameter, double mean, double stdev) {
//        this(parameter, -1, mean, stdev);
//    }
//
//    public NormalParameterPrior(Parameter parameter, int dimension, double mean, double stdev) {
//        normalDistribution = new NormalDistribution(mean, stdev);
//        setParameter(parameter);
//        setDimension(dimension);
//    }
//
//    public double getLogPriorComponent(double value) {
//        return Math.log(normalDistribution.pdf(value));
//    }
//
//    public Element createElement(Document d) {
//        Element e = d.createElement("normalPrior");
//        e.setAttribute("mean", normalDistribution.mean() + "");
//        e.setAttribute("stdev", normalDistribution.getSD() + "");
//        return e;
//    }
//
//    public final double getMean() {
//        return normalDistribution.mean();
//    }
//
//    public final double getStdev() {
//        return normalDistribution.getSD();
//    }
//
//    public String toString() {
//        return "Normal(" + formatter.format(getMean()).trim() + ", " + formatter.format(getStdev()).trim() + ")";
//    }
//
//    public String toHTML() {
//        return "<font color=\"#FF00FF\">" + toString() + "</font>";
//    }
//}
