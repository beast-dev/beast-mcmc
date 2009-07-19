/*
 * Tutorial1.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.tutorial;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.GammaDistributionModel;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.mcmc.MCMC;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.ScaleOperator;
import dr.inference.trace.TraceAnalysis;
import dr.inference.trace.TraceException;
import dr.util.Attribute;

import java.io.IOException;

/**
 * @author Alexei Drummond
 */
public class Tutorial1 {

    public static void main(String[] arg) throws IOException, TraceException {

        Variable.D shape = new Variable.D("shape", 1.0);
        shape.addBounds(new Parameter.DefaultBounds(100, 0, 1));
        Variable.D scale = new Variable.D("scale", 1.0);
        scale.addBounds(new Parameter.DefaultBounds(100, 0, 1));
        GammaDistributionModel gamma = new GammaDistributionModel(shape, scale);

        DistributionLikelihood likelihood = new DistributionLikelihood(gamma);

        Attribute.Default<double[]> d = new Attribute.Default<double[]>(
                "x", new double[]{1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4, 1, 2, 3, 4});

        likelihood.addData(d);

        MCMC mcmc = new MCMC("tutorial1:gamma");

        MCMCOperator shapeMove = new ScaleOperator(shape, 0.75);
        MCMCOperator scaleMove = new ScaleOperator(scale, 0.75);

        MCLogger logger1 = new MCLogger(100);
        logger1.add(scale);
        logger1.add(shape);

        MCLogger logger2 = new MCLogger("tutorial1.log", 100, false);
        logger2.add(scale);
        logger2.add(shape);

        mcmc.init(100000, likelihood, new MCMCOperator[]{shapeMove, scaleMove}, new Logger[]{logger1, logger2});

        mcmc.chain();

        TraceAnalysis.report("tutorial1.log");
    }
}
