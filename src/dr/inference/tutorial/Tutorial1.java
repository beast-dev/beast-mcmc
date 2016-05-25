/*
 * Tutorial1.java
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

package dr.inference.tutorial;

import dr.inference.distribution.DistributionLikelihood;
import dr.inference.distribution.NormalDistributionModel;
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
 * This tutorial in code covers:
 * (1) The creation of random variables with bounds.
 * (2) The creation of a one-dimensional normal distribution model.
 * (3) Creation of normal likelihood and data
 * (4) Initialization of MCMC chain with proposal distribution and two loggers.
 * (5) Post-run trace analysis
 *
 * @author Alexei Drummond
 */
public class Tutorial1 {

    public static void main(String[] arg) throws IOException, TraceException {

        // constructing random variable representing mean of normal distribution
        Variable.D mean = new Variable.D("mean", 1.0);
        // give mean a uniform prior [-1000, 1000]
        mean.addBounds(new Parameter.DefaultBounds(1000, -1000, 1));

        // constructing random variable representing stdev of normal distribution
        Variable.D stdev = new Variable.D("stdev", 1.0);
        // give stdev a uniform prior [0, 1000]
        stdev.addBounds(new Parameter.DefaultBounds(1000, 0, 1));

        // construct normal distribution model
        NormalDistributionModel normal = new NormalDistributionModel(mean, stdev);

        // construct a likelihood for normal distribution
        DistributionLikelihood likelihood = new DistributionLikelihood(normal);

        // construct data
        Attribute.Default<double[]> d = new Attribute.Default<double[]>(
                "x", new double[]{1, 2, 3, 4, 5, 6, 7, 8, 9});

        // add data (representing 9 independent observations) to likelihood
        likelihood.addData(d);

        // construct two "operators" to be used as the proposal distribution
        MCMCOperator meanMove = new ScaleOperator(mean, 0.75);
        MCMCOperator stdevMove = new ScaleOperator(stdev, 0.75);

        // construct a logger to log progress of MCMC run to stdout (screen)
        MCLogger logger1 = new MCLogger(100);
        logger1.add(mean);
        logger1.add(stdev);

        // construct a logger to log to a log file for later analysis
        MCLogger logger2 = new MCLogger("tutorial1.log", 100, false, 0);
        logger2.add(mean);
        logger2.add(stdev);

        // construct MCMC object
        MCMC mcmc = new MCMC("tutorial1:normal");

        // initialize MCMC with chain length, likelihood, operators and loggers
        mcmc.init(100000, likelihood, new MCMCOperator[]{meanMove, stdevMove}, new Logger[]{logger1, logger2});

        // run the mcmc
        mcmc.chain();

        // perform post-analysis
        TraceAnalysis.report("tutorial1.log");
    }
}
