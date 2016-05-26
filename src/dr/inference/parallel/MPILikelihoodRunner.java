/*
 * MPILikelihoodRunner.java
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

package dr.inference.parallel;

import dr.inference.loggers.Logger;
import dr.inference.model.AbstractModel;
import dr.inference.model.Likelihood;
import dr.util.Identifiable;
import dr.xml.*;
import mpi.MPI;

import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */

public class MPILikelihoodRunner implements Runnable, Identifiable {

	public static final String PARALLEL_CALCULATOR = "parallelCalculator";


	public MPILikelihoodRunner(String id, Likelihood likelihood) {
		this.id = id;
		this.likelihood = likelihood;
	}

	boolean terminate = false;

	public void run() {

		while (!terminate) {
			// Receive msg
			ServiceRequest request = MPIServices.getRequest(0);
			//     System.err.println("Process "+mpiRank+" received request "+request.getId());
			switch (request) {
				case calculateLikeliood:
					calculateLikelihood();
					break;
				case terminateProcess:
					terminateMe();
				default:
					break;
			}
		}

		// finalize();
	}


	private void calculateLikelihood() {
		// Receive current state
//		System.err.println("did i get here?");
		//	MPI.Finalize();
		//	System.exit(-1);
		AbstractModel model = (AbstractModel) likelihood.getModel();
//		System.err.println("trying to receive model state");
		model.receiveState(0);
		likelihood.makeDirty();
		model.fireModelChanged();
//		System.err.println("state received");

		// perform calculation
		double logLikelihood = likelihood.getLogLikelihood();
//		double logLikelihood = likelihood.
//		System.err.println("likelihood = "+logLikelihood+" on process "+mpiRank);
		MPIServices.sendDouble(logLikelihood, 0);
	}

	private void terminateMe() {
		terminate = true;
		//finalize();
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;

	}

	/**
	 * Initialize connection to MPI world
	 */
	public void init() {
		//  MPI.Init(null);
		mpiRank = MPI.COMM_WORLD.Rank();
		mpiSize = MPI.COMM_WORLD.Size();
	}

	public void finalize() {
		MPI.Finalize();
	}

	/**
	 * XML Parser
	 */

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return PARALLEL_CALCULATOR;
		}

		/**
		 * @return a tree object based on the XML element it was passed.
		 */
		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);

			MPILikelihoodRunner runner = new MPILikelihoodRunner("slave", likelihood);
			runner.init();
			//         MCMCOptions options = new MCMCOptions();
			//          OperatorSchedule opsched = (OperatorSchedule)xo.getChild(OperatorSchedule.class);
			ArrayList<Object> loggers = new ArrayList<Object>();

			//         options.setChainLength(xo.getIntegerAttribute(CHAIN_LENGTH));


			for (int i = 0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);
				if (child instanceof Logger) {
					loggers.add(child);
				}
			}

			Logger[] loggerArray = new Logger[loggers.size()];
			loggers.toArray(loggerArray);

			java.util.logging.Logger.getLogger("dr.inference").info("Creating the parallelCalculator chain:");
			/*+
					"\n  chainLength=" + options.getChainLength() +
					"\n  autoOptimize=" + options.useCoercion());*/

			runner.init();

//                MarkovChain mc = mcmc.getMarkovChain();
//                double initialScore = mc.getCurrentScore();

//                if (initialScore == Double.NEGATIVE_INFINITY) {
//                    throw new IllegalArgumentException("The initial model is invalid because it has zero likelihood!");
//                }


			return runner;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns an MPI-based likelihood calculator.";
		}

		public Class getReturnType() {
			return MPILikelihoodRunner.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
//                AttributeRule.newIntegerRule(CHAIN_LENGTH),
//                AttributeRule.newBooleanRule(COERCION, true),
//                AttributeRule.newIntegerRule(PRE_BURNIN, true),
//                new ElementRule(OperatorSchedule.class ),
				new ElementRule(Likelihood.class),
//                new ElementRule(Logger.class, 1, Integer.MAX_VALUE )
		};

	};

	// Private variables

	private String id;
	private int mpiRank;
	private int mpiSize;
	private final Likelihood likelihood;

}
