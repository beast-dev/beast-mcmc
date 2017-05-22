/*
 * MLOptimizer.java
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

package dr.inference.ml;

import dr.inference.loggers.Logger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.markovchain.MarkovChainListener;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.operators.OperatorSchedule;
import dr.util.Identifiable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A heuristic optimizer that uses the MCMC framework.
 *
 * @author Andrew Rambaut
 *
 * @version $Id: MLOptimizer.java,v 1.5 2006/06/13 03:50:54 alexei Exp $
 */
public class MLOptimizer implements Runnable, Identifiable {

	/** the likelihood function */
	private final Likelihood likelihood;

	/** the operator schedule */
	private OperatorSchedule schedule;

	private final Logger[] loggers;

	int chainLength;

	private String id = null;

	/**
	 * Constructor
	 * @param chainLength the chain length
	 * @param schedule operator schedule to be used in chain.
	 */
	public MLOptimizer(String id,
		int chainLength,
		Likelihood likelihood,
		OperatorSchedule schedule,
		Logger[] loggers) {

		this.id = id;

		mc = new MarkovChain(likelihood, schedule, new GreatDelugeCriterion(0.2), 2000, 1, MarkovChain.EVALUATION_TEST_THRESHOLD, false);
        //mc = new MarkovChain(null, likelihood, schedule, new HillClimbingCriterion(), false);

		this.chainLength = chainLength;
		this.likelihood = likelihood;
		this.loggers = loggers;

		setOperatorSchedule(schedule);

		//initialize transients
		currentState = 0;
	}

	public void run() {
		chain();
	}

	/**
	 * This method actually intiates the MCMC analysis.
	 * the site patterns have been dropped from the site model
	 * to reduce the footprint of the sample.
	 */
	public void chain() {
        currentState = 0;

		if (loggers != null) {
            for (Logger logger : loggers) {
                logger.startLogging();
            }
        }

		timer.start();

		mc.reset();
		timer.start();

			mc.addMarkovChainListener(chainListener);
			mc.runChain(getChainLength(), true/*, 0*/);
			mc.removeMarkovChainListener(chainListener);

		timer.stop();
	}

	/** @return the likelihood function. */
	public Likelihood getLikelihood() { return likelihood; }

	/** @return the timer. */
	public dr.util.Timer getTimer() { return timer; }

	/** set the operator schedule used in this MCMC analysis.
	 * MUST be called before chain!
	 */
	public void setOperatorSchedule(OperatorSchedule sched) {
		this.schedule = sched;
	}

	/** @return the operator schedule used in this MCMC analysis. */
	public OperatorSchedule getOperatorSchedule() { return schedule; }

	/** @return the length of this analysis.*/
	public final int getChainLength() { return chainLength; }

	// TRANSIENT PUBLIC METHODS *****************************************

	/** @return the current state of the MCMC analysis. */
	public final long getCurrentState() { return currentState; }

	/** @return the progress (0 to 1) of the MCMC analysis. */
	public final double getProgress() {
		return (double)currentState / chainLength;
	}

    private final MarkovChainListener chainListener = new MarkovChainListener() {
		@Override
		public void bestState(long state, MarkovChain markovChain, Model bestModel) {
			currentState = state;
		}

		@Override
		public void currentState(long state, MarkovChain markovChain, Model currentModel) {

			currentState = state;

			if (loggers != null) {
				for (Logger logger : loggers) {
					logger.log(state);
				}
			}
		}

		@Override
		public void finished(long chainLength, MarkovChain markovChain) {
			currentState = chainLength;

			if (loggers != null) {
				for (Logger logger : loggers) {
					logger.log(currentState);
					logger.stopLogging();
				}
			}

           /* if (false) {
                NumberFormatter formatter = new NumberFormatter(8);

                System.out.println();
                System.out.println("Operator analysis");
                for (int i =0; i < schedule.getOperatorCount(); i++) {

                    MCMCOperator op = schedule.getOperator(i);
                    double acceptanceProb = MCMCOperator.Utils.getAcceptanceProbability(op);
                    System.out.println(formatter.formatToFieldWidth(op.getOperatorName(), 30) + "\t" + formatter.formatDecimal(op.getMeanDeviation(), 2) + "\t" + formatter.formatDecimal(acceptanceProb, 4));
                }
                System.out.println();
            }
*/
		}
		// for receiving
        /** Called when a new new best likelihood state is found. */
        public void bestLklModel(long state, Model bestModel) {
            currentState = state;
        }

	};

	/**
	 * Creates a DOM element that represents this MCMC analysis.
	 */
	public Element createElement(Document d) {
		throw new RuntimeException("Not implemented!");
	}

	public String getId() { return id; }

	public void setId(String id) { this.id = id; }

	// PRIVATE TRANSIENTS

	private long currentState;

    private final dr.util.Timer timer = new dr.util.Timer();

	/** this markov chain does most of the work. */
	private MarkovChain mc = null;
}

