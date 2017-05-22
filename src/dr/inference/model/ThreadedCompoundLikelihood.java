/*
 * ThreadedCompoundLikelihood.java
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

package dr.inference.model;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.ReentrantLock;

import dr.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.util.NumberFormatter;

/**
 * A likelihood function which is simply the product of a set of likelihood functions.
 *
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @version $Id: CompoundLikelihood.java,v 1.19 2005/05/25 09:14:36 rambaut Exp $
 */
public class ThreadedCompoundLikelihood implements Likelihood {

	public static final boolean DEBUG = false;

	public ThreadedCompoundLikelihood() {
	}

	public ThreadedCompoundLikelihood(List<Likelihood> likelihoods) {
		for (Likelihood likelihood : likelihoods) {
			addLikelihood(likelihood);
		}
	}

	public void addLikelihood(Likelihood likelihood) {

		if (!likelihoods.contains(likelihood)) {

			likelihoods.add(likelihood);
			if (likelihood.getModel() != null) {
				compoundModel.addModel(likelihood.getModel());
			}

			likelihoodCallers.add(new LikelihoodCaller(likelihood));

			//System.err.println("LikelihoodCallers size: " + likelihoodCallers.size());
		}
	}

	@Override
	public Set<Likelihood> getLikelihoodSet() {
		Set<Likelihood> set = new HashSet<Likelihood>();
		for (Likelihood l : likelihoods) {
			set.add(l);
			set.addAll(l.getLikelihoodSet());
		}
		return set;
	}

	public int getLikelihoodCount() {
		return likelihoods.size();
	}

	public final Likelihood getLikelihood(int i) {
		return likelihoods.get(i);
	}

	// **************************************************************
	// Likelihood IMPLEMENTATION
	// **************************************************************

	public Model getModel() {
		return compoundModel;
	}

	public double getLogLikelihood() {

		double logLikelihood = 0.0;

		boolean knownLikelihoods = true;
		for (Likelihood likelihood : likelihoods) {
			if (!((BeagleTreeLikelihood)likelihood).isLikelihoodKnown()) {
				knownLikelihoods = false;
				break;
			} else {
				logLikelihood += likelihood.getLogLikelihood();
			}
		}
		
		if (knownLikelihoods) {

			if (DEBUG) {
				//System.err.println("BTLs are known; total logLikelihood = " + logLikelihood);
				//double check if the total loglikelihood will be identical by recalculating
				double backupLikelihood = logLikelihood;
				logLikelihood = 0.0;

				if (threads == null) {
					// first call so setup a thread for each likelihood...
					threads = new LikelihoodThread[likelihoodCallers.size()];
					for (int i = 0; i < threads.length; i++) {
						// and start them running...
						threads[i] = new LikelihoodThread();
						threads[i].start();
					}
				}

				for (int i = 0; i < threads.length; i++) {
					// set the caller which will be called in each thread
					threads[i].setCaller(likelihoodCallers.get(i));
				}

				for (LikelihoodThread thread : threads) {
					// now wait for the results to be set...
					Double result = thread.getResult();
					while (result == null) {
						result = thread.getResult();
					}
					logLikelihood += result;
				}

				if (backupLikelihood != logLikelihood) {
					//System.err.println("Likelihood recalculation does not return stored likelihood");
					throw new RuntimeException("Likelihood recalculation does not return stored likelihood");
				}
			}
		} else {
			//System.err.println("BTLs are not known: recalculate");

			logLikelihood = 0.0;

			//double start = System.nanoTime();
			//System.err.println("TCL getLogLikelihood()");

			if (threads == null) {
				//System.err.println("threads == null");
				// first call so setup a thread for each likelihood...
				threads = new LikelihoodThread[likelihoodCallers.size()];
				//System.err.println("LikelihoodThreads: " + threads.length);
				for (int i = 0; i < threads.length; i++) {
					// and start them running...
					threads[i] = new LikelihoodThread();
					threads[i].start();
				}
			}

			//double setStart = System.nanoTime();
			for (int i = 0; i < threads.length; i++) {
				// set the caller which will be called in each thread
				threads[i].setCaller(likelihoodCallers.get(i));
			}
			//double setEnd = System.nanoTime();
			//System.err.println("setting callers: " + (setEnd - setStart));

			//start = System.nanoTime();
			for (LikelihoodThread thread : threads) {
				//double testone = System.nanoTime();
				// now wait for the results to be set...
				Double result = thread.getResult();
				while (result == null) {
					result = thread.getResult();
				}
				logLikelihood += result;
				//double testtwo = System.nanoTime();
				//System.err.println(thread.getName() + " - result = " + result + ": " + (testtwo - testone));
			}
			//end = System.nanoTime();
			//double end = System.nanoTime();
			//System.err.println("TCL total time: " + (end - start));
		}
		return logLikelihood; // * weightFactor;

	}

	public boolean evaluateEarly() {
		return false;
	}

	public void makeDirty() {

		for (Likelihood likelihood : likelihoods) {
			likelihood.makeDirty();
		}
	}

	public String prettyName() {
		return Abstract.getPrettyName(this);
	}

	public String getDiagnosis() {
		String message = "";
		boolean first = true;

		for (Likelihood lik : likelihoods) {

			if (!first) {
				message += ", ";
			} else {
				first = false;
			}

			String id = lik.getId();
			if (id == null || id.trim().length() == 0) {
				String[] parts = lik.getClass().getName().split("\\.");
				id = parts[parts.length - 1];
			}

			message += id + "=";


			if (lik instanceof ThreadedCompoundLikelihood) {
				String d = ((ThreadedCompoundLikelihood) lik).getDiagnosis();
				if (d != null && d.length() > 0) {
					message += "(" + d + ")";
				}
			} else {

				if (lik.getLogLikelihood() == Double.NEGATIVE_INFINITY) {
					message += "-Inf";
				} else if (Double.isNaN(lik.getLogLikelihood())) {
					message += "NaN";
				} else {
					NumberFormatter nf = new NumberFormatter(6);
					message += nf.formatDecimal(lik.getLogLikelihood(), 4);
				}
			}
		}

		return message;
	}

	public String toString() {

		return Double.toString(getLogLikelihood());

	}

	public void setWeightFactor(double w) { weightFactor = w; }

	public double getWeightFactor() { return weightFactor; }

	// **************************************************************
	// Loggable IMPLEMENTATION
	// **************************************************************

	/**
	 * @return the log columns.
	 */
	public dr.inference.loggers.LogColumn[] getColumns() {
		return new dr.inference.loggers.LogColumn[]{
				new LikelihoodColumn(getId())
		};
	}

	private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
		public LikelihoodColumn(String label) {
			super(label);
		}

		public double getDoubleValue() {
			return getLogLikelihood();
		}
	}

	// **************************************************************
	// Identifiable IMPLEMENTATION
	// **************************************************************

	private String id = null;

	public void setId(String id) {
		this.id = id;
	}

	public String getId() {
		return id;
	}

	private LikelihoodThread[] threads;

	private final ArrayList<Likelihood> likelihoods = new ArrayList<Likelihood>();
	private final CompoundModel compoundModel = new CompoundModel("compoundModel");

	private final List<LikelihoodCaller> likelihoodCallers = new ArrayList<LikelihoodCaller>();

	private double weightFactor = 1.0;

	class LikelihoodCaller {

		public LikelihoodCaller(Likelihood likelihood) {
			this.likelihood = likelihood;
		}

		public double call() {
			return likelihood.getLogLikelihood();
		}

		private final Likelihood likelihood;
	}

	class LikelihoodThread extends Thread {

		public LikelihoodThread() {
		}

		public void setCaller(LikelihoodCaller caller) {
			lock.lock();
			resultAvailable = false;
			try {
				this.caller = caller;
				condition.signal();
			} finally {
				lock.unlock();
			}
		}

		/**
		 * Main run loop
		 */
		 public void run() {
			while (true) {
				lock.lock();
				try {
					while( caller == null)
						condition.await();
					result = caller.call(); // SLOW
					resultAvailable = true;
					caller = null;
				} catch (InterruptedException e){

				} finally {
					lock.unlock();
				}
			}
		 }

		 public Double getResult() {
			 Double returnValue = null;
			 if (!lock.isLocked() && resultAvailable)  { // thread is not busy and completed
				 resultAvailable = false; // TODO need to lock before changing resultAvailable?
				 returnValue = result;
			 }
			 return returnValue;
		 }

		 private LikelihoodCaller caller = null;
		 private Double result = Double.NaN;
		 private boolean resultAvailable = false;
		 private final ReentrantLock lock = new ReentrantLock();
		 private final Condition condition = lock.newCondition();
	}

	public boolean isUsed() {
		return isUsed;
	}

	public void setUsed() {
		isUsed = true;

		for (Likelihood l : likelihoods) {
			l.setUsed();
		}
	}

	private boolean isUsed = false;

}