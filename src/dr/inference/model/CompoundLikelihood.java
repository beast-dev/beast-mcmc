/*
 * CompoundLikelihood.java
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

package dr.inference.model;

import dr.util.NumberFormatter;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.logging.Logger;

/**
 * A likelihood function which is simply the product of a set of likelihood functions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: CompoundLikelihood.java,v 1.19 2005/05/25 09:14:36 rambaut Exp $
 */
public class CompoundLikelihood implements Likelihood {

	public static final String COMPOUND_LIKELIHOOD = "compoundLikelihood";
	public static final String THREADS = "threads";
	public static final String POSTERIOR = "posterior";
	public static final String PRIOR = "prior";
	public static final String LIKELIHOOD = "likelihood";

	public CompoundLikelihood(int threads) {
		if (threads > 1) {
			pool = Executors.newFixedThreadPool(threads);
//			pool = Executors.newCachedThreadPool();
		} else {
			pool = null;
		}
	}

	public void addLikelihood(Likelihood likelihood) {

		if (!likelihoods.contains(likelihood)) {

			likelihoods.add(likelihood);
			if (likelihood.getModel() != null) {
				compoundModel.addModel(likelihood.getModel());
			}

			likelihoodCallers.add(new LikelihoodCaller(likelihood));
		}
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

    // todo: remove in release
    static int DEBUG = 0;               

	public double getLogLikelihood() {
		double logLikelihood = 0.0;

		if (pool == null) {
			// Single threaded

			for (Likelihood likelihood : likelihoods) {
			   final double l = likelihood.getLogLikelihood();
				// if the likelihood is zero then short cut the rest of the likelihoods
				// This means that expensive likelihoods such as TreeLikelihoods should
				// be put after cheap ones such as BooleanLikelihoods
				if( l == Double.NEGATIVE_INFINITY )
                    return Double.NEGATIVE_INFINITY;
				logLikelihood += l;
			}
		} else {
			try {

				List<Future<Double>> results = pool.invokeAll(likelihoodCallers);

				for (Future<Double> result : results) {
					double logL = result.get();
					logLikelihood += logL;
				}

			} catch (InterruptedException e) {
				e.printStackTrace();
			} catch (ExecutionException e) {
				e.printStackTrace();
			}
		}

        if( DEBUG > 0 ) {
            int t = DEBUG; DEBUG = 0;
            System.err.println(getId() + ": " + getDiagnosis() + " = " + logLikelihood);
            DEBUG = t;
        }
		return logLikelihood;
	}

    public void makeDirty() {
        for( Likelihood likelihood : likelihoods ) {
            likelihood.makeDirty();
        }
    }

    public String getDiagnosis() {
        String message = "";
        boolean first = true;

        final NumberFormatter nf = new NumberFormatter(6);

        for( Likelihood lik : likelihoods ) {

            if( !first ) {
                message += ", ";
            } else {
                first = false;
            }

            message += lik.prettyName() + "=";

            if( lik instanceof CompoundLikelihood ) {
                final String d = ((CompoundLikelihood) lik).getDiagnosis();
                if( d != null && d.length() > 0 ) {
                    message += "{" + d + "}";
                }
            } else {

                final double logLikelihood = lik.getLogLikelihood();
                if( logLikelihood == Double.NEGATIVE_INFINITY ) {
                    message += "-Inf";
                } else if( Double.isNaN(logLikelihood) ) {
                    message += "NaN";
                } else {
                    message += nf.formatDecimal(logLikelihood, 4);
                }
            }
        }

        return message;
    }

	public String toString() {
        return getId();
        // really bad for debugging
		//return Double.toString(getLogLikelihood());
	}

    public String prettyName() {
        return Abstract.getPrettyName(this);
    }

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

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return COMPOUND_LIKELIHOOD;
		}

		public String[] getParserNames() {
			return new String[]{getParserName(), POSTERIOR, PRIOR, LIKELIHOOD};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			int threads = xo.getAttribute(THREADS, 0);

            if (xo.getName().equalsIgnoreCase(LIKELIHOOD) && System.getProperty("thread_count") != null) {
                threads = Integer.parseInt(System.getProperty("thread_count"));
                if (threads < 0 || threads > 1000) {
                    // put an upper limit here - may be unnecessary?
                    threads = 0;
                }
            }

			CompoundLikelihood compoundLikelihood = new CompoundLikelihood(threads);

			for(int i = 0; i < xo.getChildCount(); i++) {
                final Object child = xo.getChild(i);
                if( child instanceof Likelihood ) {
					compoundLikelihood.addLikelihood((Likelihood) child);
				} else {

                    throw new XMLParseException("An element (" + child + ") which is not a likelihood has been added to a "
                            + COMPOUND_LIKELIHOOD + " element");
				}
			}

			if (threads > 1) {
				Logger.getLogger("dr.evomodel").info("Likelihood is using " + threads + " threads.");
			}

			return compoundLikelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A likelihood function which is simply the product of its component likelihood functions.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
                AttributeRule.newIntegerRule(THREADS, true),
				new ElementRule(Likelihood.class, 1, Integer.MAX_VALUE)				
		};

		public Class getReturnType() {
			return CompoundLikelihood.class;
		}
	};

	private final ExecutorService pool;

	private final ArrayList<Likelihood> likelihoods = new ArrayList<Likelihood>();
	private final CompoundModel compoundModel = new CompoundModel("compoundModel");


	private final List<Callable<Double>> likelihoodCallers = new ArrayList<Callable<Double>>();

	class LikelihoodCaller implements Callable<Double> {

		public LikelihoodCaller(Likelihood likelihood) {
			this.likelihood = likelihood;
		}

		public Double call() throws Exception {
			return likelihood.getLogLikelihood();
		}

		private final Likelihood likelihood;
	}
}

