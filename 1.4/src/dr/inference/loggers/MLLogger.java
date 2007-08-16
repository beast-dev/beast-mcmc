/*
 * MLLogger.java
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

package dr.inference.loggers;

import dr.inference.model.Likelihood;
import dr.util.Identifiable;
import dr.xml.*;

import java.io.*;

/**
 * A logger that stores maximum likelihood states.
 *
 * @version $Id: MLLogger.java,v 1.21 2005/07/27 22:09:21 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class MLLogger extends MCLogger {

	public static final String LOG_ML = "logML";
	public static final String LIKELIHOOD = "ml";

	private Likelihood likelihood;
	private double bestLikelihood;
	private int bestState;
	private String[] bestValues = null;
	private int logEvery = 0;

	/**
	 * Constructor
	 */
	public MLLogger(Likelihood likelihood, LogFormatter formatter, int logEvery) {

		super(formatter, logEvery);

		this.likelihood = likelihood;
	}

	public void startLogging() {
		bestLikelihood = Double.NEGATIVE_INFINITY;
		bestState = 0;
		bestValues = new String[getColumnCount()];

		if (logEvery > 0) {
			String[] labels = new String[getColumnCount() + 1];

			labels[0] = "state";

			for (int i = 0; i < getColumnCount(); i++) {
				labels[i+1] = getColumnLabel(i);
			}

			logLabels(labels);
		}

		super.startLogging();
	}

	public void log(int state) {

		double lik;

		lik = likelihood.getLogLikelihood();

		if (lik > bestLikelihood) {

			for (int i = 0; i < getColumnCount(); i++) {
				bestValues[i] = getColumnFormatted(i);
			}

			bestState = state;
			bestLikelihood = lik;

			if (logEvery == 1) {

				String[] values = new String[getColumnCount() + 1];

				values[0] = Integer.toString(bestState);

				for (int i = 0; i < getColumnCount(); i++) {
					values[i+1] = bestValues[i];
				}

				logValues(values);
			}
		}

		if (logEvery > 1 && (state % logEvery == 0)) {

			String[] values = new String[getColumnCount() + 1];

			values[0] = Integer.toString(bestState);

			for (int i = 0; i < getColumnCount(); i++) {
				values[i+1] = bestValues[i];
			}

			logValues(values);
		}
	}

	public void stopLogging() {
		String[] values = new String[getColumnCount() + 2];

		values[0] = Integer.toString(bestState);
		values[1] = Double.toString(bestLikelihood);

		for (int i = 0; i < getColumnCount(); i++) {
			values[i+2] = bestValues[i];
		}

		if (logEvery > 0) {
			logValues(values);
		} else {
			String[] labels = new String[getColumnCount() + 2];

			labels[0] = "state";
			labels[1] = "ML";

			for (int i = 0; i < getColumnCount(); i++) {
				labels[i+2] = getColumnLabel(i);
			}

			logLabels(labels);
			logValues(values);
		}

		super.stopLogging();
	}

	public static XMLObjectParser ML_LOGGER_PARSER = new MCLogger.LoggerParser() {

		public String getParserName() { return LOG_ML; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			Likelihood likelihood = (Likelihood)xo.getSocketChild(LIKELIHOOD);

			// logEvery of zero only displays at the end
			int logEvery = 0;

			if (xo.hasAttribute(LOG_EVERY)) {
				logEvery = xo.getIntegerAttribute(LOG_EVERY);
			}

			PrintWriter pw = null;

			if (xo.hasAttribute(FILE_NAME)) {

				String fileName = xo.getStringAttribute(FILE_NAME);
				try {
					File file = new File(fileName);
					String name = file.getName();
					String parent = file.getParent();

					if (!file.isAbsolute()) {
						parent = System.getProperty("user.dir");
					}

					pw = new PrintWriter(new FileOutputStream(new File(parent, name)));
				} catch (FileNotFoundException fnfe) {
					throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
				}

			} else {
				pw = new PrintWriter(System.out);
			}

			LogFormatter formatter = new TabDelimitedFormatter(pw);

			MLLogger logger = new MLLogger(likelihood, formatter, logEvery);

			if (xo.hasAttribute(TITLE)) {
				logger.setTitle(xo.getStringAttribute(TITLE));
			}

			for (int i = 0; i < xo.getChildCount(); i++) {
				Object child = xo.getChild(i);

				if (child instanceof Columns) {

					logger.addColumns(((Columns)child).getColumns());

				} else if (child instanceof Loggable) {

					logger.add((Loggable)child);

				} else if (child instanceof Identifiable) {

					logger.addColumn(new LogColumn.Default(((Identifiable)child).getId(), child));

				} else {

					logger.addColumn(new LogColumn.Default(child.getClass().toString(), child));
				}
			}

			return logger;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************
		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			AttributeRule.newIntegerRule(LOG_EVERY, true),
			new ElementRule(LIKELIHOOD,
				new XMLSyntaxRule[] { new ElementRule(Likelihood.class) }),
			new OrRule(
				new ElementRule(Columns.class, 1, Integer.MAX_VALUE),
				new ElementRule(Loggable.class, 1, Integer.MAX_VALUE)
			)
		};

		public String getParserDescription() {
			return "Logs one or more items every time the given likelihood improves";
		}

		public Class getReturnType() { return MLLogger.class; }
	};
}
