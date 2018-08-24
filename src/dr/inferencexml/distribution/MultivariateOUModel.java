/*
 * MultivariateOUModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.distribution;

import dr.inference.distribution.GeneralizedLinearModel;
import dr.oldevomodel.substmodel.PositiveDefiniteSubstitutionModel;
import dr.oldevomodel.substmodel.SubstitutionModel;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Statistic;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
@Deprecated // GLM stuff is now in inference.glm - this is here for backwards compatibility temporarily
public class MultivariateOUModel extends GeneralizedLinearModel implements Statistic {

    private SubstitutionModel Q;
    private MatrixParameter gamma;
    private double[] time;
    private double[] deltaTime;
    private double[] design;
    private double[] W;
    private double[] initialPriorMean;
    private int K;
    private int Ksquared;
    private int numTimeSteps;
    private double[][] GminusWGWt;

    private MultivariateNormalDistribution initialPrior;

    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown;
    private double logLikelihood;
    private double storedLogLikelihood;
    private boolean conditionalPrecisionKnown = false;
    private boolean storedConditionPrecisionKnown;

    private double[] storedWt;
    private double[] Wt;
    private double[] conditionPrecisionVector;
    private double[] storedConditionPrecisionVector;

    private int[] mapTime;
    private List<Double> deltaTimeList; // todo could just use a Map<Double,Integer>

    public MultivariateOUModel(SubstitutionModel substitutionModel, Parameter dependentParam,
                               MatrixParameter gamma, double[] time, double[] design) {
        super(dependentParam);
        this.Q = substitutionModel;
        this.time = time;
        this.design = design;
        this.gamma = gamma;

        K = substitutionModel.getDataType().getStateCount();
        Ksquared = K * K;

        W = new double[Ksquared];

        initialPriorMean = new double[K]; // todo send this mean in constructor

        StringBuffer sb = new StringBuffer("Constructing a multivariate OU model:\n");
        sb.append("\tOutcome dimension = ");
        sb.append(K);
        Logger.getLogger("dr.inference.distribution").info(sb.toString());

        setupTimes();
        addVariable(gamma);
        addModel(substitutionModel);

    }

    private void setupTimes() {
        deltaTimeList = new ArrayList<Double>();
        numTimeSteps = time.length / K - 1;
        deltaTime = new double[numTimeSteps];
        mapTime = new int[numTimeSteps];
        double currentTime = time[0];
        int index = 0;
        for (int i = 0; i < numTimeSteps; i++) {
            index += K;
            deltaTime[i] = time[index] - currentTime;
            currentTime = time[index];
            if (!deltaTimeList.contains(deltaTime[i])) {
                deltaTimeList.add(deltaTime[i]);
            }
            mapTime[i] = deltaTimeList.indexOf(deltaTime[i]);
            ((PositiveDefiniteSubstitutionModel) Q).addPrecalculatedTime(-deltaTime[i]); // todo get rid of negative sign
        }
        Logger.getLogger("dr.inference.distribution").info(
                "\tTime increments: " + new Vector(deltaTime)
        );

        Wt = new double[Ksquared * deltaTimeList.size()];
        storedWt = new double[Ksquared * deltaTimeList.size()];
        conditionPrecisionVector = new double[Ksquared * deltaTimeList.size()];
        storedConditionPrecisionVector = new double[Ksquared * deltaTimeList.size()];

        calculateConditionPrecision();
    }

    private void calculateConditionPrecision() {

        int index = 0;
        double[] tempW = new double[Ksquared];
        double[][] G = gamma.getParameterAsMatrix();

        for (double deltaTime : deltaTimeList) {

            Q.getTransitionProbabilities(-deltaTime, tempW);

            System.arraycopy(tempW, 0, Wt, Ksquared * index, Ksquared);

            double[][] WG = new double[K][K]; // needs to start with zeros
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    for (int k = 0; k < K; k++)
                        WG[i][j] += tempW[i * K + k] * G[k][j];
                }
            }

            double[][] WGWt = new double[K][K]; // needs to start with zeros
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++) {
                    for (int k = 0; k < K; k++)
                        WGWt[i][j] += WG[i][k] * tempW[j * K + k];
                }
            }

            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++)
                    WGWt[i][j] = G[i][j] - WGWt[i][j];
            }

            WGWt = new Matrix(WGWt).inverse().toComponents();

            for (int i = 0; i < K; i++)
                System.arraycopy(WGWt[i], 0, conditionPrecisionVector, Ksquared * index + K * i, K);

            index++;
        }

        conditionalPrecisionKnown = true;
    }

    public final double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public double calculateLogLikelihood(double[] x) {
        return calculateLogLikelihood();
    }

    public double calculateLogLikelihood() {

        double logLikelihood = 0;
        double[] previous = new double[K];
        double[] current = new double[K];

        double[] tmpHolder;
        double[][] G = gamma.getParameterAsMatrix();
        double[] theta = dependentParam.getParameterValues();
        double[] Xbeta = null;
        boolean hasEffects = getNumberOfFixedEffects() > 0;

        if (!conditionalPrecisionKnown)
            calculateConditionPrecision();

        // Prior on initial time-point

        try {
            if (new Matrix(G).determinant() < 0.01)
                return Double.NEGATIVE_INFINITY;

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }


        int index = 0;

        if (!hasEffects) {
            for (int i = 0; i < K; i++)
                previous[i] = theta[index++];
        } else {
            Xbeta = getXBeta();
            for (int i = 0; i < K; i++) {
                previous[i] = theta[index] - Xbeta[index];
                index++;
            }
        }

        initialPrior = new MultivariateNormalDistribution(initialPriorMean, new Matrix(G).inverse().toComponents());
        logLikelihood += initialPrior.logPdf(previous);

        double save = logLikelihood;
        double save2 = 0;
        int oldMapTime = -1;

        double[][] conditionalPrecision = new double[K][K];

        for (int timeStep = 0; timeStep < numTimeSteps; timeStep++) {

            int thisMapTime = mapTime[timeStep];

            if (thisMapTime != oldMapTime) {

                System.arraycopy(Wt, Ksquared * thisMapTime, W, 0, Ksquared);

                for (int i = 0; i < K; i++)
                    System.arraycopy(conditionPrecisionVector, Ksquared * thisMapTime + K * i, conditionalPrecision[i], 0, K);

                oldMapTime = thisMapTime;
            }

            double[] mean = new double[K];
            int u = 0;
            for (int i = 0; i < K; i++) {
                for (int j = 0; j < K; j++)
                    mean[i] += W[u++] * previous[j];
            }

//			// start of removable part;
//			double[][] WG = new double[K][K];
//			for (int i = 0; i < K; i++) {
//				for (int j = 0; j < K; j++) {
//					for (int k = 0; k < K; k++)
//						WG[i][j] += W[i * K + k] * G[k][j];
//				}
//			}
//
//
//			double[][] WGWt = new double[K][K];
//			for (int i = 0; i < K; i++) {
//				for (int j = 0; j < K; j++) {
//					for (int k = 0; k < K; k++)
//						WGWt[i][j] += WG[i][k] * W[j * K + k];
//				}
//			}
//
//
//			for (int i = 0; i < K; i++) {
//				for (int j = 0; j < K; j++)
//					WGWt[i][j] = G[i][j] - WGWt[i][j];
//
//			}
//
//			double[][] oldPrecision = new Matrix(WGWt).inverse().toComponents();
//
//
//			GminusWGWt = WGWt;

            // calculate density of current time step

            // end of removable part;


            if (!hasEffects) {
                for (int i = 0; i < K; i++)
                    current[i] = theta[index++];
            } else {
                for (int i = 0; i < K; i++) {
                    current[i] = theta[index] - Xbeta[index];
                    index++;
                }
            }

            MultivariateNormalDistribution density = new MultivariateNormalDistribution(
                    mean, conditionalPrecision);

            double partialLogLikelihood = density.logPdf(current);

            if (partialLogLikelihood > 10) {

                return Double.NEGATIVE_INFINITY;

            }

            logLikelihood += partialLogLikelihood;

            // move to next point
            tmpHolder = previous;
            previous = current;
            current = tmpHolder;

        }

        if (logLikelihood > 100) {
            System.err.println("got here end");
            System.err.println("save1 = " + save);
            System.err.println("save2 = " + save2);
            System.exit(-1);
        }

        likelihoodKnown = true;

        return logLikelihood;
    }

    protected boolean confirmIndependentParameters() {
        return true;
    }

    public boolean requiresScale() {
        return true;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        conditionalPrecisionKnown = false;
        likelihoodKnown = false;
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {

        if (parameter == gamma) {
            conditionalPrecisionKnown = false;
        }
        likelihoodKnown = false;

    }

    protected void storeState() {
        System.arraycopy(Wt, 0, storedWt, 0, Wt.length);
        System.arraycopy(conditionPrecisionVector, 0, storedConditionPrecisionVector, 0, conditionPrecisionVector.length);
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedConditionPrecisionKnown = conditionalPrecisionKnown;
    }

    protected void restoreState() {
        double[] holder = Wt;
        Wt = storedWt;
        storedWt = holder;

        holder = conditionPrecisionVector;
        conditionPrecisionVector = storedConditionPrecisionVector;
        storedConditionPrecisionVector = holder;

        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        conditionalPrecisionKnown = storedConditionPrecisionKnown;
    }

    protected void acceptState() {
    }

    public String getStatisticName() {
        return getId();
    }

    public String getDimensionName(int dim) {
        return getId() + dim;
    }

    public void setDimensionNames(String[] names) {
        // do nothing
    }

    public int getDimension() {
        return W.length;
    }

    public double getStatisticValue(int dim) {

//		int x = dim / K;
//		int y = dim - x * K;

//		if( GminusWGWt != null )
//			return GminusWGWt[x][y];
        if (W != null)
            return W[dim];
        return 0;
    }

    @Override
    public double getValueSum() {
        double sum = 0.0;
        for (int i = 0; i < getDimension(); i++) {
            sum += getStatisticValue(i);
        }
        return sum;
    }

    public String getAttributeName() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public double[] getAttributeValue() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }
}
