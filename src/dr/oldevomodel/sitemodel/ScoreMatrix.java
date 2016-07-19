/*
 * ScoreMatrix.java
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

package dr.oldevomodel.sitemodel;

import dr.evolution.alignment.SitePatterns;

/**
 * @author alexei
 * @version $Id: ScoreMatrix.java,v 1.2 2005/04/08 15:35:36 rambaut Exp $
 */
public class ScoreMatrix {

    SiteModel siteModel;
    double time;
    double[] matrix;

    double[] logp;

    // base frequencies

    double[][] logOddScores;
    double[][] pscores;
    int stateCount;

    public ScoreMatrix(SiteModel siteModel, double time) {
        this.siteModel = siteModel;

        logp = siteModel.getSubstitutionModel().getFrequencyModel().getFrequencies();

        for (int i = 0; i < logp.length; i++) {
            logp[i] = Math.log(logp[i]);
        }

        stateCount = siteModel.getFrequencyModel().getFrequencyCount();
        matrix = new double[stateCount*stateCount];
        pscores = new double[stateCount][stateCount];
        logOddScores = new double[stateCount][stateCount];

        setTime(time);
    }

    public void setTime(double time) {
        this.time = time;
        siteModel.getSubstitutionModel().getTransitionProbabilities(time, matrix);

        //base frequencies

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                pscores[i][j] = Math.log(matrix[i*stateCount + j]);
                logOddScores[i][j] = pscores[i][j] - logp[i] - logp[j];
            }
        }
    }

    /**
     * @param i state in first sequences
     * @param j state in second sequences
     * @return the score of aligning state i with state j
     */
    public final double getScore(int i, int j) {
        return pscores[i][j] - logp[j];
    }

    public final double getScore(SitePatterns patterns) {

        if (patterns.getPatternLength() != 2) throw new IllegalArgumentException();

        double logL = 0.0;
        for (int i = 0; i < patterns.getPatternCount(); i++) {
            double weight = patterns.getPatternWeight(i);
            int[] pattern = patterns.getPattern(i);

            int x = pattern[0];
            int y = pattern[1];

            if (x < stateCount && y < stateCount) {
                logL += getScore(x,y) * weight;
            }

        }
        return logL;
    }
}
