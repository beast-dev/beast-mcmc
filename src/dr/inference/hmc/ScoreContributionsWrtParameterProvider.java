/*
 * ScoreContributionsWrtParameterProvider.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 */

package dr.inference.hmc;

/**
 * Optional extension for gradients that can expose additive local scores.
 * The empirical Fisher is the sum of outer products of these score rows.
 */
public interface ScoreContributionsWrtParameterProvider extends GradientWrtParameterProvider {

    double[][] getScoreContributionsLogDensity();
}
