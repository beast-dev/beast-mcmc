/*
 * MaximizationOptions.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.treelime.options;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * The Maximization tab's settings: what is maximized, how long the L-BFGS run may take, and which parameters are
 * optimized. Parameters are referred to by id, so the settings survive model changes.
 */
public class MaximizationOptions implements Serializable {

    private static final long serialVersionUID = 1L;

    public enum Objective {
        MAP("Maximum a posteriori (likelihood and priors)", "joint"),
        ML("Maximum likelihood (tree likelihoods only)", "likelihood"),
        PRIOR("Prior only (for checking)", "prior");

        private final String label;
        private final String densityId;

        Objective(String label, String densityId) {
            this.label = label;
            this.densityId = densityId;
        }

        /**
         * @return the id of the density the maximizer maximizes
         */
        public String getDensityId() {
            return densityId;
        }

        public String toString() {
            return label;
        }
    }

    public enum StepSize {
        TINY("tiny"), SMALL("small"), MEDIUM("medium"), LARGE("large");

        private final String name;

        StepSize(String name) {
            this.name = name;
        }

        public String getName() {
            return name;
        }

        public String toString() {
            return name;
        }
    }

    public enum GradientMode {
        ANALYTIC_WHERE_AVAILABLE("Analytic where available, numeric otherwise"),
        NUMERIC("Numeric for everything");

        private final String label;

        GradientMode(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }
    }

    public Objective objective = Objective.MAP;
    private GradientMode gradientMode = GradientMode.ANALYTIC_WHERE_AVAILABLE;
    public int nIterations = 1000; // at least 1: 0 means no limit to lbfgs4j
    public String maxRunTime = ""; // a duration with an optional unit, "" for none
    public long maxEvaluations = 0; // 0 for none
    public boolean printToScreen = true;
    public StepSize stepSize = StepSize.MEDIUM;
    public final Set<String> excluded = new HashSet<String>(); // ids of parameters held at their initial values
    public String reportFileName = null; // null: <stem>.max.txt
    // ids of parameters with a lower bound other than the default; a null value for none. Null in older files
    private Map<String, Double> lowerBounds = new HashMap<String, Double>();

    /**
     * The gamma shape's default: its optimum often lies at 0, where L-BFGS chases a plateau (the discrete gamma tends
     * to a model of invariant sites); the maximizer holds it here instead.
     */
    public static final double GAMMA_SHAPE_LOWER_BOUND = 0.05;

    // null in files saved before the setting existed
    public GradientMode getGradientMode() {
        return gradientMode != null ? gradientMode : GradientMode.ANALYTIC_WHERE_AVAILABLE;
    }

    public void setGradientMode(GradientMode gradientMode) {
        this.gradientMode = gradientMode;
    }

    public boolean isAnalytic() {
        return getGradientMode() == GradientMode.ANALYTIC_WHERE_AVAILABLE;
    }

    /**
     * @return the bound at which the maximizer holds the parameter once a step goes below it (holdAtLowerBound), or
     * null for none
     */
    public Double getLowerBound(Parameter parameter) {
        if (lowerBounds != null && lowerBounds.containsKey(parameter.getName())) {
            return lowerBounds.get(parameter.getName());
        }
        return getDefaultLowerBound(parameter);
    }

    public static Double getDefaultLowerBound(Parameter parameter) {
        return parameter.getBaseName().equals("alpha") ? GAMMA_SHAPE_LOWER_BOUND : null;
    }

    /**
     * @param bound null for none
     */
    public void setLowerBound(Parameter parameter, Double bound) {
        if (lowerBounds == null) {
            lowerBounds = new HashMap<String, Double>();
        }
        lowerBounds.put(parameter.getName(), bound);
    }

    public String getReportFileName(String fileNameStem) {
        return reportFileName != null ? reportFileName : fileNameStem + ".max.txt";
    }

    /**
     * How a parameter enters the maximization: with a transform to the real line, or held fixed and why.
     */
    public enum Treatment {
        LOG("log", null),
        LOGIT("logit", null),
        NONE("none", null),
        NODE_HEIGHTS("node-height ratios: the root and every internal node", null),
        FIXED_BY_PRIOR(null, "fixed in the Priors tab"),
        EXCLUDED(null, "not selected"),
        STATISTIC(null, "a statistic, not a parameter"),
        DISCRETE(null, "discrete"),
        NODE_HEIGHT(null, "node heights need the analytic gradient mode"),
        SIMPLEX(null, "sums to a constant: needs a simplex transform, not available yet"),
        LINKED(null, "linked parameters are not optimized yet"),
        UNSUPPORTED(null, "no suitable transform");

        private final String transform;
        private final String reason;

        Treatment(String transform, String reason) {
            this.transform = transform;
            this.reason = reason;
        }

        public boolean isOptimized() {
            return transform != null;
        }

        /**
         * @return the BEAST transform type, or null when the parameter is held fixed
         */
        public String getTransform() {
            return transform;
        }

        /**
         * @return why the parameter is held fixed, or null when it is optimized
         */
        public String getReason() {
            return reason;
        }

        public String toString() {
            return transform != null ? transform : reason;
        }
    }

    /**
     * @param ignoreSelection true to classify as if every parameter were selected (for the Optimize column)
     */
    public Treatment classify(Parameter parameter, boolean ignoreSelection) {
        if (parameter.isStatistic) {
            return Treatment.STATISTIC;
        }
        if (parameter.isDiscrete) {
            return Treatment.DISCRETE;
        }
        if (parameter.isFixed()) {
            return Treatment.FIXED_BY_PRIOR;
        }
        if (parameter.isNodeHeight) {
            // the root height stands for its tree's whole node-height block
            if (!isAnalytic()) {
                return Treatment.NODE_HEIGHT;
            }
            if (!ignoreSelection && excluded.contains(parameter.getName())) {
                return Treatment.EXCLUDED;
            }
            return Treatment.NODE_HEIGHTS;
        }
        if (parameter.isMaintainedSum) {
            return Treatment.SIMPLEX;
        }
        if (parameter.isLinked) {
            return Treatment.LINKED;
        }
        if (!ignoreSelection && excluded.contains(parameter.getName())) {
            return Treatment.EXCLUDED;
        }
        if (parameter.isZeroOne) {
            return Treatment.LOGIT;
        }
        if (parameter.isNonNegative) {
            return Treatment.LOG;
        }
        if (parameter.isInRealSpace()) {
            return Treatment.NONE;
        }
        return Treatment.UNSUPPORTED;
    }

    /**
     * @return the parameters to optimize, grouped by transform in the order log, logit, none (the order the generator
     * writes them)
     */
    public List<Parameter> selectOptimized(List<Parameter> parameters) {
        List<Parameter> selected = new ArrayList<Parameter>();
        for (Treatment treatment : new Treatment[]{Treatment.LOG, Treatment.LOGIT, Treatment.NONE}) {
            for (Parameter parameter : parameters) {
                if (classify(parameter, false) == treatment) {
                    selected.add(parameter);
                }
            }
        }
        return selected;
    }
}
