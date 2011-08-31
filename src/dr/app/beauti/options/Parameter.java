/*
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.beauti.options;

import dr.app.beauti.enumTypes.PriorScaleType;
import dr.app.beauti.enumTypes.PriorType;
import dr.math.distributions.Distribution;

import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Parameter {

    private String prefix = null;
    private boolean priorEdited;

    private boolean meanInRealSpace = false;

    // Required para
    private final String baseName;
    private final String description;

    // final Builder para
    public final String taxaId;
    public final boolean isNodeHeight;
    public final boolean isStatistic;
    public final boolean isCached;

    private final PartitionOptions options;

    // editable Builder para
    public PriorScaleType scaleType;
    public double initial;

    public boolean isFixed;
    public boolean isDiscrete;

    public boolean priorFixed;

    public PriorType priorType;

    public double lower;
    public double upper;
    public double mean;
    public double stdev;
    public double shape;
    public double scale;
    public double offset;
    public double uniformUpper; // http://code.google.com/p/beast-mcmc/issues/detail?id=491
    public double uniformLower;

    public static class Builder {
        // Required para
        private final String baseName;
        private final String description;

        // Optional para - initialized to default values
        private PriorScaleType scaleType = PriorScaleType.NONE;
        private double initial = Double.NaN;

        private String taxaId = null;
        private boolean isNodeHeight = false;
        private boolean isStatistic = false;
        private boolean isCached = false;
        private PartitionOptions options = null;

        private PriorType priorType = PriorType.NONE_TREE_PRIOR;
        private double upper = Double.NaN;
        private double lower = Double.NaN;
        public double mean = 0.0;
        public double stdev = 1.0;
        public double shape = 1.0;
        public double scale = 1.0;
        public double offset = 0.0;
        public double uniformUpper = 0.0;
        public double uniformLower = 0.0;

        private boolean isDiscrete = false;
        private boolean isFixed = false;

        private boolean priorFixed = false;


        public Builder(String name, String description) {
            this.baseName = name;
            this.description = description;
        }

        public Builder scaleType(PriorScaleType scaleType) {
            this.scaleType = scaleType;
            return this;
        }

        public Builder initial(double initial) {
            this.initial = initial;
            return this;
        }

        public Builder taxaId(String taxaId) {
            this.taxaId = taxaId;
            return this;
        }

        public Builder isNodeHeight(boolean isNodeHeight) {
            this.isNodeHeight = isNodeHeight;
            return this;
        }

        public Builder isCached(boolean isCached) {
            this.isCached = isCached;
            return this;
        }

        public Builder isStatistic(boolean isStatistic) {
            this.isStatistic = isStatistic;
            return this;
        }

        public Builder partitionOptions(PartitionOptions options) {
            this.options = options;
            return this;
        }

        public Builder prior(PriorType priorType) {
            this.priorType = priorType;
            return this;
        }

        public Builder isDiscrete(boolean isDiscrete) {
            this.isDiscrete = isDiscrete;
            return this;
        }

        public Builder isFixed(boolean isFixed) {
            this.isFixed = isFixed;
            return this;
        }

        public Builder priorFixed(boolean priorFixed) {
            this.priorFixed = priorFixed;
            return this;
        }

        public Builder upper(double upper) {
            this.upper = upper;
            return this;
        }

        public Builder lower(double lower) {
            this.lower = lower;
            return this;
        }

        public Builder uniformUpper(double uniformUpper) {
            this.uniformUpper = uniformUpper;
            return this;
        }

        public Builder uniformLower(double uniformLower) {
            this.uniformLower = uniformLower;
            return this;
        }

        public Builder mean(double mean) {
            this.mean = mean;
            return this;
        }

        public Builder stdev(double stdev) {
            this.stdev = stdev;
            return this;
        }

        public Builder shape(double shape) {
            this.shape = shape;
            return this;
        }

        public Builder scale(double scale) {
            this.scale = scale;
            return this;
        }

        public Builder offset(double offset) {
            this.offset = offset;
            return this;
        }

        public Parameter build() {
            return new Parameter(this);
        }

        public Parameter build(Map<String, Parameter> map) {
            final Parameter parameter = new Parameter(this);
            map.put(baseName, parameter);
            return parameter;
        }
    }

    private Parameter(Builder builder) {
        baseName = builder.baseName;
        description = builder.description;
        scaleType = builder.scaleType;
        initial = builder.initial;
        taxaId = builder.taxaId;
        isNodeHeight = builder.isNodeHeight;
        isStatistic = builder.isStatistic;
        isCached = builder.isCached;
        options = builder.options;
        priorType = builder.priorType;
        isDiscrete = builder.isDiscrete;
        isFixed = builder.isFixed;
        priorFixed = builder.priorFixed;
        upper = builder.upper;
        lower = builder.lower;
        mean = builder.mean;
        stdev = builder.stdev;
        shape = builder.shape;
        scale = builder.scale;
        offset = builder.offset;
        uniformUpper = builder.uniformUpper;
        uniformLower = builder.uniformLower;

        // ExponentialDistribution(1.0 / mean)
        if (priorType == PriorType.EXPONENTIAL_PRIOR && mean == 0) mean = 1;
        if (priorType == PriorType.LOGNORMAL_PRIOR && meanInRealSpace && mean <= 0) mean = 0.01;
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    private String getFullName() {
        if (prefix != null) return prefix + baseName;
        return baseName;
    }

    public String getBaseName() {
        return baseName;
    }

    public String getName() {
        if (taxaId != null) {
            return "tmrca(" + taxaId + ")";
        } else {
            return getFullName();
        }
    }

    public String getXMLName() { // only for BeautiTemplate
        if (taxaId != null) {
            return "tmrca_" + taxaId;
        } else {
            return getFullName();
        }
    }

    public String getDescription() {
        if (taxaId != null) {
            return "tmrca statistic for taxon set " + taxaId;
//        } else if (prefix != null) {
//            return description + " of partition " + prefix;
        }
        return description;
    }

    public double getPriorExpectationMean() {
        double expMean = 1.0;
        Distribution dist = priorType.getDistributionClass(this);
        if (dist != null) {
            expMean = dist.mean();

            if (expMean == 0) {
                expMean = dist.quantile(0.975);
            }

            if (expMean == 0) {
                expMean = 1.0;
            }
        }

        return expMean;
    }

    public PartitionOptions getOptions() {
        return options;
    }

    public void setPriorEdited(boolean priorEdited) {
        this.priorEdited = priorEdited;
    }

    public boolean isPriorEdited() {
        return priorEdited;
    }

    public boolean isPriorImproper() {
        if (
            // 1/x is an improper prior but we probably don't want to flag it as
            // such (or we want to make a more explicit distinction about when it
            // might be appropriate:
            /* priorType == PriorType.ONE_OVER_X_PRIOR || */
                (priorType == PriorType.UNIFORM_PRIOR && (Double.isInfinite(uniformUpper) || Double.isInfinite(uniformLower)))) {
            return true;
        }
        return false;
    }

    public boolean isMeanInRealSpace() {
        return meanInRealSpace;
    }

    public void setMeanInRealSpace(boolean meanInRealSpace) {
        this.meanInRealSpace = meanInRealSpace;
    }

}
