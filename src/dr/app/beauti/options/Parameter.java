/*
 * Parameter.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.PriorScaleType;
import dr.app.beauti.types.PriorType;
import dr.evolution.util.Taxa;
import dr.math.distributions.Distribution;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Parameter implements Serializable {
    private static final long serialVersionUID = -8008521511485034329L;

    public static final double UNIFORM_MAX_BOUND = 1.0E100;

    private String prefix = null;
    private boolean priorEdited;

    private boolean meanInRealSpace = false;

    // Required para
    private String baseName;
    private final String description;

    private int dimensionWeight = 1;

    private final List<Parameter> subParameters = new ArrayList<Parameter>();

    // final Builder para
    public String taxaId; // needs to change TMRCA stat name. Issue 520
    public final boolean isNodeHeight;
    public final boolean isStatistic;
    public final boolean isDiscrete;
    public final boolean isHierarchical;
    public final boolean isCMTCRate;
    public final boolean isNonNegative;
    public final boolean isZeroOne;
    public final boolean isCached;
    public final boolean isAdaptiveMultivariateCompatible;
    public boolean isMaintainedSum;
    public boolean isCalibratedYule = false;
//    public final double lower;
//    public final double upper;

    private PartitionOptions options;

    public final PriorScaleType scaleType;

    public final boolean isPriorFixed;
    public PriorType priorType;

    private Parameter parent;

    private Taxa taxonSet = null;

    public double getInitial() {
        return initial;
    }

    public void setInitial(double initial) {
        this.initial = initial;
    }

    // Editable fields
    private boolean isFixed;
    private double initial;
    public double maintainedSum;
    public double dimension;
    public boolean isTruncated;
    public double truncationUpper;
    public double truncationLower;
    public double mean;
    public double stdev;
    public double shape;
    public double shapeB;
    public double scale;
    public double offset;
    public double precision;
    public double uniformUpper;
    public double uniformLower;

    public boolean isLinked;
    public String linkedName;

    public static class Builder {
        // Required para
        private final String baseName;
        private final String description;

        // Optional para - initialized to default values
        private PriorScaleType scaleType = PriorScaleType.NONE;

        private String taxaId = null;
        private boolean isNodeHeight = false;
        private boolean isDiscrete = false;
        private boolean isHierarchical = false;
        private boolean isCMTCRate = false;
        private boolean isNonNegative = false;
        private boolean isZeroOne = false;
        private boolean isMaintainedSum = false;
        private boolean isStatistic = false;
        private boolean isCached = false;

        private PartitionOptions options = null;

        private Taxa taxonSet = null;

        private PriorType priorType = PriorType.NONE_TREE_PRIOR;
        private boolean isPriorFixed = false;

        private boolean isAdaptiveMultivariateCompatible = false;

        private double maintainedSum = 1.0;
        private double dimension = 1;
        private double initial = Double.NaN;
        //        private double upper = Double.NaN;
//        private double lower = Double.NaN;
        private boolean isTruncated = false;
        public double truncationUpper = Double.POSITIVE_INFINITY;
        public double truncationLower = Double.NEGATIVE_INFINITY;
        public double mean = 0.0;
        public double stdev = 1.0;
        public double shape = 1.0;
        public double shapeB = 3.0;
        public double scale = 1.0;
        public double offset = 0.0;
        public double precision = 1.0;

        // the uniform distribution has explicit bounds (ignores the truncations):
        public double uniformUpper = UNIFORM_MAX_BOUND;
        public double uniformLower = -UNIFORM_MAX_BOUND;

        private boolean isFixed = false;


        public Builder(String name, String description) {
            this.baseName = name;
            this.description = description;
        }

        public Builder duplicate(Parameter source) {
            scaleType = source.scaleType;
            taxaId = source.taxaId;
            isNodeHeight = source.isNodeHeight;
            isDiscrete = source.isDiscrete;
            isHierarchical = source.isHierarchical;
            isCMTCRate = source.isCMTCRate;
            isNonNegative = source.isNonNegative;
            isZeroOne = source.isZeroOne;
            isMaintainedSum = source.isMaintainedSum;
            isStatistic = source.isStatistic;
            isCached = source.isCached;
            options = source.options;
            priorType = source.priorType;
            isPriorFixed = source.isPriorFixed;
            isAdaptiveMultivariateCompatible = source.isAdaptiveMultivariateCompatible;
            initial = source.initial;
            dimension = source.dimension;
            maintainedSum = source.maintainedSum;
            isTruncated = source.isTruncated;
            truncationUpper = source.truncationUpper;
            truncationLower = source.truncationLower;
            mean = source.mean;
            stdev = source.stdev;
            shape = source.shape;
            shapeB = source.shapeB;
            scale = source.scale;
            offset = source.offset;
            precision = source.precision;
            uniformUpper = source.uniformUpper;
            uniformLower = source.uniformLower;
            isFixed = source.isFixed;
            return this;
        }

        public Builder scaleType(PriorScaleType scaleType) {
            this.scaleType = scaleType;
            return this;
        }

        public Builder initial(double initial) {
            this.initial = initial;
            return this;
        }

        public Builder dimension(int dimension) {
            this.dimension = dimension;
            return this;
        }

        public Builder maintainedSum(double maintainedSum) {
            this.maintainedSum = maintainedSum;
            this.isMaintainedSum = true;
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

        public Builder taxonSet(Taxa taxonSet) {
            this.taxonSet = taxonSet;
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

        public Builder isHierarchical(boolean isHierarchical) {
            this.isHierarchical = isHierarchical;
            return this;
        }

        public Builder isNonNegative(boolean isNonNegative) {
            this.isNonNegative = isNonNegative;
            if (isNonNegative) {
                this.uniformLower = 0.0;
            }
            return this;
        }

        public Builder isZeroOne(boolean isZeroOne) {
            this.isZeroOne = isZeroOne;
            if (isZeroOne) {
                this.uniformLower = 0.0;
                this.uniformUpper = 1.0;
            }
            return this;
        }

        public Builder isMaintainedSum(boolean isMaintainedMean) {
            this.isMaintainedSum = isMaintainedMean;
            return this;
        }

        public Builder isAdaptiveMultivariateCompatible(boolean isAdaptiveMultivariateCompatible) {
            this.isAdaptiveMultivariateCompatible = isAdaptiveMultivariateCompatible;
            return this;
        }


        public Builder isCMTCRate(boolean isCMTCRate) {
            this.isCMTCRate = isCMTCRate;
            return this;
        }

        public Builder isFixed(boolean isFixed) {
            this.isFixed = isFixed;
            return this;
        }

        public Builder isPriorFixed(boolean priorFixed) {
            this.isPriorFixed = priorFixed;
            return this;
        }

        public Builder uniformUpper(double upper) {
            this.uniformUpper = upper;
            return this;
        }

        public Builder uniformLower(double lower) {
            this.uniformLower = lower;
            return this;
        }

        public Builder truncationUpper(double truncationUpper) {
            this.isTruncated = true;
            this.truncationUpper = truncationUpper;
            return this;
        }

        public Builder truncationLower(double truncationLower) {
            this.isTruncated = true;
            this.truncationLower = truncationLower;
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

        public Builder precision(double precision) {
            this.precision = precision;
            return this;
        }

        public Builder shape(double shape) {
            this.shape = shape;
            return this;
        }

        public Builder shapeB(double shapeB) {
            this.shapeB = shapeB;
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
        maintainedSum = builder.maintainedSum;
        dimension = builder.dimension;
        taxaId = builder.taxaId;
        isNodeHeight = builder.isNodeHeight;
        isStatistic = builder.isStatistic;
        isCached = builder.isCached;
        options = builder.options;
        priorType = builder.priorType;
        isDiscrete = builder.isDiscrete;
        isFixed = builder.isFixed;
        isHierarchical = builder.isHierarchical;
        isCMTCRate = builder.isCMTCRate;
        isNonNegative = builder.isNonNegative;
        isZeroOne = builder.isZeroOne;
        isMaintainedSum = builder.isMaintainedSum;
        isPriorFixed = builder.isPriorFixed;
        isAdaptiveMultivariateCompatible = builder.isAdaptiveMultivariateCompatible;

        taxonSet = builder.taxonSet;
        
//        upper = builder.upper;
//        lower = builder.lower;
        isTruncated = builder.isTruncated;
        truncationUpper = builder.truncationUpper;
        truncationLower = builder.truncationLower;
        mean = builder.mean;
        stdev = builder.stdev;
        shape = builder.shape;
        shapeB = builder.shapeB;
        scale = builder.scale;
        offset = builder.offset;
        uniformUpper = builder.uniformUpper;
        uniformLower = builder.uniformLower;

        // ExponentialDistribution(1.0 / mean)
        if (priorType == PriorType.EXPONENTIAL_PRIOR && mean == 0) mean = 1;
        if (priorType == PriorType.LOGNORMAL_PRIOR && meanInRealSpace && mean <= 0) mean = 0.01;
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    public String getPrefix() {
        return prefix;
    }

    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    private String getFullName() {
        if (prefix != null) return prefix + getBaseName();
        return getBaseName();
    }

    public String getBaseName() {
        if (taxaId != null) {
            return taxaId; // do not why use taxaId for tmrca, they seems duplicated
        }
        return baseName;
    }

    public String getName() {
        if (taxaId != null) {
            return "tmrca(" + getFullName() + ")";
        } else {
            return getFullName();
        }
    }

    public void setName(String name) {
        this.baseName = name;
    }


    public String getXMLName() { // only for BeautiTemplate
        if (taxaId != null) {
            return "tmrca_" + getFullName();
        } else {
            return getFullName();
        }
    }

    public String getDescription() {
        if (taxaId != null && options != null) {
            return description + " " + taxaId + " on tree " + options.getName();
        }
        return description;
    }

    public double getPriorExpectationMean() {
        double expMean = 1.0;
        Distribution dist = priorType.getDistributionInstance(this);
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

    public void setOptions(PartitionOptions options) { // need to set, which keeps consistent to taxonSetsTreeModel
        this.options = options;
    }

    public Taxa getTaxonSet() {
        return taxonSet;
    }

    public void setTaxonSet(Taxa taxonSet) {
        this.taxonSet = taxonSet;
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
                (priorType == PriorType.NONE_IMPROPER) ||
                        (priorType == PriorType.UNIFORM_PRIOR && (Double.isInfinite(getLowerBound()) || Double.isInfinite(getUpperBound())))) {
            return true;
        }
        return false;
    }

    public double getLowerBound() {
        double lower = Double.NEGATIVE_INFINITY;

        if (priorType == PriorType.UNIFORM_PRIOR || priorType == PriorType.DISCRETE_UNIFORM_PRIOR) {
            lower = uniformLower;
        }

        if (isNonNegative || isZeroOne) {
            if (lower < 0) lower = 0.0;
        }

        if (isTruncated && !Double.isInfinite(truncationLower)) {
            lower = truncationLower;
        }

        return lower;
    }

    public double getUpperBound() {
        double upper = Double.POSITIVE_INFINITY;

        if (isZeroOne) {
            if (upper > 1) {
                upper = 1.0;
            }
        }

        if (priorType == PriorType.UNIFORM_PRIOR || priorType == PriorType.DISCRETE_UNIFORM_PRIOR) {
            upper = uniformUpper;
        }

        if (isTruncated && !Double.isInfinite(truncationUpper)) {
            upper = truncationUpper;
        }

        return upper;
    }

    public boolean isFixed() {
        return priorType == PriorType.NONE_FIXED;
    }

    public void setFixed(boolean isFixed) {
        if (isFixed) {
            priorType = PriorType.NONE_FIXED;
        }
    }

    public boolean isInRealSpace() {
        return meanInRealSpace;
    }

    public void setMeanInRealSpace(boolean meanInRealSpace) {
        this.meanInRealSpace = meanInRealSpace;
    }

    public int[] getParameterDimensionWeights() {
//        if (getOptions() != null && getOptions() instanceof PartitionSubstitutionModel) {
//            return ((PartitionSubstitutionModel)getOptions()).getPartitionCodonWeights();
//        }
        if (getSubParameters().size() > 0) {
            int[] weights = new int[getSubParameters().size()];
            for (int i = 0; i < weights.length; i++) {
                weights[i] = getSubParameters().get(i).getDimensionWeight();
            }
            return weights;
        }
        return new int[] { dimensionWeight };
    }

    public int getDimensionWeight() {
        return dimensionWeight;
    }

    public void setDimensionWeight(int dimensionWeight) {
        this.dimensionWeight = dimensionWeight;
    }

    public void setParent(Parameter parent) {
        this.parent = parent;
    }

    public Parameter getParent() {
        return parent;
    }

    public void addSubParameter(Parameter parameter) {
        parameter.setParent(this);
        subParameters.add(parameter);
    }

    public void clearSubParameters() {
        subParameters.clear();
    }

    public List<Parameter> getSubParameters() {
        return subParameters;
    }

    @Override
    public String toString() {
        return getName();
    }
}
