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
import dr.evolution.util.TaxonList;
import dr.math.distributions.ExponentialDistribution;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Parameter {
    
    private String prefix = null;
    private boolean priorEdited;
        
    // Required para
    private final String baseName;
    private final String description;
    
    // final Builder para
    public final TaxonList taxa;
    public final boolean isNodeHeight;
    public final boolean isStatistic;
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
    public double stdev = 1.0;
    public double mean = 1.0;
    public double offset = 0.0;    
    public double shape = 1.0;
    public double scale = 1.0;    
    
    public static class Builder { 
        // Required para
        private final String baseName;
        private final String description;
        
        // Optional para - initialized to default values
        private PriorScaleType scaleType = PriorScaleType.NONE;
        private double initial = Double.NaN;
        
        private TaxonList taxa = null;
        private boolean isNodeHeight = false;
        private boolean isStatistic = false;
        private PartitionOptions options = null;    
        
        private PriorType priorType = PriorType.NONE;
        private double upper = Double.NaN;
        private double lower = Double.NaN;
        
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
        
        public Builder taxa(TaxonList taxa) {
            this.taxa = taxa;
            return this;
        }
        
        public Builder isNodeHeight(boolean isNodeHeight) {
            this.isNodeHeight = isNodeHeight;
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
        
        public Builder upper(double upper) {
            this.upper = upper;
            return this;
        }
        
        public Builder lower(double lower) {
            this.lower = lower;
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
        
       
        public Parameter build() {
            return new Parameter(this);
        }               
    }
    
    private Parameter(Builder builder) {
        baseName = builder.baseName;
        description = builder.description;
        scaleType = builder.scaleType;
        initial = builder.initial;        
        taxa = builder.taxa;
        isNodeHeight = builder.isNodeHeight;
        isStatistic = builder.isStatistic;
        options = builder.options;            
        priorType = builder.priorType;
        upper = builder.upper;
        lower = builder.lower;        
        isDiscrete = builder.isDiscrete;
        isFixed = builder.isFixed;        
        priorFixed = builder.priorFixed;
    } 
    
    
    /**
     * A constructor for "special" parameters which are not user-configurable
     *
     * @param name        the name
     * @param description the description
     */
    public Parameter(String name, String description) {
        this.baseName = name;
        this.description = description;
        this.scaleType = PriorScaleType.NONE;
        this.isNodeHeight = false;
        this.isStatistic = false;
        this.taxa = null;
        this.options = null;
        this.priorType = PriorType.NONE;
        this.initial = Double.NaN;
        this.lower = Double.NaN;
        this.upper = Double.NaN;
    }

    public Parameter(String name, String description, PriorScaleType scale,
                     double initial, double lower, double upper) {
        this.baseName = name;
        this.description = description;
        this.initial = initial;
        this.isNodeHeight = false;
        this.isStatistic = false;

        this.taxa = null;
        this.options = null;
        
        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scaleType = scale;
        this.setPriorEdited(false);
        this.lower = lower;
        this.upper = upper;
//
//        uniformLower = lower;
//        uniformUpper = upper;
    }

    public Parameter(TaxonList taxa, String description) {
        this.taxa = taxa;
        this.baseName = taxa.getId();
        this.description = description;

        this.options = null;
        
        this.isNodeHeight = true;
        this.isStatistic = true;
        this.priorType = PriorType.NONE;
        this.scaleType = PriorScaleType.TIME_SCALE;
        this.setPriorEdited(false);
        this.lower = 0.0;
        this.upper = Double.MAX_VALUE;

//        uniformLower = lower;
//        uniformUpper = upper;
    }

    public Parameter(String name, String description, boolean isDiscrete) {
        this.taxa = null;
        this.options = null;

        this.baseName = name;
        this.description = description;

        this.isNodeHeight = false;
        this.isStatistic = true;
        this.isDiscrete = isDiscrete;
        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scaleType = PriorScaleType.NONE;
        this.setPriorEdited(false);
        this.initial = Double.NaN;
        this.lower = Double.NaN;
        this.upper = Double.NaN;
    }

    public Parameter(String name, String description, double lower, double upper) {
        this.taxa = null;
        this.options = null;

        this.baseName = name;
        this.description = description;

        this.isNodeHeight = false;
        this.isStatistic = true;
        this.isDiscrete = false;
        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scaleType = PriorScaleType.NONE;
        this.setPriorEdited(false);
        this.initial = Double.NaN;
        this.lower = lower;
        this.upper = upper;

//        uniformLower = lower;
//        uniformUpper = upper;
    }
    
    public Parameter(PartitionOptions options, String name, String description, PriorScaleType scale, double initial,
            double lower, double upper) {
        this.baseName = name;
        this.description = description;
        this.initial = initial;
        this.isNodeHeight = false;
        this.isStatistic = false;

        this.taxa = null;
        this.options = options;

        this.priorType = PriorType.UNIFORM_PRIOR;
        this.scaleType = scale;
        this.setPriorEdited(false);
        this.lower = lower;
        this.upper = upper;

//        uniformLower = lower;
//        uniformUpper = upper;
    }

    public Parameter(PartitionOptions options, String name, String description, boolean isNodeHeight,
                     double initial, double lower, double upper) {
        this.baseName = name;
        this.description = description;
        this.initial = initial;

        this.taxa = null;        
        this.options = options;

        this.isNodeHeight = isNodeHeight;
        this.isStatistic = false;
        this.priorType = PriorType.NONE;
        this.scaleType = PriorScaleType.TIME_SCALE;
        this.setPriorEdited(false);
        this.lower = lower;
        this.upper = upper;

//        uniformLower = lower;
//        uniformUpper = upper;
    }

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
        if (taxa != null) {
            return "tmrca(" + taxa.getId() + ")";
        } else {
            return getFullName();
        }
    }

    public String getXMLName() {
        if (taxa != null) {
            return "tmrca_" + taxa.getId();
        } else {
            return getFullName();
        }
    }

    public String getDescription() {
        if (taxa != null) {
            return description + taxa.getId();
        } else if (prefix != null) {
            return description + " of partition " + prefix;
        }
        return description;
    }

    public double getPriorExpectation() {

        switch (priorType) {
            case LOGNORMAL_PRIOR:
                return LogNormalDistribution.mean(mean, stdev) + offset;
            case NORMAL_PRIOR:
                return NormalDistribution.mean(mean, stdev);
            case EXPONENTIAL_PRIOR:
                return ExponentialDistribution.mean(mean) + offset;
        }
        return 1.0;
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

}
