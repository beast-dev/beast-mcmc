/*
 * Operator.java
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

import dr.app.beauti.types.OperatorType;

import java.io.Serializable;
import java.util.Map;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class Operator implements Serializable {

    private static final long serialVersionUID = -1165043783660638155L;
    private String prefix = null;

    // final
    private String baseName;
    private final String description;
    private final OperatorType operatorType;
    private Parameter parameter1;
    private Parameter parameter2;
    private final PartitionOptions options;
    private final String tag;

    // editable
    private double tuning;
    private double weight;
    private boolean tuningEdited = false;
    private boolean isUsed = true;

    private boolean autoOptimize = true;

    private String idref;

    public static class Builder {
        // Required para
        private final String baseName;
        private final String description;
        private final Parameter parameter1;

        private final OperatorType operatorType;
        private final double tuning;
        private final double weight;

        // Optional para - initialized to default values
        private Parameter parameter2 = null;
        private PartitionOptions options = null;
        private String tag = null;
        private String idref = null;

        private boolean autoOptimize = true;

        public Builder(String name, String description, Parameter parameter, OperatorType type, double tuning, double weight) {
            this.baseName = name;
            this.description = description;
            this.parameter1 = parameter;
            this.operatorType = type;
            this.tuning = tuning;
            this.weight = weight;
        }


        public Builder parameter2(Parameter parameter2) {
            this.parameter2 = parameter2;
            return this;
        }

        public Builder partitionOptions(PartitionOptions options) {
            this.options = options;
            return this;
        }

        public Builder tag(String tag) {
            this.tag = tag;
            return this;
        }

        public Builder idref(String idref) {
            this.idref = idref;
            return this;
        }

        public Builder autoOptimize(boolean autoOptimize) {
            this.autoOptimize = autoOptimize;
            return this;
        }

        public Operator build() {
            return new Operator(this);
        }

        public Operator build(Map<String, Operator> map) {
            final Operator operator = new Operator(this);
            map.put(baseName, operator);
            return operator;
        }
    }

    private Operator(Builder builder) {
        baseName = builder.baseName;
        description = builder.description;
        parameter1 = builder.parameter1;
        operatorType = builder.operatorType;
        tuning = builder.tuning;
        weight = builder.weight;
        parameter2 = builder.parameter2;
        options = builder.options;
        tag = builder.tag;
        idref = builder.idref;
        isUsed = true;
        tuningEdited = false;
        autoOptimize = builder.autoOptimize;
    }

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    public String getDescription() {
        if (description == null || description.length() == 0) {
            String prefix = "";
            if (operatorType == OperatorType.SCALE || operatorType == OperatorType.SCALE_ALL) {
                prefix = "Scales the ";
            } else if (operatorType == OperatorType.RANDOM_WALK) {
                prefix = "A random-walk on the ";
            }
            return prefix + parameter1.getDescription();
        } else {
            return description;
        }
    }

    public boolean isTunable() {
        return tuning > 0;
    }

    public String getTag() {
        return tag;
    }

    public String getIdref() {
        return idref;
    }

    public double getWeight() {
        return weight;
    }

    public void setWeight(double weight) {
        this.weight = weight;
    }

    public double getTuning() {
        return tuning;
    }

    public void setTuning(double tuning) {
        this.tuning = tuning;
        tuningEdited = true;
    }

    public boolean isTuningEdited() {
        return tuningEdited;
    }

    public OperatorType getOperatorType() {
        return operatorType;
    }

    public boolean isUsed() {
        return isParameterFixed() ? false : isUsed;
    }

    public void setUsed(boolean used) {
        this.isUsed = used;
    }

    public boolean isAutoOptimize() {
        return autoOptimize;
    }

    public boolean isParameterFixed() {
        return parameter1 != null && parameter1.isFixed();
    }

    public Parameter getParameter1() {
        return parameter1;
    }

    public void setParameter1(Parameter parameter1) {
        this.parameter1 = parameter1;
    }

    public Parameter getParameter2() {
        return parameter2;
    }

    public void setParameter2(Parameter parameter2) {
        this.parameter2 = parameter2;
    }

    public void setPrefix(String prefix) {
        this.prefix = prefix;
    }

    public PartitionOptions getOptions() {
        return options;
    }

    public String getPrefix() {
        return prefix;
    }

    public String getName() {
        String name = baseName;
        if (prefix != null) {
            name = prefix + baseName;
        }
        return name;
    }

    public void setName(String name) {
        this.baseName = name;
    }

    public String getBaseName() {
        return baseName;
    }

}
