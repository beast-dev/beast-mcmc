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

import dr.app.beauti.types.OperatorType;

import java.io.Serializable;
import java.util.Map;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public class Operator implements Serializable {

    private String prefix = null;

    // final
    private final String baseName;
    private final String description;
    public final OperatorType operatorType;
    public final Parameter parameter1;
    public final Parameter parameter2;
    private final PartitionOptions options;
    public final String tag;

    // editable
    public double tuning;
    public double weight;
    public boolean tuningEdited;
    public boolean inUse;
    public String idref;

    private ClockModelGroup clockModelGroup = null;

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

        private boolean inUse = true;
        private boolean tuningEdited = false;

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

        public Builder isInUse(boolean inUse) {
            this.inUse = inUse;
            return this;
        }

        public Builder tuningEdited(boolean tuningEdited) {
            this.tuningEdited = tuningEdited;
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
        inUse = builder.inUse;
        tuningEdited = builder.tuningEdited;
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

    public String getBaseName() {
        return baseName;
    }

    public ClockModelGroup getClockModelGroup() {
        return clockModelGroup;
    }

    public void setClockModelGroup(ClockModelGroup clockModelGroup) {
        this.clockModelGroup = clockModelGroup;
    }

}
