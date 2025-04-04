/*
 * PartitionOptions.java
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

package dr.app.beauti.options;

import dr.app.beauti.types.PriorScaleType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class PartitionOptions extends ModelOptions {

    private String name;
    protected final BeautiOptions options;

//    protected double[] avgRootAndRate = new double[]{1.0, 1.0};

    public PartitionOptions(BeautiOptions options) {
        this.options = options;
    }

    public PartitionOptions(BeautiOptions options, String name) {
        this.options = options;
        this.name = name;
    }

    protected void createParameterTree(PartitionOptions options, String name, String description, boolean isNodeHeight) {
        new Parameter.Builder(name, description).isNodeHeight(isNodeHeight).scaleType(PriorScaleType.TIME_SCALE)
                .isNonNegative(true).initial(Double.NaN).partitionOptions(options).build(parameters);
    }

    public Parameter getParameter(String name) {

        Parameter parameter = parameters.get(name);

        if (parameter == null) {
            throw new IllegalArgumentException("Parameter with name, " + name + ", is unknown");
        }

        parameter.setPrefix(getPrefix());

//        autoScale(parameter); // not include clock rate, and treeModel.rootHeight

        return parameter;
    }

    public boolean hasParameter(String name) {
        return parameters.get(name) != null;
    }

    public Operator getOperator(String name) {

        Operator operator = operators.get(name);

        if (operator == null) throw new IllegalArgumentException("Operator with name, " + name + ", is unknown");

        operator.setPrefix(getPrefix());

        return operator;
    }

    public String getName() {
        return name;

    }

    public void setName(String name) {
        this.name = name;
    }

    public String toString() {
        return getName();
    }

    public DataType getDataType() {
        if (options.getDataPartitions(this).size() == 0) {
            return Nucleotides.INSTANCE;
        }
        return options.getDataPartitions(this).get(0).getDataType();
    }

    public BeautiOptions getOptions() {
        return options;
    }

}
