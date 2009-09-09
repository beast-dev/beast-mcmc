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


/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Walter Xie
 */
public abstract class PartitionOptions extends ModelOptions {   

    abstract public String getPrefix(); 
    
    abstract public Class<?> getPartitionClassType(); 
    
    public Parameter createParameter(PartitionOptions options, String name, String description, PriorScaleType scale,
            double value, double lower, double upper) {
        final Parameter parameter = new Parameter.Builder(name, description).scaleType(scale).prior(PriorType.UNIFORM_PRIOR).initial(value)        
                .lower(lower).upper(upper).partitionOptions(options).build();
        parameters.put(name, parameter);
        return parameter;
    }

    public void createParameter(PartitionOptions options, String name, String description, boolean isNodeHeight, double value,
            double lower, double upper) {
        final Parameter parameter = new Parameter.Builder(name, description).isNodeHeight(isNodeHeight).scaleType(PriorScaleType.TIME_SCALE)
                .initial(value).lower(lower).upper(upper).partitionOptions(options).build();
        parameters.put(name, parameter);
    }
}
