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

import dr.app.beauti.enumTypes.OperatorType;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class Operator {

    public Operator(String name,
                    String description,
                    Parameter parameter,
                    OperatorType operatorType,
                    double tuning,
                    double weight) {

        this.baseName = name;
        this.description = description;
        this.parameter1 = parameter;
        this.parameter2 = null;
        this.tag = null;

        this.type = operatorType;
        this.tuningEdited = false;
        this.tuning = tuning;
        this.weight = weight;

        this.inUse = true;
    }

    public Operator(String name, String description,
                    Parameter parameter1, Parameter parameter2,
                    OperatorType operatorType, double tuning, double weight) {
        this.baseName = name;
        this.description = description;
        this.parameter1 = parameter1;
        this.parameter2 = parameter2;
        this.tag = null;
        
        this.type = operatorType;
        this.tuningEdited = false;
        this.tuning = tuning;
        this.weight = weight;

        this.inUse = true;
    }
    
    public Operator(String name, String description,
		            Parameter parameter, String tag, String idref,
		            OperatorType operatorType, double tuning, double weight) {
		this.baseName = name;
		this.description = description;
		this.parameter1 = parameter;
		this.parameter2 = null;	
		
		this.tag = tag;
		this.idref = idref;
		
		this.type = operatorType;
		this.tuningEdited = false;
		this.tuning = tuning;
		this.weight = weight;
		
		this.inUse = true;
	}

    public String getDescription() {
        if (description == null || description.length() == 0) {
            String prefix = "";
            if (type == OperatorType.SCALE || type == OperatorType.SCALE_ALL) {
                prefix = "Scales the ";
            } else if (type == OperatorType.RANDOM_WALK) {
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

    private final String baseName;
    public String prefix = null;

    public final String description;

    public final OperatorType type;
    public boolean tuningEdited;
    public double tuning;
    public double weight;
    public boolean inUse;
    
    public final String tag;
    public String idref;
    
    public final Parameter parameter1;
    public final Parameter parameter2;    
}
