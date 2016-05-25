/*
 * Layer.java
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

package dr.app.mapper.application.mapper;

import dr.inference.trace.Trace;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class Layer {
    public Layer(String name, Trace trace, LayerType type, List<String> labels) {
        this.name = name;
        this.trace = trace;
        this.type = type;
        this.labels = labels;
    }

    public String getName() {
        return name;
    }

    public Trace getTrace() {
        return trace;
    }

    public LayerType getType() {
        return type;
    }

    public List<String> getLabels() {
        return labels;
    }

    private String name;
    private Trace trace;
    private LayerType type;

    private List<String> labels;
}
