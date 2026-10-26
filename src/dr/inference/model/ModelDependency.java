/*
 * ModelDependency.java
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

package dr.inference.model;

import java.util.ArrayList;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Set;

/**
 * Finds the models below a root model that hold a given parameter directly as a variable.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class ModelDependency {

    private ModelDependency() {
    }

    /**
     * @return every model reachable from root (depth first) with a variable that is parameter, or a
     * CompoundParameter containing it
     */
    public static List<Model> findHolders(Model root, Parameter parameter) {
        List<Model> holders = new ArrayList<>();
        Set<Model> visited = Collections.newSetFromMap(new IdentityHashMap<Model, Boolean>());
        visit(root, parameter, visited, holders);
        return holders;
    }

    private static void visit(Model model, Parameter parameter, Set<Model> visited, List<Model> holders) {
        if (model == null || !visited.add(model)) {
            return;
        }
        for (int i = 0; i < model.getVariableCount(); ++i) {
            if (holds(model.getVariable(i), parameter)) {
                holders.add(model);
                break;
            }
        }
        for (int i = 0; i < model.getModelCount(); ++i) {
            visit(model.getModel(i), parameter, visited, holders);
        }
    }

    private static boolean holds(Variable variable, Parameter parameter) {
        if (variable == parameter) {
            return true;
        }
        if (variable instanceof CompoundParameter) {
            CompoundParameter compound = (CompoundParameter) variable;
            for (int i = 0; i < compound.getParameterCount(); ++i) {
                if (holds(compound.getParameter(i), parameter)) {
                    return true;
                }
            }
        }
        return false;
    }
}
