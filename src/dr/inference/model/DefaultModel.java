/*
 * DefaultModel.java
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

package dr.inference.model;

import dr.inference.loggers.LogColumn;
import dr.inferencexml.model.DefaultModelParser;

/**
 * @author Marc Suchard
 */
public class DefaultModel extends AbstractModelLikelihood {

    public DefaultModel() {
        super(DefaultModelParser.DUMMY_MODEL);
    }

    public DefaultModel(String str) {
        super(str);
    }

    public DefaultModel(Parameter parameter) {
        super(DefaultModelParser.DUMMY_MODEL);
        addVariable(parameter);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {

    }

    protected void storeState() {

    }

    protected void restoreState() {

    }

    protected void acceptState() {

    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return 0;
    }

    public void makeDirty() {

    }

    public LogColumn[] getColumns() {
        return new LogColumn[0];
    }

}


