/*
 * LogAdditiveCtmcRateProvider.java
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

package dr.evomodel.substmodel;

import dr.inference.loggers.LogColumn;
import dr.inference.model.*;

public interface LogAdditiveCtmcRateProvider extends Model, Likelihood {

    double[] getXBeta();

    Parameter getLogRateParameter();

    LogColumn[] getColumns();

    default double[] getRates() {
        double[] rates = getXBeta();
        for (int i = 0; i < rates.length; ++i) {
            rates[i] = Math.exp(rates[i]);
        }
        return rates;
    }

    interface Integrated extends LogAdditiveCtmcRateProvider { }

    interface DataAugmented extends LogAdditiveCtmcRateProvider {

        Parameter getLogRateParameter();

        class Basic extends AbstractModelLikelihood implements DataAugmented {

            private final Parameter logRateParameter;

            public Basic(String name, Parameter logRateParameter) {
                super(name);
                this.logRateParameter = logRateParameter;

                addVariable(logRateParameter);
            }

            public Parameter getLogRateParameter() { return logRateParameter; }

            @Override
            public double[] getXBeta() { // TODO this function should _not_ exponentiate

                final int fieldDim = logRateParameter.getDimension();
                double[] rates = new double[fieldDim];

                for (int i = 0; i < fieldDim; ++i) {
                    rates[i] = Math.exp(logRateParameter.getParameterValue(i));
                }
                return rates;
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) { }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index,
                                                      Parameter.ChangeType type) { }

            @Override
            protected void storeState() { }

            @Override
            protected void restoreState() { }

            @Override
            protected void acceptState() { }

            @Override
            public Model getModel() { return this; }

            @Override
            public double getLogLikelihood() { return 0; }

            @Override
            public void makeDirty() { }
        }
    }
}
