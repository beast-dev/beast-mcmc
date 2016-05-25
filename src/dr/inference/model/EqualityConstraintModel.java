/*
 * EqualityConstraintModel.java
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

import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Andrew Rambaut
 */
public class EqualityConstraintModel extends AbstractModel {

    public static final String EQUALITY_CONSTRAINT = "equalityConstraint";

    public EqualityConstraintModel(String name, List<Variable> parameterList) {
        super(name);
        this.parameterList = parameterList;
        for (Variable p : parameterList) {
            addVariable(p);
        }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // Do nothing
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        if (noReentry) {
            return; // Already processing a member variable
        }

        noReentry = true;
        if (type == Variable.ChangeType.ALL_VALUES_CHANGED) {
            Object[] newValues = variable.getValues();
            for (Variable p : parameterList) {
                if (p != variable) {
                    for (int i = 0; i < newValues.length; ++i) {
                        p.setValue(i, newValues[i]);
                    }
                }
            }
        } else if (type == Variable.ChangeType.VALUE_CHANGED) {
            Object newValue = variable.getValue(index);
            for (Variable p : parameterList) {
                if (p != variable) {
                    p.setValue(index, newValue);
                }
            }
        } else {
            throw new IllegalArgumentException("Variable dimensional updates are not yet implemented.");
        }
        noReentry = false;
    }

    protected void storeState() {
        // Do nothing
    }

    protected void restoreState() {
        // Do nothing
    }

    protected void acceptState() {
        // Do nothing
    }

    private static void setEqual(List<Variable> parameterList) {

        final int dim = parameterList.get(0).getSize();


        if (parameterList.get(0) instanceof Parameter) {
            // Can take an average
            double[] total = new double[dim];

            for (Variable v : parameterList) {
                Parameter p = (Parameter) v;
                for (int j = 0; j < dim; ++j) {
                    total[j] += p.getParameterValue(j);
                }
            }

            for (int j = 0; j < dim; ++j) {
                total[j] /= parameterList.size();
            }

            for (Variable v : parameterList) {
                Parameter p = (Parameter) v;
                for (int j = 0; j < dim; ++j) {
                    p.setParameterValue(j, total[j]);
                }
            }
        }
    }

    private static int checkDimensions(List<Parameter> parameterList) {
        int dim = parameterList.get(0).getDimension();
        for (Parameter p : parameterList) {
            if (p.getDimension() != dim) {
                return -1;
            }
        }
        return dim;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int dim = -1;
            double[] firstValues = null;
            List<Variable> parameterList = new ArrayList<Variable>();

            for (int i = 0; i < xo.getChildCount(); ++i) {
                Parameter param = (Parameter) xo.getChild(i);
                if (i == 0) {
                    dim = param.getDimension();
                    firstValues = param.getParameterValues();
                } else {
                    if (param.getDimension() != dim) {
                        throw new XMLParseException("All parameters must have the same dimension.");
                    }
                    for (int j = 0; j < dim; ++j) {
                        param.setParameterValue(j, firstValues[j]);
                    }
                }
                parameterList.add(param);
            }
            return new EqualityConstraintModel("equalityConstraint", parameterList);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(Parameter.class, 2, Integer.MAX_VALUE),
        };

        public String getParserDescription() {
            return "Forces a set of variables to have equal values";
        }

        public Class getReturnType() {
            return EqualityConstraintModel.class;
        }

        public String getParserName() {
            return EQUALITY_CONSTRAINT;
        }
    };

    final private List<Variable> parameterList;
    private boolean noReentry = false;
}
