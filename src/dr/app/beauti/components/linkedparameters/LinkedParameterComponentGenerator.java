/*
 * LinkedParameterComponentGenerator.java
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

package dr.app.beauti.components.linkedparameters;

import dr.app.beauti.generator.BaseComponentGenerator;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.Parameter;
import dr.app.beauti.util.XMLWriter;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 * @version $Id$
 */
public class LinkedParameterComponentGenerator extends BaseComponentGenerator {

    public LinkedParameterComponentGenerator(final BeautiOptions options) {
        super(options);
    }

    public boolean usesInsertionPoint(final InsertionPoint point) {
        LinkedParameterComponentOptions comp = (LinkedParameterComponentOptions)
                options.getComponentOptions(LinkedParameterComponentOptions.class);

        switch (point) {
            case BEFORE_OPERATORS:
            case IN_FILE_LOG_PARAMETERS:
                return !comp.isEmpty();
            default:
                return false;
        }
    }

    protected void generate(final InsertionPoint point, final Object item, final String prefix, final XMLWriter writer) {
        LinkedParameterComponentOptions comp = (LinkedParameterComponentOptions)
                options.getComponentOptions(LinkedParameterComponentOptions.class);

        switch (point) {
            case BEFORE_OPERATORS:
                for (LinkedParameter linkedParameter : comp.getLinkedParameterList()) {
                    generateJointParameter(linkedParameter, comp.getDependentParameters(linkedParameter), writer);
                }
                break;
            case IN_FILE_LOG_PARAMETERS:
                for (LinkedParameter linkedParameter : comp.getLinkedParameterList()) {
                    writer.writeIDref(ParameterParser.PARAMETER, linkedParameter.getName());
                }
                break;
            default:
                throw new IllegalArgumentException("This insertion point is not implemented for " + this.getClass().getName());
        }
    }

    private void generateJointParameter(LinkedParameter linkedParameter, List<Parameter> parameters, XMLWriter writer) {

        writer.writeOpenTag("jointParameter",
                new Attribute[]{
                        new Attribute.Default<String>(XMLParser.ID, linkedParameter.getName())
                }
        );

        for (Parameter parameter : parameters) {
            writer.writeTag(ParameterParser.PARAMETER, new Attribute.Default<String>(XMLParser.IDREF, parameter.getName()), true);
        }

        writer.writeCloseTag("jointParameter");
    }


    protected String getCommentLabel() {
        return "Linked parameter";
    }

}