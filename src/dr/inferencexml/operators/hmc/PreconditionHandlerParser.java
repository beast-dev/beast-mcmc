/*
 * PreconditionHandlerParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.operators.hmc;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Parameter;
import dr.inference.model.PriorPreconditioningProvider;
import dr.inference.operators.hmc.MassPreconditionScheduler;
import dr.inference.operators.hmc.MassPreconditioner;
import dr.inference.operators.hmc.MassPreconditioningOptions;
import dr.inference.operators.hmc.PreconditionHandler;
import dr.util.Transform;
import dr.xml.*;

import static dr.util.Transform.Util.parseTransform;

/**
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class PreconditionHandlerParser extends AbstractXMLObjectParser {

    private final static String PRECONDITIONING = "preconditioning";
    private final static String PRECONDITIONING_SCHEDULE = "preconditioningSchedule";
    private final static String PRECONDITIONING_UPDATE_FREQUENCY = "preconditioningUpdateFrequency";
    final static String PRECONDITIONING_MAX_UPDATE = "preconditioningMaxUpdate";
    final static String PRECONDITIONING_DELAY = "preconditioningDelay";
    private final static String PRECONDITIONING_MEMORY = "preconditioningMemory";
    private final static String PRECONDITIONER = "preconditioner";
    private final static String PRECONDITIONING_GUESS_INIT_MASS = "guessInitialMass";
    private final static String BOUNDS = "bounds";
    private final static String CONSTANT = "constant";

    public static PreconditionHandler parsePreconditionHandler(XMLObject xo) throws XMLParseException {
        MassPreconditioner.Type preconditioningType;
        if (xo.hasChildNamed(PRECONDITIONER)) {
            preconditioningType = MassPreconditioner.Type.PRIOR_DIAGONAL;
        } else {
            preconditioningType = parsePreconditioning(xo);
        }

        MassPreconditionScheduler.Type preconditionSchedulerType = parsePreconditionScheduler(xo, preconditioningType);

        int preconditioningUpdateFrequency = xo.getAttribute(PRECONDITIONING_UPDATE_FREQUENCY, 0);

        int preconditioningMaxUpdate = xo.getAttribute(PRECONDITIONING_MAX_UPDATE, 0);

        int preconditioningDelay = xo.getAttribute(PRECONDITIONING_DELAY, 0);

        int preconditioningMemory = xo.getAttribute(PRECONDITIONING_MEMORY, 0);

        boolean guessInitialMass = xo.getAttribute(PRECONDITIONING_GUESS_INIT_MASS, false);

        Transform transform = parseTransform(xo);

        GradientWrtParameterProvider derivative =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);

        Parameter eigenLowerBound, eigenUpperBound;
        if (xo.hasChildNamed(BOUNDS)) {
            eigenLowerBound = xo.getChild(BOUNDS).getAllChildren(Parameter.class).get(0);
            eigenUpperBound = xo.getChild(BOUNDS).getAllChildren(Parameter.class).get(1);
        } else {
            eigenLowerBound = new Parameter.Default(1E-2);
            eigenUpperBound = new Parameter.Default(1E2);
        }

        Parameter preconditioningAddedConstant;
        if (xo.hasChildNamed(CONSTANT)) {
            preconditioningAddedConstant = (Parameter) xo.getChild(CONSTANT).getChild(Parameter.class);
        } else {
            preconditioningAddedConstant = new Parameter.Default(0.0);
        }

        MassPreconditioningOptions preconditioningOptions =
                new MassPreconditioningOptions.Default(preconditioningUpdateFrequency, preconditioningMaxUpdate,
                        preconditioningDelay, preconditioningMemory, guessInitialMass, eigenLowerBound, eigenUpperBound, preconditioningAddedConstant);

        MassPreconditioner preconditioner;

        if (xo.hasChildNamed(PRECONDITIONER)) {
            PriorPreconditioningProvider priorPreconditioningProvider = (PriorPreconditioningProvider) xo.getChild(PRECONDITIONER).getChild(PriorPreconditioningProvider.class);
            if (priorPreconditioningProvider !=  null) {
                preconditioner = new MassPreconditioner.PriorPreconditioner(priorPreconditioningProvider, transform);
            } else {
                throw new XMLParseException("Unknown preconditioner specified");
            }
        } else {
            preconditioner = preconditioningType.factory(derivative, transform, preconditioningOptions);
        }
        PreconditionHandler preconditionHandler = new PreconditionHandler(preconditioner, preconditioningOptions, preconditionSchedulerType);

        return preconditionHandler;
    }
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return parsePreconditionHandler(xo);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {

            new ElementRule(Transform.MultivariableTransformWithParameter.class, true),

            new XORRule(
                    new ElementRule(GradientWrtParameterProvider.class),
                    new ElementRule(PRECONDITIONER, new XMLSyntaxRule[]{
                            new XORRule(
                                    new ElementRule(MassPreconditioner.class),
                                    new ElementRule(PriorPreconditioningProvider.class)
                            ),
                    })
            ),

            new ElementRule(BOUNDS, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class, 2, 2)
            }, true),

            new ElementRule(CONSTANT, Parameter.class, "A small constant added to the diagonal of the preconditioning mass matrix."),

    };

    static MassPreconditioner.Type parsePreconditioning(XMLObject xo) throws XMLParseException {

        return MassPreconditioner.Type.parseFromString(
                xo.getAttribute(PRECONDITIONING, MassPreconditioner.Type.NONE.getName())
        );
    }

    static MassPreconditionScheduler.Type parsePreconditionScheduler(XMLObject xo,
                                                                     MassPreconditioner.Type preconditioningType) throws XMLParseException {
        if (preconditioningType == MassPreconditioner.Type.NONE) {
            return MassPreconditionScheduler.Type.NONE;
        } else {
            return MassPreconditionScheduler.Type.parseFromString(
                    xo.getAttribute(PRECONDITIONING_SCHEDULE, MassPreconditionScheduler.Type.DEFAULT.getName()));
        }
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return PreconditionHandler.class;
    }

    @Override
    public String getParserName() {
        return PRECONDITIONING;
    }
}
