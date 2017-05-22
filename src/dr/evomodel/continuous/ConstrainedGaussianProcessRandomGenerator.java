/*
 * ConstrainedGaussianProcess.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.continuous;

import dr.evolution.tree.NodeRef;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inference.operators.EllipticalSliceOperator;
import dr.inferencexml.operators.EllipticalSliceOperatorParser;
import dr.math.KroneckerOperation;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 * @author Max R. Tolkoff
 */
public class ConstrainedGaussianProcessRandomGenerator implements GaussianProcessRandomGenerator {

    private final GaussianProcessRandomGenerator generator;
    private final boolean translationInvariant;
    private final boolean rotationInvariant;

    public ConstrainedGaussianProcessRandomGenerator(GaussianProcessRandomGenerator generator,
                                                     boolean translationInvariant, boolean rotationInvariant) {
        this.generator = generator;
        this.translationInvariant = translationInvariant;
        this.rotationInvariant = rotationInvariant;
    }

    @Override
    public Likelihood getLikelihood() {
        throw new RuntimeException("Not yet implemented");
//        return generator.getLikelihood();
    }

    @Override
    public int getDimension() {
        return generator.getDimension();
    }

    @Override
    public double[][] getPrecisionMatrix() {
        throw new RuntimeException("Not yet implemented");
//        return generator.getPrecisionMatrix();
    }

    @Override
    public Object nextRandom() {

        double[] draw = (double[]) generator.nextRandom();
        EllipticalSliceOperator.transformPoint(draw, translationInvariant, rotationInvariant,
                2); // TODO Generalize for dim != 2
        return draw;
    }

    @Override
    public double logPdf(Object x) {
        throw new RuntimeException("Not yet implemented");
//        return generator.logPdf(x);

    }

    public boolean isTranslationInvariant() {
        return translationInvariant;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String CONSTAINED_GAUSSIAN_PROCESS = "constrainedGaussianProcess";
        public static final String TRANSLATION_INVARIANT = EllipticalSliceOperatorParser.TRANSLATION_INVARIANT;
        public static final String ROTATION_INVARIANT = EllipticalSliceOperatorParser.ROTATION_INVARIANT;

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(GaussianProcessRandomGenerator.class),
                AttributeRule.newBooleanRule(TRANSLATION_INVARIANT, true),
                AttributeRule.newBooleanRule(ROTATION_INVARIANT, true),
        };

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean translationInvariant = xo.getAttribute(TRANSLATION_INVARIANT, false);
            boolean rotationInvariant = xo.getAttribute(ROTATION_INVARIANT, false);

            GaussianProcessRandomGenerator generator =
                    (GaussianProcessRandomGenerator) xo.getChild(GaussianProcessRandomGenerator.class);

            return new ConstrainedGaussianProcessRandomGenerator(generator, translationInvariant, rotationInvariant);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Returns a random draw of traits given a trait model and a prior";
        }

        @Override
        public Class getReturnType() {
            return ConstrainedGaussianProcessRandomGenerator.class;
        }

        @Override
        public String getParserName() {
            return CONSTAINED_GAUSSIAN_PROCESS;
        }
    };
}
