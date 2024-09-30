/*
 * LatentFactorModelSimulator.java
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

package dr.evomodel.continuous;

import dr.inference.model.FastMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class LatentFactorModelSimulator {

    public static final String SIMULATE_LATENT_FACTOR_MODEL = "simulateLatentFactorModel";
    public static final String PRECISION = "precision";


    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return SIMULATE_LATENT_FACTOR_MODEL; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            final String name = xo.hasId() ? xo.getId() : null;

            FullyConjugateMultivariateTraitLikelihood treeSample = (FullyConjugateMultivariateTraitLikelihood)xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
            MatrixParameterInterface Loadings = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);

            GaussianProcessFromTree randomTips = new GaussianProcessFromTree(treeSample);

            int ntaxa = treeSample.getTreeModel().getExternalNodeCount();
            int ntrait = Loadings.getRowDimension();

            FastMatrixParameter data = new FastMatrixParameter(name, ntrait, ntaxa, 0);

            if(treeSample.getNumData() * treeSample.getDimTrait() != Loadings.getColumnDimension()){
                throw new RuntimeException("Number of factors in tree and loadings matrix must be identical:\n Factor Matrix: "
                        + treeSample.getNumData() * treeSample.getDimTrait() +
                        "\n Loadings Matrix: "
                        + Loadings.getColumnDimension()
                );
            }

            Parameter precision = (Parameter) xo.getChild(PRECISION).getChild(Parameter.class);
            int nfac = Loadings.getColumnDimension();

            double[] Factors = randomTips.nextRandomFast();

            for (int i = 0; i < ntrait; i++) {
                for (int j = 0; j < ntaxa; j++) {
                    double sum = 0;
                    for (int k = 0; k < nfac; k++) {
                        sum += Factors[j * nfac + k] * Loadings.getParameterValue(i, k);
                    }
                    double draw = (Double)
                            (new NormalDistribution(sum, 1 / Math.sqrt(precision.getParameterValue(i)))).nextRandom();
                    data.setParameterValue(i, j, draw);
                }

            }

            return data;
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
                new ElementRule(MatrixParameterInterface.class),
                new ElementRule(PRECISION,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}
                        )
        };

        public String getParserDescription() {
            return "Simulates a data matrix for a Latent Factor Model";
        }

        public Class getReturnType() { return MatrixParameter.class; }
    };

}
