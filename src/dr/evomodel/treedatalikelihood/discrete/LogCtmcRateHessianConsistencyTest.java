/*
 * LogCtmcRateHessianConsistencyTest.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.Reportable;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * Compares the two ways LogCtmcRateHessian can compute the empirical Fisher information: the spectral route,
 * which uses the branch log-rate scores, and the generic route through LogAdditiveSubstitutionModelHessian, which
 * is used when the parameter is not the full set of K(K-1) rates. On a K(K-1)-dimensional parameter both routes
 * must agree, whatever the normalization weights are.
 */
public class LogCtmcRateHessianConsistencyTest implements Reportable {

    private final LogCtmcRateHessian spectral;
    private final LogAdditiveSubstitutionModelHessian generic;

    public LogCtmcRateHessianConsistencyTest(LogCtmcRateHessian hessian,
                                             TreeDataLikelihood treeDataLikelihood,
                                             String traitName) {
        final TreeTrait branchLogRateScores =
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchLogRateScoreTraitName(traitName));
        final TreeTrait branchDifferentials =
                treeDataLikelihood.getTreeTrait(SpectralExactGradientDelegate.getBranchDifferentialTraitName(traitName));

        if (branchLogRateScores == null || branchDifferentials == null ||
                hessian.getParameter().getDimension() != hessian.stateCount * (hessian.stateCount - 1)) {
            throw new IllegalArgumentException("The spectral route is not used for this Hessian, " +
                    "so there is nothing to compare the generic route with");
        }

        this.spectral = hessian;
        this.generic = new LogAdditiveSubstitutionModelHessian(hessian, branchDifferentials);
    }

    @Override
    public String getReport() {

        // evaluating the gradient runs the likelihood and the pre-order pass that fill the branch traits
        spectral.getGradientLogDensity();

        final double[] spectralDiagonal = spectral.getDiagonalHessianLogDensity().clone();
        final double[] genericDiagonal = generic.getDiagonalHessianLogDensity().clone();

        final double[][] spectralMatrix = copy(spectral.getHessianLogDensity());
        final double[][] genericMatrix = copy(generic.getHessianLogDensity());

        double maxSpectralDiagonal = 0.0;
        for (double value : spectralDiagonal) {
            maxSpectralDiagonal = Math.max(maxSpectralDiagonal, Math.abs(value));
        }
        if (maxSpectralDiagonal < MIN_INFORMATION) {
            throw new IllegalStateException("The spectral Hessian carries no information (only the ridge), " +
                    "so the comparison with the generic route would be vacuous");
        }

        double maxDiagonalDifference = 0.0;
        for (int i = 0; i < spectralDiagonal.length; ++i) {
            maxDiagonalDifference = Math.max(maxDiagonalDifference, Math.abs(spectralDiagonal[i] - genericDiagonal[i]));
        }

        double maxMatrixDifference = 0.0;
        for (int i = 0; i < spectralMatrix.length; ++i) {
            for (int j = 0; j < spectralMatrix[i].length; ++j) {
                maxMatrixDifference = Math.max(maxMatrixDifference, Math.abs(spectralMatrix[i][j] - genericMatrix[i][j]));
            }
        }

        StringBuilder sb = new StringBuilder();
        sb.append("LogCtmcRateHessian spectral route vs generic route\n");
        sb.append("maxAbsDiagonalDifference: ").append(maxDiagonalDifference).append('\n');
        sb.append("maxAbsMatrixDifference: ").append(maxMatrixDifference).append('\n');
        appendVector(sb, "spectralDiagonal", spectralDiagonal);
        appendVector(sb, "genericDiagonal", genericDiagonal);

        return sb.toString();
    }

    private static double[][] copy(double[][] matrix) {
        final double[][] copy = new double[matrix.length][];
        for (int i = 0; i < matrix.length; ++i) {
            copy[i] = matrix[i].clone();
        }
        return copy;
    }

    private static void appendVector(StringBuilder sb, String label, double[] values) {
        sb.append(label).append(":");
        for (double value : values) {
            sb.append(' ').append(value);
        }
        sb.append('\n');
    }

    private static final double MIN_INFORMATION = 1.0E-6;

    private static final String PARSER_NAME = "logCtmcRateHessianConsistencyTest";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final String traitName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
            final LogCtmcRateHessian hessian = (LogCtmcRateHessian) xo.getChild(LogCtmcRateHessian.class);
            final TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            try {
                return new LogCtmcRateHessianConsistencyTest(hessian, treeDataLikelihood, traitName);
            } catch (IllegalArgumentException e) {
                throw new XMLParseException(e.getMessage());
            }
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                    new ElementRule(LogCtmcRateHessian.class),
                    new ElementRule(TreeDataLikelihood.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Checks that the spectral and generic routes of LogCtmcRateHessian agree.";
        }

        @Override
        public Class getReturnType() {
            return LogCtmcRateHessianConsistencyTest.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }
    };
}
