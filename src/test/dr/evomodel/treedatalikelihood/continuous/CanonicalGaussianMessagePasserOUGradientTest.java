/*
 * CanonicalGaussianMessagePasserOUGradientTest.java
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

package test.dr.evomodel.treedatalikelihood.continuous;

/**
 * Validates canonical OU tree-passer gradients against numerical finite differences.
 */
public class CanonicalGaussianMessagePasserOUGradientTest extends CanonicalOUMessagePasserTestSupport {

    public CanonicalGaussianMessagePasserOUGradientTest(final String name) {
        super(name);
    }

    public void testGradientQConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (conjugate root, fully observed)");
        checkGradientQ(buildOUSetup("canonGradQ_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientAConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂A analytic vs numerical (conjugate root, fully observed)");
        checkGradientA(buildOUSetup("canonGradA_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientMuConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂mu analytic vs numerical (conjugate root, fully observed)");
        checkGradientMu(buildOUSetup("canonGradMu_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientBranchLengthsConjugateRootFullyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂t analytic vs numerical (conjugate root, fully observed)");
        checkGradientBranchLengths(buildOUSetup("canonGradT_full", buildFullyObservedTips()), "fully observed");
    }

    public void testGradientQConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (conjugate root, partially observed)");
        checkGradientQ(buildOUSetup("canonGradQ_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientQFixedRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂Q analytic vs numerical (fixed root, partially observed)");
        checkGradientQ(buildOUSetup("canonGradQ_fixed_partial", buildPartiallyObservedTips(), rootPriorInf),
                "fixed root, partially observed");
    }

    public void testGradientAConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂A analytic vs numerical (conjugate root, partially observed)");
        checkGradientA(buildOUSetup("canonGradA_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientMuConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂mu analytic vs numerical (conjugate root, partially observed)");
        checkGradientMu(buildOUSetup("canonGradMu_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testGradientBranchLengthsConjugateRootPartiallyObservedMatchesNumerical() {
        System.out.println("\nTest: canonical OU ∂logL/∂t analytic vs numerical (conjugate root, partially observed)");
        checkGradientBranchLengths(buildOUSetup("canonGradT_partial", buildPartiallyObservedTips()), "partially observed");
    }

    public void testOrthogonalBlockConjugateRootFullyObservedNativeSelectionGradientMatchesNumerical() {
        System.out.println("\nTest: canonical OU orthogonal-block native ∂logL/∂A matches numerical (fully observed)");
        checkOrthogonalBlockNativeSelectionGradient(
                buildOrthogonalBlockOUSetup("canonOrthGrad_full", buildFullyObservedTips()),
                "fully observed");
    }

    public void testOrthogonalBlockConjugateRootPartiallyObservedNativeSelectionGradientMatchesNumerical() {
        System.out.println("\nTest: canonical OU orthogonal-block native ∂logL/∂A matches numerical (partially observed)");
        checkOrthogonalBlockNativeSelectionGradient(
                buildOrthogonalBlockOUSetup("canonOrthGrad_partial", buildPartiallyObservedTips()),
                "partially observed");
    }

    public void testOrthogonalBlockConjugateRootPartiallyObservedGradientQMatchesDense() {
        System.out.println("\nTest: canonical OU orthogonal-block ∂logL/∂Q matches dense (partially observed)");
        checkOrthogonalBlockGradientQMatchesDense(
                buildOrthogonalBlockOUSetup("canonOrthGradQ_partial", buildPartiallyObservedTips()),
                "partially observed");
    }
}
