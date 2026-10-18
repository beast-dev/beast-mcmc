/*
 * MultiTypeBirthDeathLikelihoodTest.java
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

package test.dr.evomodel.speciation;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.util.Units;
import dr.evomodel.speciation.EfficientSpeciationLikelihood;
import dr.evomodel.speciation.MultiTypeBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * Tests for the PROTOTYPE {@link MultiTypeBirthDeathSerialSamplingModel}.
 * <p>
 * "testReducesToSingleType" is the main correctness check: with numTypes=1, the multi-type
 * model's log-likelihood must match {@link NewBirthDeathSerialSamplingModel}'s closed-form
 * answer on the same tree/parameters, i.e. the general K machinery must exactly reproduce
 * the existing exact K=1 model as a special case.
 * <p>
 * "testUnreachableSecondTypeDilutesConditioning" documents (and checks the size of) a
 * subtlety noted in the model's class-level Javadoc: declaring a second type that the data
 * can never visit still changes the conditionOnSurvival normalization (it has its own
 * survival probability, diluting the denominator), even though it leaves the unconditioned
 * likelihood exactly unchanged. This was caught as an apparent "bug" while building this
 * prototype and turned out to be a real, if slightly unsatisfying, modeling consequence of
 * the "uniform sum over declared types" origin-marginalization convention documented on the
 * class -- see the design note "multitype-BDS-fast-approximations.md" and the class Javadoc.
 * <p>
 * "testTwoTypeSanity" exercises real two-type dynamics (migration and cross-type birth) and
 * only checks internal sanity (finite output), since no independent ground truth is
 * available for that case yet.
 *
 * @author Marc A. Suchard
 * @author Claude
 */
public class MultiTypeBirthDeathLikelihoodTest extends TestCase {

    private static final String NEWICK = "((1:1.0,2:1.0):1.0,3:2.0);";
    private static final double LAMBDA = 1.5;
    private static final double MU = 0.4;
    private static final double PSI = 0.3;
    private static final double R = 1.0;
    private static final double RHO = 0.2;
    private static final double ORIGIN = 4.0;

    public MultiTypeBirthDeathLikelihoodTest(String name) {
        super(name);
    }

    private DefaultTreeModel makeTree(String[] tipTaxonIds, String[] tipTypes) throws Exception {
        NewickImporter importer = new NewickImporter(NEWICK);
        FlexibleTree flexTree = (FlexibleTree) importer.importTree(null);
        for (int i = 0; i < flexTree.getExternalNodeCount(); ++i) {
            NodeRef node = flexTree.getExternalNode(i);
            String id = flexTree.getNodeTaxon(node).getId();
            for (int k = 0; k < tipTaxonIds.length; ++k) {
                if (tipTaxonIds[k].equals(id)) {
                    flexTree.getNodeTaxon(node).setAttribute("state", tipTypes[k]);
                }
            }
        }
        return new DefaultTreeModel(flexTree);
    }

    private double referenceLogL(DefaultTreeModel tree, boolean condition) {
        NewBirthDeathSerialSamplingModel refModel = new NewBirthDeathSerialSamplingModel(
                new Parameter.Default(LAMBDA), new Parameter.Default(MU), new Parameter.Default(PSI),
                new Parameter.Default(R), new Parameter.Default(RHO), new Parameter.Default(ORIGIN),
                condition, 1, Double.POSITIVE_INFINITY, Units.Type.YEARS);
        refModel.setupTimeline(null);
        // The reference model's likelihood lives behind the interval-processing (Efficient)
        // machinery; a plain SpeciationLikelihood cannot drive it (calculateTreeLogLikelihood
        // is intentionally unimplemented on NewBirthDeathSerialSamplingModel itself).
        EfficientSpeciationLikelihood refLikelihood = new EfficientSpeciationLikelihood(tree, refModel, null, "ref");
        return refLikelihood.getLogLikelihood();
    }

    public void testReducesToSingleType() throws Exception {

        // tip "3" is at height 0 (present-day sample); tips "1"/"2" coalesce at height 1.
        DefaultTreeModel tree = makeTree(new String[]{"1", "2", "3"}, new String[]{"A", "A", "A"});

        double refLogL = referenceLogL(tree, true);

        MultiTypeBirthDeathSerialSamplingModel mtModel = new MultiTypeBirthDeathSerialSamplingModel(
                "mtbd1", 1, new String[]{"A"},
                new Parameter.Default(new double[]{LAMBDA}),
                new Parameter.Default(new double[]{0.0}),
                new Parameter.Default(new double[]{MU}),
                new Parameter.Default(new double[]{PSI}),
                new Parameter.Default(new double[]{R}),
                new Parameter.Default(new double[]{RHO}),
                new Parameter.Default(ORIGIN),
                true, 1, Double.POSITIVE_INFINITY, "state", 4000, Units.Type.YEARS);

        double mtLogL = new SpeciationLikelihood(tree, mtModel, "mt1").getLogLikelihood();

        System.out.println("Reference (K=1 exact) logL = " + refLogL);
        System.out.println("Multi-type (K=1)      logL = " + mtLogL);

        // RK4 at 4000 steps/unit time vs. the exact closed form: expect agreement to ~1e-6.
        assertEquals(refLogL, mtLogL, 1e-6);
    }

    public void testUnreachableSecondTypeDilutesConditioning() throws Exception {

        DefaultTreeModel tree = makeTree(new String[]{"1", "2", "3"}, new String[]{"A", "A", "A"});

        // K=2, but type "B" is unreachable from an all-type-A tree: no B->A or A->B birth,
        // no migration. B still has its own (nonzero) death/sampling rates, hence its own
        // nonzero survival probability 1-E_B(origin).
        Parameter birthRate = new Parameter.Default(new double[]{LAMBDA, 0.0, 0.0, 0.0});
        Parameter migrationRate = new Parameter.Default(new double[]{0.0, 0.0, 0.0, 0.0});
        Parameter deathRate = new Parameter.Default(new double[]{MU, MU});
        Parameter samplingRate = new Parameter.Default(new double[]{PSI, PSI});
        Parameter treatment = new Parameter.Default(new double[]{R, R});
        Parameter samplingProb = new Parameter.Default(new double[]{RHO, RHO});
        Parameter originParam = new Parameter.Default(ORIGIN);
        String[] typeLabels = {"A", "B"};

        // Unconditioned (conditionOnSurvival=false): type B never contributes to the sum
        // over types at any node (it can never appear), so this MUST match the K=1 model
        // exactly regardless of the extra, unreachable type.
        double refUnconditioned = referenceLogL(tree, false);
        MultiTypeBirthDeathSerialSamplingModel mtModelUnconditioned = new MultiTypeBirthDeathSerialSamplingModel(
                "mtbdA", 2, typeLabels, birthRate, migrationRate, deathRate, samplingRate,
                treatment, samplingProb, originParam, false, 1, Double.POSITIVE_INFINITY,
                "state", 4000, Units.Type.YEARS);
        double mtUnconditioned = new SpeciationLikelihood(tree, mtModelUnconditioned, "mtA").getLogLikelihood();

        System.out.println("Unconditioned: reference = " + refUnconditioned + ", multi-type (K=2) = " + mtUnconditioned);
        assertEquals(refUnconditioned, mtUnconditioned, 1e-6);

        // Conditioned (conditionOnSurvival=true): now they should legitimately DIFFER, because
        // this model's conditioning normalizer sums (1-E_k(origin)) over BOTH declared types
        // with equal weight, and type B's own survival probability dilutes the denominator
        // even though B can never explain this tree. Documented in the class Javadoc.
        double refConditioned = referenceLogL(tree, true);
        MultiTypeBirthDeathSerialSamplingModel mtModelConditioned = new MultiTypeBirthDeathSerialSamplingModel(
                "mtbdB", 2, typeLabels, birthRate, migrationRate, deathRate, samplingRate,
                treatment, samplingProb, originParam, true, 1, Double.POSITIVE_INFINITY,
                "state", 4000, Units.Type.YEARS);
        double mtConditioned = new SpeciationLikelihood(tree, mtModelConditioned, "mtB").getLogLikelihood();

        System.out.println("Conditioned:   reference = " + refConditioned + ", multi-type (K=2) = " + mtConditioned);

        double diff = refConditioned - mtConditioned;
        assertTrue("declaring an unreachable second type should strictly dilute (weaken) the " +
                        "conditioned likelihood here, not leave it unchanged or improve it",
                diff > 1e-3);
    }

    public void testTwoTypeSanity() throws Exception {

        DefaultTreeModel tree = makeTree(new String[]{"1", "2", "3"}, new String[]{"A", "B", "A"});

        int K = 2;
        String[] typeLabels = {"A", "B"};
        Parameter birthRate = new Parameter.Default(new double[]{1.2, 0.1, 0.1, 1.0}); // A->A,B ; B->B,A small cross rate
        Parameter migrationRate = new Parameter.Default(new double[]{0.0, 0.3, 0.3, 0.0});
        Parameter deathRate = new Parameter.Default(new double[]{0.3, 0.3});
        Parameter samplingRate = new Parameter.Default(new double[]{0.2, 0.2});
        Parameter treatment = new Parameter.Default(new double[]{1.0, 1.0});
        Parameter samplingProb = new Parameter.Default(new double[]{0.2, 0.2});
        Parameter originParam = new Parameter.Default(4.0);

        MultiTypeBirthDeathSerialSamplingModel mtModel = new MultiTypeBirthDeathSerialSamplingModel(
                "mtbd2", K, typeLabels, birthRate, migrationRate, deathRate, samplingRate,
                treatment, samplingProb, originParam, true, 1, Double.POSITIVE_INFINITY,
                "state", 2000, Units.Type.YEARS);

        Likelihood mtLikelihood = new SpeciationLikelihood(tree, mtModel, "mt2");
        double logL = mtLikelihood.getLogLikelihood();

        System.out.println("Two-type sanity check logL = " + logL);

        assertTrue("log-likelihood should be finite", !Double.isInfinite(logL) && !Double.isNaN(logL));
    }

    public static Test suite() {
        return new TestSuite(MultiTypeBirthDeathLikelihoodTest.class);
    }
}
