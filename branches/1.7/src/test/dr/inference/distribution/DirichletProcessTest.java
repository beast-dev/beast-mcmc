/*
 * BinomialLikelihood.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.distribution;

import dr.inference.distribution.DirichletProcessLikelihood;
import dr.inference.model.*;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * A test suit for the Dirichlet Process likelihood.
 *
 * @author Andrew Rambaut
 * @author Trevor Bedford
 * @version $Id: BinomialLikelihood.java,v 1.5 2005/05/24 20:25:59 rambaut Exp $
 */

public class DirichletProcessTest extends TestCase {

    public DirichletProcessTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();
    }

    public void testDirichletProcess() {
        // these are expected results calculated manually...
        testDirichletProcess(new double[] {1, 1, 1}, 1.0, -Math.log(6.0));  // log 1/6
        testDirichletProcess(new double[] {2, 1}, 1.0, -Math.log(6.0));  // log 1/6
        testDirichletProcess(new double[] {3}, 1.0, -Math.log(3.0));   // log 1/3

        // check the unoccupied groups are not counted
        testDirichletProcess(new double[] {2, 1, 0}, 1.0, -Math.log(6.0));  // log 1/6
        testDirichletProcess(new double[] {3, 0, 0}, 1.0, -Math.log(3.0));   // log 1/3

        // these are results calculated by the method but confirmed against an independent implementation...
        testDirichletProcess(new double[] {1, 1, 1, 1, 1}, 0.5, -6.851184927493743);
        testDirichletProcess(new double[] {2, 1, 1, 1, 0}, 0.5, -6.158037746933798);
        testDirichletProcess(new double[] {3, 1, 1, 0, 0}, 0.5, -4.771743385813907);
        testDirichletProcess(new double[] {4, 1, 0, 0, 0}, 0.5, -2.9799839165858524);
        testDirichletProcess(new double[] {5, 0, 0, 0, 0}, 0.5, -0.9005423749060166);

        testDirichletProcess(new double[] {1, 1, 1, 1, 1}, 1.0, -4.787491742782046);
        testDirichletProcess(new double[] {2, 1, 1, 1, 0}, 1.0, -4.787491742782046);
        testDirichletProcess(new double[] {3, 1, 1, 0, 0}, 1.0, -4.0943445622221);
        testDirichletProcess(new double[] {4, 1, 0, 0, 0}, 1.0, -2.995732273553991);
        testDirichletProcess(new double[] {5, 0, 0, 0, 0}, 1.0, -1.6094379124341005);

        testDirichletProcess(new double[] {1, 1, 1, 1, 1}, 2.0, -3.1135153092103747);
        testDirichletProcess(new double[] {2, 1, 1, 1, 0}, 2.0, -3.80666248977032);
        testDirichletProcess(new double[] {3, 1, 1, 0, 0}, 2.0, -3.80666248977032);
        testDirichletProcess(new double[] {4, 1, 0, 0, 0}, 2.0, -3.401197381662156);
        testDirichletProcess(new double[] {5, 0, 0, 0, 0}, 2.0, -2.7080502011022105);

        testDirichletProcess(new double[] {1, 1, 1, 1, 1}, 5.0, -1.5765840875630222);
        testDirichletProcess(new double[] {2, 1, 1, 1, 0}, 5.0, -3.1860219999971227);
        testDirichletProcess(new double[] {3, 1, 1, 0, 0}, 5.0, -4.102312731871278);
        testDirichletProcess(new double[] {4, 1, 0, 0, 0}, 5.0, -4.613138355637268);
        testDirichletProcess(new double[] {5, 0, 0, 0, 0}, 5.0, -4.836281906951478);
    }

    private void testDirichletProcess(double[] eta, double chi, double expectedLogL) {
        Parameter etaParameter = new Parameter.Default("eta", eta);
        Parameter chiParameter = new Parameter.Default("chi", chi, 0.0, Double.MAX_VALUE);

        DirichletProcessLikelihood dirichlet = new DirichletProcessLikelihood(etaParameter, chiParameter);

        assertEquals(expectedLogL, dirichlet.getLogLikelihood(), 1E-10);
    }


    public static Test suite() {
        return new TestSuite(DirichletProcessTest.class);
    }
}