/*
 * CompoundLikelihoodTest.java
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

package test.dr.inference.model;

import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import junit.framework.TestCase;

import java.util.Arrays;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

/**
 * Failures in a thread-pooled CompoundLikelihood must reach the caller, as they do when evaluated serially.
 */
public class CompoundLikelihoodTest extends TestCase {

    public void testPooledSumMatchesSerial() {
        CompoundLikelihood pooled = new CompoundLikelihood(-1, Arrays.<Likelihood>asList(
                new Fixed(-1.0), new Fixed(-2.0), new Fixed(-4.0)));
        CompoundLikelihood serial = new CompoundLikelihood(0, Arrays.<Likelihood>asList(
                new Fixed(-1.0), new Fixed(-2.0), new Fixed(-4.0)));

        assertEquals(3, pooled.getThreadCount());
        assertEquals(0, serial.getThreadCount());
        assertEquals(-7.0, pooled.getLogLikelihood(), 0.0);
        assertEquals(-7.0, serial.getLogLikelihood(), 0.0);
    }

    public void testPooledExceptionPropagates() {
        IllegalStateException failure = new IllegalStateException("partition failed");
        CompoundLikelihood pooled = new CompoundLikelihood(-1, Arrays.<Likelihood>asList(
                new Fixed(-1.0), new Throwing(failure), new Fixed(-4.0)));

        try {
            double logL = pooled.getLogLikelihood();
            fail("expected the partition's exception, got logL = " + logL);
        } catch (IllegalStateException e) {
            assertSame(failure, e);
        }
    }

    public void testPooledInterruptPropagates() {
        CompoundLikelihood pooled = new CompoundLikelihood(-1, Arrays.<Likelihood>asList(
                new Fixed(-1.0), new Blocking(new CountDownLatch(1))));

        Thread.currentThread().interrupt();
        try {
            double logL = pooled.getLogLikelihood();
            fail("expected an exception on interrupt, got logL = " + logL);
        } catch (RuntimeException e) {
            assertTrue(e.getCause() instanceof InterruptedException);
        } finally {
            assertTrue("interrupt status should be restored", Thread.interrupted());
        }
    }

    private static class Fixed extends Likelihood.Abstract {
        Fixed(double logL) {
            super(null);
            this.logL = logL;
        }

        protected double calculateLogLikelihood() {
            return logL;
        }

        private final double logL;
    }

    private static class Throwing extends Likelihood.Abstract {
        Throwing(RuntimeException failure) {
            super(null);
            this.failure = failure;
        }

        protected double calculateLogLikelihood() {
            throw failure;
        }

        private final RuntimeException failure;
    }

    private static class Blocking extends Likelihood.Abstract {
        Blocking(CountDownLatch latch) {
            super(null);
            this.latch = latch;
        }

        protected double calculateLogLikelihood() {
            try {
                latch.await(10, TimeUnit.SECONDS); // released only by the pool cancelling this task
            } catch (InterruptedException e) {
                // cancelled
            }
            return 0.0;
        }

        private final CountDownLatch latch;
    }
}
