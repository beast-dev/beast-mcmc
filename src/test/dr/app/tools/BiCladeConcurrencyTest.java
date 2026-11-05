/*
 * BiCladeConcurrencyTest.java
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

package test.dr.app.tools;

import dr.app.tools.treeannotator.BiClade;
import junit.framework.TestCase;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;

/**
 * TreeAnnotator collects clade heights on a thread pool: concurrent addChildHeightValues on one clade must keep every
 * value, with the left and right values paired.
 *
 * @author Xiang Ji
 */
public class BiCladeConcurrencyTest extends TestCase {

    public void testConcurrentChildHeightValues() throws Exception {
        final int threads = 8;
        final int perThread = 20000;
        for (int round = 0; round < 5; ++round) {
            final BiClade clade = new BiClade("clade", 2);
            ExecutorService pool = Executors.newFixedThreadPool(threads);
            List<Future<?>> futures = new ArrayList<>();
            for (int t = 0; t < threads; ++t) {
                final double left = t;
                futures.add(pool.submit(() -> {
                    for (int i = 0; i < perThread; ++i) {
                        clade.addChildHeightValues(left, left + 0.5);
                    }
                }));
            }
            for (Future<?> future : futures) {
                future.get();
            }
            pool.shutdown();
            assertEquals(threads * perThread, clade.getLeftHeightValues().size());
            assertEquals(threads * perThread, clade.getRightHeightValues().size());
            for (int i = 0; i < clade.getLeftHeightValues().size(); ++i) {
                assertEquals(clade.getLeftHeightValues().get(i) + 0.5, clade.getRightHeightValues().get(i), 0.0);
            }
        }
    }
}
