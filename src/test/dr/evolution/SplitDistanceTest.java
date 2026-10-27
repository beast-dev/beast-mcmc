/*
 * SplitDistanceTest.java
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

package test.dr.evolution;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.treemetrics.SplitDistance;
import junit.framework.TestCase;

/**
 * Symmetric rooted and unrooted Robinson-Foulds distances (plan 7 M0).
 */
public class SplitDistanceTest extends TestCase {

    private static final String TREE = "((((a:1,b:1):1,c:2):1,d:3):2,((e:1,f:1):2,g:3):2);";

    public SplitDistanceTest(String name) {
        super(name);
    }

    // the same topology with other branch lengths and child order
    public void testSameTopology() throws Exception {
        Tree t = tree(TREE);
        Tree u = tree("((g:4,(f:2,e:2):2):1,(d:3,(c:2,(b:1,a:1):1):1):2);");
        assertEquals(0, SplitDistance.rooted(t, u));
        assertEquals(0, SplitDistance.unrooted(t, u));
    }

    // one NNI below the root changes one clade and one split: 2 both ways
    public void testOneNni() throws Exception {
        Tree t = tree(TREE);
        Tree u = tree("((((a:1,c:1):1,b:2):1,d:3):2,((e:1,f:1):2,g:3):2);");
        assertEquals(2, SplitDistance.rooted(t, u));
        assertEquals(2, SplitDistance.rooted(u, t));
        assertEquals(2, SplitDistance.unrooted(t, u));
        // n = 7: at most 2 (n - 2) = 10 rooted, 2 (n - 3) = 8 unrooted
        assertEquals(0.2, SplitDistance.normalizedRooted(t, u), 1E-12);
        assertEquals(0.25, SplitDistance.normalizedUnrooted(t, u), 1E-12);
    }

    // a root-child NNI moves the root only: the unrooted tree is the same
    public void testRerooting() throws Exception {
        Tree t = tree(TREE);
        // the root's children were (abcd) and (efg); now (abc) and (d(efg))
        Tree u = tree("(((a:1,b:1):1,c:2):3,(d:3,((e:1,f:1):1,g:2):1):2);");
        assertEquals(0, SplitDistance.unrooted(t, u));
        assertEquals(2, SplitDistance.rooted(t, u));
        // a root on a tip's edge
        Tree v = tree("(g:5,(((a:1,b:1):1,c:2):1,(d:2,(e:1,f:1):1):1):1);");
        assertTrue(SplitDistance.unrooted(t, v) > 0);
        Tree w = tree("(g:5,((((a:1,b:1):1,c:2):1,d:3):1,(e:1,f:1):3):1);");
        assertEquals(0, SplitDistance.unrooted(t, w));
    }

    // completely different resolutions reach the maximum
    public void testMaximum() throws Exception {
        Tree t = tree("(((((a:1,b:1):1,c:2):1,d:3):1,e:4):1,f:5);");
        Tree u = tree("(((((a:1,c:1):1,e:2):1,b:3):1,f:4):1,d:5);");
        assertEquals(SplitDistance.rooted(t, u) / 8.0, SplitDistance.normalizedRooted(t, u), 1E-12);
        assertTrue(SplitDistance.normalizedRooted(t, u) <= 1.0);
        Tree x = tree("((a:1,b:1):2,((c:1,d:1):1,(e:1,f:1):1):1);");
        Tree y = tree("((a:1,c:1):2,((b:1,e:1):1,(d:1,f:1):1):1);");
        // x's splits ab, cd, ef, cdef=ab; y's ac, be, df: no split in common
        assertEquals(1.0, SplitDistance.normalizedUnrooted(x, y), 1E-12);
    }

    public void testDifferentTaxaAreRefused() throws Exception {
        try {
            SplitDistance.rooted(tree(TREE), tree("((((a:1,b:1):1,c:2):1,d:3):2,((e:1,f:1):2,h:3):2);"));
            fail("different taxa");
        } catch (IllegalArgumentException expected) {
        }
    }

    private static Tree tree(String newick) throws Exception {
        return new NewickImporter(newick).importTree(null);
    }
}
