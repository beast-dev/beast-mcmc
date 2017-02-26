/*
 * AlloppSpeciesNetworkModelTEST.java
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

package test.dr.evomodel.alloppnet.speciation;

import dr.evomodel.alloppnet.speciation.AlloppMulLabTree;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings.ApSpInfo;
import dr.math.MathUtils;
import dr.math.distributions.GammaDistribution;
import org.junit.Before;
import org.junit.Test;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesBindings.Individual;
import dr.evomodel.alloppnet.speciation.AlloppSpeciesNetworkModel;



import static org.junit.Assert.assertEquals;


/**
 *
 * @author Graham Jones
 *         Date: 01/06/2011
 */


public class AlloppSpeciesNetworkModelTEST {

    @Before
    public void setUp() throws Exception {
    }




    @Test
    public void testAlloppSpeciesNetworkModel() {

        //testGammaQuantile();
        testNetworkToMulLabTree();
        testLhoodMulLabTree();
        // grjtodo-soon would like a test of gene tree in network compatibility.
    }


	/*
	 *  ****************** TESTING MulLabTree conversion ******************
	 *
	 */


    // trying to reproduce weird bug when quantile returned 4.9e-324
    public void testGammaQuantile() {
        MathUtils.setSeed(42);
        for (int i = 0; i < 1000000; i++) {
            double s = MathUtils.uniform(.999, 1.001);
            double b = MathUtils.uniform(6e-5, 7e-5);
            final GammaDistribution gamma = new GammaDistribution(s,b);
            double q = MathUtils.uniform(3e-7, 4e-7);
            double x = gamma.quantile(q);
            assert x > 1e-20;
        }
    }


    public void testNetworkToMulLabTree() {
        ApSpInfo[] apspecies = new ApSpInfo[5];
        apspecies[0] = new ApSpInfo("a", 2, new Individual[0]);
        apspecies[1] = new ApSpInfo("b", 4, new Individual[0]);
        apspecies[2] = new ApSpInfo("c", 4, new Individual[0]);
        apspecies[3] = new ApSpInfo("d", 4, new Individual[0]);
        apspecies[4] = new ApSpInfo("e", 2, new Individual[0]);

        AlloppSpeciesBindings testASB = new AlloppSpeciesBindings(apspecies);
        AlloppSpeciesNetworkModel testASNM = new AlloppSpeciesNetworkModel(testASB);

        System.out.println("Tests of Network To MulLabTree conversion with dips a,e and tets b,c,d");
        String newick;
        newick = testASNM.testExampleNetworkToMulLabTree(1);
        System.out.println(newick + "\n\n");
        assertEquals(0,  newick.compareTo("((((b0,c0),d0),a),(((b1,c1),d1),e))"));
        newick = testASNM.testExampleNetworkToMulLabTree(2);
        System.out.println(newick + "\n\n");
        assertEquals(0,  newick.compareTo("(((((b0,c0),d0),a),((b1,c1),d1)),e)"));
        newick = testASNM.testExampleNetworkToMulLabTree(3);
        System.out.println(newick + "\n\n");
        assertEquals(0,  newick.compareTo("(((((b0,c0),d0),((b1,c1),d1)),a),e)"));
        newick = testASNM.testExampleNetworkToMulLabTree(4);
        System.out.println(newick + "\n\n");
        assertEquals(0,  newick.compareTo("((((b0,c0),a),(d0,d1)),((b1,c1),e))"));
        newick = testASNM.testExampleNetworkToMulLabTree(5);
        System.out.println(newick + "\n\n");
        assertEquals(0,  newick.compareTo("(((((a,b0),c0),c1),(d0,d1)),(b1,e))"));
        System.out.println("");
        System.out.println("");
    }


    public void testLhoodMulLabTree() {
        ApSpInfo[] apspecies = new ApSpInfo[3];
        apspecies[0] = new ApSpInfo("a", 2, new Individual[0]);
        apspecies[1] = new ApSpInfo("b", 2, new Individual[0]);
        apspecies[2] = new ApSpInfo("z", 4, new Individual[0]);

        AlloppSpeciesBindings testASB = new AlloppSpeciesBindings(apspecies);
        AlloppMulLabTree testamlt = new AlloppMulLabTree(testASB);
        double  llhood = testamlt.testGeneTreeInMULTreeLogLikelihood();
        assertEquals(llhood, -13.562218135041552713, 1e-10);
        System.out.println(llhood);
    }


    /*   This R code produces -13.562218135041552713

#     0   2   3   1
#     a   z0  z1  b b'
#     | # | # | # | |
#     | # | + | # | |  hyb 0.005
#     | 4 | # | # | |  (a,z0) in multree 0.01
#      \ /  # | # | |
#       v   # | # | |   (a,z0) in gene tree  0.015
#       |   # | 5 | |   (z1,b) in multree    0.02
#       |   #  \ / /
#       |   #   v /      (z1,b) in gene tree    0.025
#       |   #   | |
#        \  #  / /
#         | 6  | |      root multree 0.03
#         |    |/       (z1,b),b'    0.035
#         |   /
#          \ /
#           v         root gene tree 0.045
#           |


# here, numbers 1-9 show branch or part-branch indices
#     a   z0  z1  b b'
#     | # 2 # 3 # | |
#     1 # | + | # | |  hyb 0.005
#     | # 4 # 5 # |6|  (a,z0) in multree 0.01
#      \ /  # | # | |
#       v   # | # | |   (a,z0) in gene tree  0.015
#       |   # | # | |   (z1,b) in multree    0.02
#       |   #  \ / /
#      7|   #   v /      (z1,b) in gene tree    0.025
#       |   #  8| |
#        \  #  / /
#         | #  | |      root multree 0.03
#         |    |/       (z1,b),b'    0.035
#         | 9  /
#          \ /
#           v         root gene tree 0.045
#           |




brs <- matrix(0, nrow=9, ncol=3)
tms <- matrix(0, nrow=9, ncol=4)
brs[1,] <- c(.003,  .001,    1);     tms[1,] <- c(.0,         .01,  -1,-1)
brs[2,] <- c(.003,  .001,    1);     tms[2,] <- c(.00,        .005, -1,-1)
brs[3,] <- c(.003,  .001,    1);     tms[3,] <- c(.00,        .005, -1,-1)
brs[4,] <- c(.001,  .001,    1);     tms[4,] <- c(.005,       .01,  -1,-1)
brs[5,] <- c(.001,  .001,    1);     tms[5,] <- c(.005,       .01,  -1,-1)
brs[6,] <- c(.003,  .001,    2);     tms[6,] <- c(.0,         .02,  -1,-1)
brs[7,] <- c(.002,  .001,    2);     tms[7,] <- c(.01,  .015, .03,  -1)
brs[8,] <- c(.002,  .001,    3);     tms[8,] <- c(.02,  .025, .03,  -1)
brs[9,] <- c(.002,  .002,    3);     tms[9,] <- c(.03,  .035, .045, 999)



branchllhood <- function(b, tm) {
  llhood <- 0
  tm <- tm[tm>=0]
  tm <- tm - tm[1]
  k <- length(tm)
  n <- b[3]
  p0 <- b[1]
  pk <- b[2]
  for (i in 1:(k-1)) {
    pi0 <- linear.interp.pop(tm[1], tm[k], p0, pk, tm[i])
    pi1 <- linear.interp.pop(tm[1], tm[k], p0, pk, tm[i+1])
    if (i < k-1) {
      x <- onecoalllhood(n, tm[i], tm[i+1], pi0, pi1)
      n <- n-1
      } else {
      x <- nocoalllhood(n, tm[i], tm[i+1], pi0, pi1)
      }
    llhood <- llhood + x
    }
  llhood
  }

linear.interp.pop <- function(t0, t1, p0, p1, x)
{
(p0*(t1-x) + p1*(x-t0))/(t1-t0)
}

onecoalllhood <- function(n, t0, t1, p0, p1) {
  x <- -log(p1) - integrate(pop.integrand, n=n, t0=t0, t1=t1, p0=p0, p1=p1, lower=t0, upper=t1)$value
  cat("onecoalllhood", "n", n, "t0", t0, "t1", t1, "p0", p0, "p1", p1, "llhood", x, "\n")
  x
  }


nocoalllhood <- function(n, t0, t1, p0, p1) {
  x <- - integrate(pop.integrand, n=n, t0=t0, t1=t1, p0=p0, p1=p1, lower=t0, upper=t1)$value
  cat("nocoalllhood", "n", n, "t0", t0, "t1", t1, "p0", p0, "p1", p1, "llhood", x, "\n")
  x
  }


pop.integrand <- function(x, n, t0, t1, p0, p1)
{
(n*(n-1)/2) * (t1-t0) / (p0*(t1-x) + p1*(x-t0))
}


llhood <- 0
for (i in 1:dim(brs)[1]) {
  x <- branchllhood(brs[i,], tms[i,])
  cat("branch", i, "llhood", x, "\n")
  llhood <- llhood + x
  }
llhood

     */

}

