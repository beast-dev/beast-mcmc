package test.dr.evomodel.speciation;

import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;

import dr.evolution.tree.*;
import dr.evolution.util.Taxon;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesBindings.*;
import dr.evomodel.speciation.AlloppMSCoalescent;
import dr.evomodel.tree.TreeModel;


/**
 * 
 * @author Graham Jones
 *         Date: 01/06/2011
 */



/*
 * This comment describes how trees are made from XML in  the usual way.
 * It should probably be elsewhere, though it is relevant to 
 * test code which makes trees.
 * 
 * Usually, a tree is made by CoalescentSimulatorParser from XML like

   <!-- Generate a random starting tree under the coalescent process     -->
   <!-- No calibration                                              -->
   	<coalescentTree id="startingTree" rootHeight="0.018">
   		<taxa idref="taxa"/>
    		<constantSize idref="constant"/>
    	</coalescentTree>

 * which calls code in CoalescentSimulator which
 * "Simulates a set of coalescent intervals given a demographic model".
 * 
 * TreeModelParser reads XML like
   
   	<!-- Generate a tree model           -->
   		<treeModel id="treeModel">
   			<coalescentTree idref="startingTree"/>
   			<rootHeight>
   				<parameter id="treeModel.rootHeight"/>
   			</rootHeight>
   			<nodeHeights internalNodes="true">
   				<parameter id="treeModel.internalNodeHeights"/>
   			</nodeHeights>
   			<nodeHeights internalNodes="true" rootNode="true">
   				<parameter id="treeModel.allInternalNodeHeights"/>
   			</nodeHeights>
   		</treeModel>
   		
 * which uses the <coalescentTree id="startingTree"> made by 
 * CoalescentSimulator
 * 
 * The only calls to TreeModel(String, Tree) or TreeModel(Tree)
 * are from parseXMLObject() in TreeModelParser or in test code.
 * 
 */




public class AlloppSpeciesNetworkModelTEST {

	@Before
	public void setUp() throws Exception {
	}

	/*
	 * The test is too strict in that it tests against a specific newick 
	 * string when there many (128) possible correct ones.
	 */
	@Test
	public void testAlloppSpeciesNetworkModel() {
		testNetworkToMulLabTree();
		testLogLhoodGTreeInNetwork2tets5indivs();
	}
	
	
	public class NetworkToMulLabTreeTEST {
		public int testcase;
		NetworkToMulLabTreeTEST(int testcase) {
			this.testcase = testcase;
		}
	}
	
	
	public class LogLhoodGTreeInNetworkTEST {
		public TreeModel gtreemodels[];
		public double popfactors[];
		public double popvalues[];
		public double heights[];
		LogLhoodGTreeInNetworkTEST(TreeModel gtreemodels[], double popfactors[], 
				                   double popvalues[], double heights[]) {
			this.gtreemodels = gtreemodels;
			this.popfactors = popfactors;
			this.popvalues = popvalues;
			this.heights = heights;
		}
	}

	
	
	public void testNetworkToMulLabTree() {
		ApSpInfo[] apspecies = new ApSpInfo[5];
		apspecies[0] = new ApSpInfo("a", 2, new Individual[0]);
		apspecies[1] = new ApSpInfo("b", 4, new Individual[0]);
		apspecies[2] = new ApSpInfo("c", 4, new Individual[0]);
		apspecies[3] = new ApSpInfo("d", 4, new Individual[0]);
		apspecies[4] = new ApSpInfo("e", 2, new Individual[0]);

		NetworkToMulLabTreeTEST nmltTEST = new NetworkToMulLabTreeTEST(-1);
		AlloppSpeciesBindings testASB = new AlloppSpeciesBindings(apspecies, nmltTEST);
		AlloppSpeciesNetworkModel testASNM = new AlloppSpeciesNetworkModel(testASB, nmltTEST);
		String newick;
		newick = testASNM.testExampleNetworkToMulLabTree(new NetworkToMulLabTreeTEST(1));
		assertEquals(0,  newick.compareTo("((a0,((b0,c0),d0)),(e0,((b1,c1),d1)))"));
		newick = testASNM.testExampleNetworkToMulLabTree(new NetworkToMulLabTreeTEST(2));
		assertEquals(0,  newick.compareTo("(((a0,((b0,c0),d0)),((b1,c1),d1)),e0)"));
		newick = testASNM.testExampleNetworkToMulLabTree(new NetworkToMulLabTreeTEST(3));
		assertEquals(0,  newick.compareTo("((a0,(((b0,c0),d0),((b1,c1),d1))),e0)"));
		newick = testASNM.testExampleNetworkToMulLabTree(new NetworkToMulLabTreeTEST(4));
		assertEquals(0,  newick.compareTo("(((a0,(b0,c0)),(d0,d1)),(e0,(b1,c1)))"));
		newick = testASNM.testExampleNetworkToMulLabTree(new NetworkToMulLabTreeTEST(5));
		assertEquals(0,  newick.compareTo("(((((a0,b0),c0),c1),(d0,d1)),(e0,b1))"));
		}
	
	
	
	public void testLogLhoodGTreeInNetwork2tets5indivs() {
		Taxon a1taxa[] = new Taxon[2];
		Taxon a2taxa[] = new Taxon[2];
		Taxon b1taxa[] = new Taxon[2];
		Taxon b2taxa[] = new Taxon[2];
		Taxon b3taxa[] = new Taxon[2];
		AlloppSpeciesBindings.Individual aIndivs[] = new AlloppSpeciesBindings.Individual[2];
		AlloppSpeciesBindings.Individual bIndivs[] = new AlloppSpeciesBindings.Individual[3];
		
		ApSpInfo[] apsp = new ApSpInfo[2];
		a1taxa[0] = new Taxon("a1A");
		a1taxa[1] = new Taxon("a1B");
		aIndivs[0] = new AlloppSpeciesBindings.Individual("a1", a1taxa);
		a2taxa[0] = new Taxon("a2A");
		a2taxa[1] = new Taxon("a2B");
		aIndivs[1] = new AlloppSpeciesBindings.Individual("a2", a2taxa);
		apsp[0] = new ApSpInfo("a", 4, aIndivs);
		
		b1taxa[0] = new Taxon("b1A");
		b1taxa[1] = new Taxon("b1B");
		bIndivs[0] = new AlloppSpeciesBindings.Individual("b1", b1taxa);
		b2taxa[0] = new Taxon("b2A");
		b2taxa[1] = new Taxon("b2B");
		bIndivs[1] = new AlloppSpeciesBindings.Individual("b2", b2taxa);	
		b3taxa[0] = new Taxon("b3A");
		b3taxa[1] = new Taxon("b3B");
		bIndivs[2] = new AlloppSpeciesBindings.Individual("b3", b3taxa);	
		apsp[1] = new ApSpInfo("b", 4, bIndivs);	
        SimpleNode[] nodes = new SimpleNode[19];
        for (int n=0; n < 19; n++) {
            nodes[n] = new SimpleNode();
        }
        		
/*
This is R code which produces 17.87787855
last two cols are gene tree node times, and species tree node and hybridization times 


#	000 010 100 110 120 101 111 121 001 011    sp,indiv,seq indices
#	a1A a2A b1A b2A b3A b1B b2B b3B a1B a2B    names
#	 0   1   2   5   6   8   9  11  13  14
#	 |	 | # |    \ /  #  \ /   | #  \ /
#	 |   | # |     7   #   10  /  #   15     0.01
#	 |   | # |     |   #    \ /   #  /
#    \   | # |     |   #    12   #  /        0.02
#	   \  \#/      |   #    |    # /                  0.03 tetra split
#       \  3       |   #    |     /          0.04
#	     \ |       |   #    |    /
#	       4       |    #   |   /            0.05
#	        \       \   +   |  /                      0.06 hyb
#	          \      |  #   | /
#	            \    |  #   16               0.08
#	             \    \ #  /                          0.09 diploid split
#	              \    |  /
#                  \   17                    0.11
#	                \  /
#	                 18                      0.12


# 	 0   1   2   5   6   8   9  11  13  14
#	 .01 2   #   .03 3   #  .03 3   #  .01 2  
#		       #           #          #          0.01
#		br1    #    br2    #    br6   #   br5
#	         #           #         #           0.02
#	 .02 2   #  .04 2    # .04 1   #  .02 1             0.03 tetra split
#        .06 4         #      .06  2         0.04
#		                    #          
#		       br3          #       br7          0.05
#		     .05,.06 2      +    .05,.06 2                0.06 hyb
#		       br4          #       br8
#		                    #                    0.08
#		        .06  2      #   .06 1                     0.09 diploid split
#		                 .12  3 
#                                            0.11
#		                   br9
#		                 .12  1                  0.12



brs <- matrix(0, nrow=9, ncol=3)
tms <- matrix(0, nrow=9, ncol=4)
brs[1,] <- c(.01,  .02,    2);     tms[1,] <- c(.00,              .03, -1,-1)
brs[2,] <- c(.03,  .04,    3);     tms[2,] <- c(.00,  .01,        .03, -1)
brs[3,] <- c(.06,  .05,    4);     tms[3,] <- c(.03,  .04, .05,   .06)
brs[4,] <- c(.06,  .06,    2);     tms[4,] <- c(.06,              .09, -1,-1)

brs[5,] <- c(.01,  .02,    2);     tms[5,] <- c(.00,  .01,        .03, -1)
brs[6,] <- c(.03,  .04,    3);     tms[6,] <- c(.00,  .01, .02,   .03)
brs[7,] <- c(.06,  .05,    2);     tms[7,] <- c(.03,              .06, -1,-1)
brs[8,] <- c(.06,  .06,    2);     tms[8,] <- c(.06,  .08,        .09, -1)

brs[9,] <- c(.12,  .12,    3);     tms[9,] <- c(.09,  .11, .12,   .18)



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
        setNode(nodes[0],  apsp[0].taxonFromIndSeq(0, 0), 0.00, null,      null);
        setNode(nodes[1],  apsp[0].taxonFromIndSeq(1, 0), 0.00, null,      null);
        setNode(nodes[2],  apsp[1].taxonFromIndSeq(0, 0), 0.00, null,      null);
        setNode(nodes[3],  new Taxon(""),                 0.04, nodes[1],  nodes[2]);
        setNode(nodes[4],  new Taxon(""),                 0.05, nodes[0],  nodes[3]);
        setNode(nodes[5],  apsp[1].taxonFromIndSeq(1, 0), 0.00, null,      null);
        setNode(nodes[6],  apsp[1].taxonFromIndSeq(2, 0), 0.00, null,      null);
        setNode(nodes[7],  new Taxon(""),                 0.01, nodes[5],  nodes[6]);
        setNode(nodes[8],  apsp[1].taxonFromIndSeq(0, 1), 0.00, null,      null);
        setNode(nodes[9],  apsp[1].taxonFromIndSeq(1, 1), 0.00, null,      null);
        setNode(nodes[10], new Taxon(""),                 0.01, nodes[8],  nodes[9]);
        setNode(nodes[11], apsp[1].taxonFromIndSeq(2, 1), 0.00, null,      null);
        setNode(nodes[12], new Taxon(""),                 0.02, nodes[10], nodes[11]);
        setNode(nodes[13], apsp[0].taxonFromIndSeq(0, 1), 0.00, null,      null);
        setNode(nodes[14], apsp[0].taxonFromIndSeq(1, 1), 0.00, null,      null);
        setNode(nodes[15], new Taxon(""),                 0.01, nodes[13], nodes[14]);
        setNode(nodes[16], new Taxon(""),                 0.08, nodes[12], nodes[15]);
        setNode(nodes[17], new Taxon(""),                 0.11, nodes[7],  nodes[16]);
        setNode(nodes[18], new Taxon(""),                 0.12, nodes[4],  nodes[17]);
        SimpleTree gtree = new SimpleTree(nodes[18]);
        TreeModel gtreemodels[] = new TreeModel[1];
        gtreemodels[0] = new TreeModel(gtree);
        double popfactors[] = new double[1];
        popfactors[0] = 2.0;
 		double popvalues[] = new double[6];
		for (int i = 0; i < popvalues.length; i++) { 
			popvalues[i] = 0.01 * (i+1);
		}
		double heights[] = new double[3];
		heights[0] = 0.03; heights[1] = 0.06; heights[2] = 0.09;       
        LogLhoodGTreeInNetworkTEST llgtnTEST = new 
                   LogLhoodGTreeInNetworkTEST(gtreemodels, popfactors, popvalues, heights);
		AlloppSpeciesBindings testASB = new AlloppSpeciesBindings(apsp, gtreemodels, llgtnTEST);
		AlloppSpeciesNetworkModel testASNM = new AlloppSpeciesNetworkModel(testASB, llgtnTEST);
		AlloppMSCoalescent testASMSC = new AlloppMSCoalescent(testASB, testASNM, llgtnTEST);
		double llhd = testASMSC.getLogLikelihood();
		String mullabtext = testASNM.mullabTreeAsText();
		System.out.println(mullabtext);
		assertEquals(17.87787855, llhd, 1e-6);
	}
	
	
	
	
	
	
	
	private void setNode(SimpleNode node, Taxon taxon, double height, SimpleNode lft, SimpleNode rgt) {
		node.setTaxon(taxon);
		node.setHeight(height);
		if (lft != null) { node.addChild(lft); }
		if (rgt != null) { node.addChild(rgt); }
	}
}
