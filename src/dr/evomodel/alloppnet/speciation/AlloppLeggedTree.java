/*
 * AlloppLeggedTree.java
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

package dr.evomodel.alloppnet.speciation;

import java.util.ArrayList;
import java.util.Formatter;
import java.util.Locale;
import java.util.Stack;


import dr.evomodel.alloppnet.util.AlloppMisc;
import jebl.util.FixedBitSet;

import dr.evolution.tree.NodeRef;
import dr.evomodel.alloppnet.tree.SlidableTree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;


/**
 *
 * A 'tree with legs' for a single ploidy level in an allopolyploid network.
 *
 * @author Graham Jones
 *         Date: 01/05/2011
 */



/*
 * class AlloppLeggedTree
 *
 * This is a `tree with legs', which is a homoploid
 * species tree which is attached to a tree of lower ploidy
 * via its legs, as part of a AlloppSpeciesNetworkModel.
 *
 * altnodes is an array of ALTNode's which implements the homoploid species tree.
 * rootn is the index of the root node.
 *
 * The fields diphistlftleg, diphistrgtleg are indices giving the
 * hyb-tips in the AlloppDiploidHistory.
 *
 */




public class AlloppLeggedTree implements  SlidableTree  {

    private ALTNode[] altnodes;
    private int rootn;

    private int diphistlftleg;
    private int diphistrgtleg;



    /*
    * parent, child[] implement the tree topology.
    *
    * height is the node height; can be > 0 for hyb-tips.
    *
    * union is a spseq-union. At diploid tips they are used normally.
    * Hyb-tips have unions derived from the tetratree and leg index (0 or 1).
    * Unions are taken towards the root.
    *
    * For a hyb-tip, tettree specifies the index of the tetratree whose
    * root comes from it. It is not used for other nodes.
    */
    private class ALTNode extends AlloppNode.Abstract implements AlloppNode, NodeRef {
        private int anc;
        private int lft;
        private int rgt;
        private double height;
        private Taxon taxon;
        private FixedBitSet union;
        private int nodeNumber;

        // dud constuctor
        ALTNode(int nn) {
            anc = -1;
            lft = -1;
            rgt = -1;
            height = -1.0;
            taxon = new Taxon("");
            union = null;
            nodeNumber = nn;
        }

        // copy constructor
        public ALTNode(ALTNode node) {
            anc = node.anc;
            lft = node.lft;
            rgt = node.rgt;
            nodeNumber = node.nodeNumber;
            copyNonTopologyFields(node);
        }

        // no-topology constructor. Copies all fields of argument node,
        // except the anc, lft, rgt fields which are set to `unknown',
        // and the index field, which is set to nn
        public ALTNode(int nn, ALTNode node) {
            anc = -1;
            lft = -1;
            rgt = -1;
            nodeNumber = nn;
            copyNonTopologyFields(node);
        }



        private void copyNonTopologyFields(ALTNode node) {
            height = node.height;
            taxon = new Taxon(node.taxon.getId());
            if (node.union == null) {
                union = null;
            } else {
                union = new FixedBitSet(node.union);
            }
        }

        @Override
        public AlloppNode getChild(int ch) {
            return ch==0 ? altnodes[lft] : altnodes[rgt];
        }

        @Override
        public AlloppNode getAnc() {
            return altnodes[anc];
        }


        @Override
        public double getHeight() {
            return height;
        }

        @Override
        public Taxon getTaxon() {
            return taxon;
        }

        @Override
        public FixedBitSet getUnion() {
            return union;
        }

        @Override
        public void setChild(int ch, AlloppNode newchild) {
            int newch = ((ALTNode)newchild).nodeNumber;
            if (ch == 0) {
                lft = newch;
            } else {
                rgt = newch;
            }
        }

        @Override
        public void setAnc(AlloppNode anc) {
            this.anc = ((ALTNode)anc).nodeNumber;
        }

        @Override
        public void setTaxon(String name) {
            this.taxon = new Taxon(name);

        }


        @Override
        public void setHeight(double height) {
            this.height = height;

        }

        @Override
        public void setUnion(FixedBitSet union) {
            this.union = union;
        }

        @Override
        public void addChildren(AlloppNode c0, AlloppNode c1) {
            lft = ((ALTNode)c0).nodeNumber;
            altnodes[lft].anc = nodeNumber;
            rgt = ((ALTNode)c1).nodeNumber;
            altnodes[rgt].anc = nodeNumber;
        }

        @Override
        public String asText(int indentlen) {
            StringBuilder s = new StringBuilder();
            Formatter formatter = new Formatter(s, Locale.US);
            if (lft < 0) {
                formatter.format("%s ", taxon.getId());
            } else {
                formatter.format("%s ", "+");
            }
            while (s.length() < 20-indentlen) {
                formatter.format("%s", " ");
            }
            formatter.format("%s ", AlloppMisc.nonnegIn8Chars(height));
            return s.toString();
        }


        @Override
        public int nofChildren() {
            return (lft < 0) ? 0 : 2;
        }

        @Override
        public int getNumber() {
            return nodeNumber;
        }

        @Override
        public void setNumber(int n) {
            nodeNumber = n;

        }



    }


    /*
     * Constructor makes a random starting (homoploid) tree.
     * Legs are left to AlloppDiploidHistory to do.
     */
    public AlloppLeggedTree(Taxon[] taxa, double rate) {
        int noftets = taxa.length;

        // Make array of dud nodes
        altnodes = new ALTNode[2 * noftets - 1];
        for (int i = 0; i < altnodes.length; i++) {
            altnodes[i] = new ALTNode(i);
        }
        ArrayList<Integer> tojoin = new ArrayList<Integer>(noftets);
        for (int n = 0; n < noftets; n++) {
            altnodes[n].setTaxon(taxa[n].getId());
            altnodes[n].setHeight(0.0);
            tojoin.add(n);
        }
        double treeheight = 0.0;
        for (int i = 0; i < noftets-1; i++) {
            int numtojoin = tojoin.size();
            int j = MathUtils.nextInt(numtojoin);
            Integer child0 = tojoin.get(j);
            tojoin.remove(j);
            int k = MathUtils.nextInt(numtojoin-1);
            Integer child1 = tojoin.get(k);
            tojoin.remove(k);
            altnodes[noftets+i].addChildren(altnodes[child0], altnodes[child1]);
            altnodes[noftets+i].setHeight(treeheight + randomnodeheight(numtojoin*rate));
            treeheight = altnodes[noftets+i].getHeight();
            tojoin.add(noftets+i);
        }
        diphistlftleg = -1;
        diphistrgtleg = -1;
        rootn = altnodes.length - 1;
    }



    /**
     * copy constructor
     */
    public AlloppLeggedTree(AlloppLeggedTree x) {
        altnodes = new ALTNode[x.altnodes.length];
        for (int n = 0; n < altnodes.length; n++) {
            altnodes[n] = new ALTNode(x.altnodes[n]);
        }
        rootn = x.rootn;
        this.diphistlftleg = x.diphistlftleg;
        this.diphistrgtleg = x.diphistrgtleg;

    }


    /*
      * Constructor. Makes a merged tree from two trees. ttree2 has the more ancient hyb time,
      * so the merged tree gets its legs.
      */
    public AlloppLeggedTree(AlloppLeggedTree ttree1, AlloppLeggedTree ttree2, double hybHeight) {
        altnodes = new ALTNode[1 + ttree1.altnodes.length + ttree2.altnodes.length];
        diphistlftleg = ttree2.diphistlftleg;
        diphistrgtleg = ttree2.diphistrgtleg;
        int nextn = copySubtree(0, (ALTNode)ttree1.getSlidableRoot());
        int lft = nextn - 1;
        nextn = copySubtree(nextn, (ALTNode)ttree2.getSlidableRoot());
        int rgt = nextn - 1;
        assert nextn == altnodes.length - 1;
        altnodes[nextn] = new ALTNode(nextn);
        altnodes[nextn].addChildren(altnodes[lft], altnodes[rgt]);
        altnodes[nextn].setHeight(hybHeight);
        rootn = nextn;
    }


    /*
    * Constructor. Makes a tree from  subtree of tetTree
    * does not fill in legs
    */
    public AlloppLeggedTree(AlloppLeggedTree tetTree, AlloppNode sub) {
        ALTNode node = (ALTNode)sub;
        int ntips = tetTree.noftipsSubtree(node);
        altnodes = new ALTNode[2*ntips-1];
        for (int n = 0; n < altnodes.length; n++) {
            altnodes[n] = new ALTNode(n);
        }
        int nextn = copySubtree(0, node);
        rootn = nextn - 1;
    }


    /*
    * Constructor for testing.
    */
    public AlloppLeggedTree(Taxon[] taxa) {
        int nTaxa = taxa.length;
        assert(nTaxa <= 4);
        int nNodes = 2 * nTaxa - 1;
        altnodes = new ALTNode[nNodes];
        for (int n = 0; n < nNodes; n++) {
            altnodes[n] = new ALTNode(n);
        }

        for (int t = 0; t<nTaxa; t++) {
            altnodes[t].setTaxon(taxa[t].getId());
            altnodes[t].setHeight(0.0);
        }
        if (nTaxa == 2) {
            altnodes[2].setHeight(1.0);
            altnodes[2].addChildren(altnodes[0], altnodes[1]);
        }
        if (nTaxa == 3) {
            altnodes[3].setHeight(altnodes[0].getHeight() + 1.0);
            altnodes[3].addChildren(altnodes[0], altnodes[1]);
            altnodes[4].setHeight(altnodes[3].getHeight() + 1.0);
            altnodes[4].addChildren(altnodes[2], altnodes[3]);
        }
        if (nTaxa == 4) {
            altnodes[4].setHeight(altnodes[0].getHeight() + 1.0);
            altnodes[4].addChildren(altnodes[0], altnodes[1]);
            altnodes[5].setHeight(altnodes[4].getHeight() + 1.0);
            altnodes[5].addChildren(altnodes[2], altnodes[4]);
            altnodes[6].setHeight(altnodes[5].getHeight() + 1.0);
            altnodes[6].addChildren(altnodes[3], altnodes[5]);
        }
        rootn = altnodes.length - 1;
    }



    // SlidableTree implementation

    @Override
    public NodeRef getSlidableRoot() {
        assert altnodes[rootn].anc < 0;
        return altnodes[rootn];
    }



    @Override
    public void replaceSlidableRoot(NodeRef root) {
        rootn = root.getNumber();
        altnodes[rootn].anc = -1;
    }



    @Override
    public int getSlidableNodeCount() {
        return altnodes.length;
    }



    @Override
    public double getSlidableNodeHeight(NodeRef node) {
        return altnodes[node.getNumber()].getHeight();
    }

    @Override
    public Taxon getSlidableNodeTaxon(NodeRef node) {
        return altnodes[node.getNumber()].getTaxon();
    }

    @Override
    public void setSlidableNodeHeight(NodeRef node, double height) {
        altnodes[node.getNumber()].height = height;

    }



    @Override
    public boolean isExternalSlidable(NodeRef node) {
        return (altnodes[node.getNumber()].lft < 0);
    }



    @Override
    public NodeRef getSlidableChild(NodeRef node, int j) {
        int n = node.getNumber();
        return j == 0 ? altnodes[ altnodes[n].lft ] : altnodes[ altnodes[n].rgt ];
    }



    @Override
    public void replaceSlidableChildren(NodeRef node, NodeRef lft, NodeRef rgt) {
        int nn = node.getNumber();
        int lftn = lft.getNumber();
        int rgtn = rgt.getNumber();
        assert altnodes[nn].lft >= 0;
        altnodes[nn].lft = lftn;
        altnodes[nn].rgt = rgtn;
        altnodes[lftn].anc = altnodes[nn].nodeNumber;
        altnodes[rgtn].anc = altnodes[nn].nodeNumber;
    }



    String asText(int tt) {
        String header = "Tetraploid tree " + String.valueOf(tt) + "     height" + System.getProperty("line.separator");
        String s = "";
        Stack<Integer> x = new Stack<Integer>();
        return header + AlloppNode.Abstract.subtreeAsText(altnodes[rootn], s, x, 0, "");

    }




    boolean leggedtreeOK() {
        int nroots = 0;
        for (int i = 0; i < altnodes.length; i++) {
            if (altnodes[i].anc < 0) {
                nroots++;
            }
        }
        if (nroots != 1) {
            return false;
        }
        for (int i = 0; i < altnodes.length; i++) {
            int nparents = 0;
            for (int j = 0; j < altnodes.length; j++) {
                if (altnodes[j].lft == i) { nparents++; }
                if (altnodes[j].rgt == i) { nparents++; }
            }
            if (altnodes[i].anc < 0  &&  nparents != 0) {
                return false;
            }
            if (altnodes[i].anc >= 0  &&  nparents != 1) {
                return false;
            }
        }
        for (int i = 0; i < altnodes.length; i++) {
            if (altnodes[i].getNumber() != i) {
                return false;
            }
        }
        for (int i = 0; i < altnodes.length; i++) {
            if (altnodes[i].lft >= 0) {
                if (altnodes[i].rgt < 0) {
                    return false;
                }
                int lft = altnodes[i].lft;
                int rgt = altnodes[i].rgt;
                if (altnodes[lft].anc != i) {
                    return false;
                }
                if (altnodes[rgt].anc != i) {
                    return false;
                }
                if (altnodes[i].height <= altnodes[lft].height) {
                    return false;
                }
                if (altnodes[i].height <= altnodes[rgt].height) {
                    return false;
                }
            } else {
                if (altnodes[i].height != 0) {
                    return false;
                }
            }
        }
        if (altnodes[rootn].anc >= 0) {
            return false;
        }
        return true;
    }





    public int scaleAllHeights(double scale) {
        int count = 0;
        for (int n = 0; n < altnodes.length; n++) {
            if (altnodes[n].nofChildren() > 0) {
                altnodes[n].height *= scale;
                count++;
            }
        }
        return count;
    }


    public ArrayList<Taxon> getSpeciesTaxons() {
        ArrayList<Taxon> sptxs = new ArrayList<Taxon>();
        for (int n = 0; n < altnodes.length; n++) {
            if (altnodes[n].nofChildren() == 0) {
                Taxon taxon = altnodes[n].getTaxon();
                sptxs.add(taxon);
            }
        }
        assert sptxs.size() == getExternalNodeCount();
        return sptxs;
    }




    public void fillinTipUnions(AlloppSpeciesBindings apsp, int leg) {
        for (int n = 0; n < altnodes.length; n++) {
            if (altnodes[n].nofChildren() == 0) {
                altnodes[n].setUnion(apsp.taxonseqToTipUnion(altnodes[n].taxon, leg));
            }
        }
    }

    public double getRootHeight() {
        return altnodes[rootn].height;
    }



    public int getExternalNodeCount() {
        return (altnodes.length + 1) / 2;
    }

    public int getInternalNodeCount() {
        return (altnodes.length - 1) / 2;
    }




    public void collectInternalHeights(ArrayList<Double> heights) {
        for (int n = 0; n < altnodes.length; n++) {
            if (altnodes[n].nofChildren() > 0) {
                heights.add(altnodes[n].height);
            }
        }
    }


    public void setDiphistLftLeg(int lftleg) {
        diphistlftleg = lftleg;
    }

    public void setDiphistRgtLeg(int rgtleg) {
        diphistrgtleg = rgtleg;
    }

    public int getDiphistLftLeg() {
        return diphistlftleg;
    }

    public int getDiphistRgtLeg() {
        return diphistrgtleg;
    }




    /***********************************************************************/
    /************************** private ************************************/
    /***********************************************************************/

    private int copySubtree(int nextn, ALTNode node) {
        if (node.nofChildren() == 0) {
            altnodes[nextn] = new ALTNode(nextn, node);
            nextn++;
        } else {
            nextn = copySubtree(nextn, (ALTNode)node.getChild(0));
            int lft = nextn - 1;
            nextn = copySubtree(nextn, (ALTNode)node.getChild(1));
            int rgt = nextn - 1;
            altnodes[nextn] = new ALTNode(nextn, node);
            altnodes[nextn].anc = -1;
            altnodes[nextn].nodeNumber = nextn;
            altnodes[nextn].addChildren(altnodes[lft], altnodes[rgt]);
            nextn++;
        }
        return nextn;
    }

    private int noftipsSubtree(ALTNode node) {
        int ntips = 0;
        if (node.lft >= 0) {
            int lftntips = noftipsSubtree(altnodes[node.lft]);
            int rgtntips = noftipsSubtree(altnodes[node.rgt]);
            ntips = lftntips + rgtntips;
        } else {
            ntips = 1;
        }
        return ntips;
    }


    private double randomnodeheight(double rate) {
        return MathUtils.nextExponential(rate) + 1e-6/rate;
        // 1e-6/rate to avoid very tiny heights
    }










}
