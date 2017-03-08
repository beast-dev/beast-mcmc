/*
 * MulMSCoalescent.java
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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ScaledDemographic;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Units;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import jebl.util.FixedBitSet;

import java.util.Arrays;

/**
 * Compute coalescent log-likelihood of a set of gene trees embedded inside one multiply labelled species tree.
 *
 * @author Joseph Heled, Graham Jones
 *         Date: 29/12/2011
 */


/*
 * Very similar to MultiSpeciesCoalescent. (10/01/2012)
 * 
 */

public class MulMSCoalescent extends Likelihood.Abstract implements Units {
    private final MulSpeciesTreeModel spTree;
    private final MulSpeciesBindings mulspb;
    private boolean checkCompatibility;  //grj not using this yet
    private final boolean[] compatibleCheckRequired;  //grj not using this yet

    public MulMSCoalescent(MulSpeciesBindings mulspb, MulSpeciesTreeModel sptree) {
        super(sptree);
        spTree = sptree;
        this.mulspb = mulspb;

        spTree.addModelRestoreListener(this);

        // recompute on any change in geneTree -
        // possible optimization: keep track which tree changed.
        final MulSpeciesBindings.GeneTreeInfo[] trees = mulspb.getGeneTrees();
        for(MulSpeciesBindings.GeneTreeInfo geneTree : trees) {
            geneTree.tree.addModelListener(this);
        }

        compatibleCheckRequired = new boolean[trees.length];
        Arrays.fill(compatibleCheckRequired, false);
        checkCompatibility = false;
    }

    // override this for efficiency, otherwise the overridden makeDirty, which results in additional overhead is called
    public void modelRestored(Model model) {
        super.makeDirty();
    }

    // Upon a direct "make dirty" enable all compatibility checks, since the last call to calculateLogLikelihood may have
    // found a non compatible tree and returned -inf. This case is not explicitly saved.
    public void makeDirty() {
        super.makeDirty(); 
        checkCompatibility = true;
        for(int i = 0; i < mulspb.getGeneTrees().length; i++) {
            compatibleCheckRequired[i] = true;
        }
    }

    protected double calculateLogLikelihood() {
        if( checkCompatibility ) {
            boolean compatibility = true;

            for(int i = 0; i < compatibleCheckRequired.length; ++i) {
                if( compatibleCheckRequired[i] ) {

                    if( !spTree.isCompatible(mulspb.getGeneTrees()[i]) ) {
                        compatibility = false;
                    }
                    compatibleCheckRequired[i] = false;
                    //System.out.println("check compatibility:" + species.getGeneTrees()[i].tree.getId() + " - " + compatibility );
                }
            }
            if( !compatibility ) {
                return Double.NEGATIVE_INFINITY;
            }
            checkCompatibility = false;
        }
        
        // grjtodo-oneday. For now, always do checks. Should make them when do a seq reassign, plus as above, otherwise not
        for (int i  = 0; i < mulspb.getGeneTrees().length; i++) {
        	if (!spTree.isCompatible(mulspb.getGeneTrees()[i])) {
        		return Double.NEGATIVE_INFINITY;
        	}
        }

        double logl = 0;
        int[] info = {0, 0};
        for( MulSpeciesBindings.GeneTreeInfo geneTree : mulspb.getGeneTrees() ) {
            final double v = treeLogLikelihood(geneTree, spTree.getRoot(), info, geneTree.popFactor());
            assert ! Double.isNaN(v);
//           if( Double.isNaN(v) ) {
//               double x = 0;              
//           }
            logl += v;
        }
        ccc += 1;
        return logl;
    }
    
    
    
    int ccc = 0;
    private final boolean verbose = false;

    
    
    private double treeLogLikelihood(MulSpeciesBindings.GeneTreeInfo geneTree, NodeRef node, int[] info, double popFactor) {
        // number of lineages remaining at node
        int nLineages;
        // location in coalescent list (optimization)
        int indexInClist = 0;
        // accumulated log-likelihood inBranchh from node to it's parent
        double like = 0;

        final double t0 = spTree.getNodeHeight(node);

        final MulSpeciesBindings.CoalInfo[] cList = geneTree.getCoalInfo();

        if( verbose ) {
            if( spTree.isRoot(node) ) {
                System.err.println("gtree:" + geneTree.tree.getId());
                System.err.println("t0 " + t0);
                for(int k = 0; k < cList.length; ++k) {
                    System.err.println(k + " " + cList[k].ctime + " " + cList[k].sinfo[0] + " " + cList[k].sinfo[1]);
                }
            }
        }

        if( spTree.isExternal(node) ) {
            nLineages = geneTree.nLineages(spTree.speciesIndex(node));
            indexInClist = 0;
        } else {
            //assert spTree.getChildCount(node) == 2;

            nLineages = 0;
            for(int nc = 0; nc < 2; ++nc) {
                final NodeRef child = spTree.getChild(node, nc);
                like += treeLogLikelihood(geneTree, child, info, popFactor);
                nLineages += info[0];
                indexInClist = Math.max(indexInClist, info[1]);
            }

            // The root of every gene tree (last coalescent point) should be always above
            // root of species tree
            if( indexInClist >=  cList.length)  {
            	System.err.println("ERROR in evomodel.speciation.MulMSCoalescent.treeLogLikelihood()");
            }
            assert indexInClist < cList.length;

            // Skip over (presumably, not tested by assert) non interesting coalescent
            // events to the first event before speciation point

            while( cList[indexInClist].ctime < t0 ) {
                ++indexInClist;
            }
        }

        final boolean isRoot = spTree.isRoot(node);

        // Upper limit
        // use of (t0 + spTree.getBranchLength(node)) caused problem since there was a tiny difference
        // between those (supposedly equal) values. we should track where the discrepancy comes from.
        final double stopTime = isRoot ? Double.MAX_VALUE : spTree.getNodeHeight(spTree.getParent(node));

        // demographic function is 0 based (relative to node height)
        // time away from node
        double lastTime = 0.0;

        // demographic function across branch
        DemographicFunction demog = spTree.getNodeDemographic(node);
        if( popFactor > 0 ) {
            demog = new ScaledDemographic(demog, popFactor);
        }

//        if(false) {
//            final double duration = isRoot ? cList[cList.length - 1].ctime - t0 : stopTime;
//            double demographicAtCoalPoint = demog.getDemographic(duration);
//            double intervalArea = demog.getIntegral(0, duration);
//            if( demographicAtCoalPoint < 1e-12 * (duration/intervalArea) ) {
//               return Double.NEGATIVE_INFINITY;
//            }
//        }
        // Species sharing this branch
        FixedBitSet subspeciesSet = spTree.spSet(node);

        if( verbose ) {
            System.err.println(TreeUtils.uniqueNewick(spTree, node) + " nl " + nLineages
                    + " " + subspeciesSet + " t0 - st " + t0 + " - " + stopTime);
        }

        while( nLineages > 1 ) {
//            if ( !( indexInClist < cList.length ) ) {
//                System.err.println( indexInClist ) ;
//            }
            assert ( indexInClist < cList.length );

            final double nextT = cList[indexInClist].ctime;

            // while rare they can be equal
            if( nextT >= stopTime ) {
                break;
            }

            if( nonEmptyIntersection(cList[indexInClist].sinfo, subspeciesSet) ) {
                final double time = nextT - t0;
                if( time > 0 ) {
                    final double interval = demog.getIntegral(lastTime, time);
                    lastTime = time;

                    final int nLineageOver2 = (nLineages * (nLineages - 1)) / 2;
                    like -= nLineageOver2 * interval;

                    final double pop = demog.getDemographic(time);  assert( pop > 0 );
                    like -= Math.log(pop);
                }
                --nLineages;
            }
            ++indexInClist;
        }


        if( nLineages > 1 ) {
            // add term for No coalescent until root
            final double interval = demog.getIntegral(lastTime, stopTime - t0);

            final int nLineageOver2 = (nLineages * (nLineages - 1)) / 2;

            like -= nLineageOver2 * interval;
        }

        info[0] = nLineages;
        info[1] = indexInClist;
        if( verbose ) {
            System.err.println(TreeUtils.uniqueNewick(spTree, node) + " stopTime " + stopTime +
                    " nl " + nLineages + " icl " + indexInClist);
        }
        return like;
    }

    public void modelChangedEvent(Model model, Object object, int index) {
        //super.modelChangedEvent(model, object, index);
        // not the above for efficiency, otherwise the overridden makeDirty, which results in additional overhead is called.
        super.makeDirty();
        
        if( model == spTree ) {
          if( object == spTree && index != -1 ) {
            // Species tree scaling
            checkCompatibility = true;
            Arrays.fill(compatibleCheckRequired, true);
          }
        } else {

            final MulSpeciesBindings.GeneTreeInfo[] trees = mulspb.getGeneTrees();
            for(int i = 0; i < mulspb.getGeneTrees().length; i++) {
                if( trees[i].tree == model ) {
                    checkCompatibility = true;
                    compatibleCheckRequired[i] = true;
                    break;
                }
            }
        }
    }

    private boolean nonEmptyIntersection(FixedBitSet[] sinfo, FixedBitSet subspeciesSet) {
        for( FixedBitSet nodeSpSet : sinfo ) {
            if( nodeSpSet.intersectCardinality(subspeciesSet) == 0 ) {
                return false;
            }
        }
        return true;
    }


    public Type getUnits() {
        return spTree.getUnits();
    }

    public void setUnits(Type units) {
        assert false;
    }

}
