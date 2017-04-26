/*
 * MulSpeciesBindings.java
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

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.alloppnet.parsers.MulSpeciesBindingsParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.*;
import dr.math.MathUtils;
import dr.util.HeapSort;
import jebl.util.FixedBitSet;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Binds taxa in gene trees to multiply labelled species.
 * 
 * Extension of Joseph Heled's SpeciesBindings to deal with multiply labelled trees.
 * It combines the SpeciesBindings code with functions from AlloppSpeciesBindings.
 *
 * @author Joseph Heled, Graham Jones
 *         Date: 21/12/2011
 */



public class MulSpeciesBindings extends AbstractModel  implements Loggable {
    // grj all gene trees
    private final GeneTreeInfo[] geneTreeInfos;

    // grj Species definition    
    private final AlloppSpeciesBindings.ApSpInfo[] apspecies;
    private final Taxon[] taxa;
    private final Map<Taxon, Integer> taxon2index = new HashMap<Taxon, Integer>();
	private final int spsq[][];
	private final int numberOfSpSeqs;

    
	// jh
    private final double[][] popTimesPair;
    private boolean dirty_pp;

    private final double[][] popTimesSingle;
    private boolean dirty_sg;
    private final boolean verbose = false;
    
    // grj
	private class SpeciesIndivPair {
		public int spIndex;
		public int ivIndex;
		
		public SpeciesIndivPair(int spIndex, int ivIndex) {
			this.spIndex = spIndex;
			this.ivIndex = ivIndex;    
		}
	}


	// mostly grj
    public MulSpeciesBindings(AlloppSpeciesBindings.ApSpInfo[] apspecies, TreeModel[] geneTrees, double[] popFactors) {
        super(MulSpeciesBindingsParser.MUL_SPECIES);

        this.apspecies = apspecies;
        
        // make the flattened arrays
        int n = 0;
        for (AlloppSpeciesBindings.ApSpInfo apspi : apspecies) {
            n += apspi.individuals.length;
        }
        AlloppSpeciesBindings.Individual [] indivs = new AlloppSpeciesBindings.Individual[n];
        n = 0;
        for (AlloppSpeciesBindings.ApSpInfo apspi : apspecies) {
            for (int i = 0; i < apspi.individuals.length; i++, n++) {
                indivs[n] = apspi.individuals[i];
            }
        }
        int t = 0;
        for (AlloppSpeciesBindings.Individual indiv : indivs) {
            t += indiv.taxa.length;
        }  
        taxa = new Taxon[t];
        t = 0;
        for (AlloppSpeciesBindings.Individual indiv : indivs) {
            for (int j = 0; j < indiv.taxa.length; j++, t++) {
                taxa[t] = indiv.taxa[j];
            }
        }
        // set up maps to indices
        for (int i = 0; i < taxa.length; i++) {
            taxon2index.put(taxa[i], i);
            }
        spsq = new int[apspecies.length][];
        int spsqindex = 0;
        for (int sp = 0; sp < apspecies.length; sp++) {
        	spsq[sp] = new int[apspecies[sp].ploidylevel/2];
        	for (int seq = 0; seq < spsq[sp].length; seq++, spsqindex++) {
        		spsq[sp][seq] = spsqindex;
        	}
        }
        numberOfSpSeqs = spsqindex;

        geneTreeInfos = new GeneTreeInfo[geneTrees.length];
        for (int i = 0; i < geneTrees.length; i++) {
        	final TreeModel gtm = geneTrees[i];
        	addModel(gtm);
        	geneTreeInfos[i] = new GeneTreeInfo(gtm, popFactors[i]);
        }
        
        // like SpeciesBindings but using number of species-sequence pairs, not number of species
        popTimesSingle = new double[numberOfSpSeqs][];
        for (int ns = 0; ns < popTimesSingle.length; ++ns) {
            popTimesSingle[ns] = new double[allCoalPointsCount(ns)];
        }
        dirty_sg = true;

        popTimesPair = new double[(numberOfSpSeqs * (numberOfSpSeqs - 1)) / 2][];
        {
            final int nps = allPairCoalPointsCount();
            for (int ns = 0; ns < popTimesPair.length; ++ns) {
                popTimesPair[ns] = new double[nps];
            }
        }

        dirty_pp = true;

        addStatistic(new SpeciesLimits());
    }




    
    
    // grj
	public int numberOfGeneTrees() {
		return geneTreeInfos.length;
	}	

    
    // grj
    public int nSpSeqs() {
        return numberOfSpSeqs;
    }
    
    // grj
	public String apspeciesName(int i) {
		return apspecies[i].name;
	}
	
	
	//grj
	public int spseqindex2sp(int spsqindex) {
		return spseqindex2spandseq(spsqindex)[0];
	}
	
	//grj
	public int spseqindex2seq(int spsqindex) {
		return spseqindex2spandseq(spsqindex)[1];
	}
	

    
    // grj
	public SpeciesIndivPair apspeciesId2speciesindiv(String apspId) {
		int sp = -1;
		int iv = -1;
		for (int s = 0; s < apspecies.length; s++) {
			for (int i = 0; i < apspecies[s].individuals.length; i++) {
				for (int t = 0; t < apspecies[s].individuals[i].taxa.length; t++) {
					Taxon taxon = apspecies[s].individuals[i].taxa[t];
					if (taxon.getId().compareTo(apspId) == 0) {
						sp = s;
						iv = i;
					}
				}
			}
		}
		assert sp != -1;
		SpeciesIndivPair x = new SpeciesIndivPair(sp, iv);

		return x;
	}
	
	// grj
	public void permuteOneSpeciesOneIndivForOneGene() {
		int i = MathUtils.nextInt(geneTreeInfos.length);
		geneTreeInfos[i].permuteOneSpeciesOneIndiv();
		geneTreeInfos[i].wasChanged();
	}
	
	
	// grj
	public void permuteSetOfIndivsForOneGene() {
		int i = MathUtils.nextInt(geneTreeInfos.length);
		geneTreeInfos[i].permuteSetOfIndivs();
		geneTreeInfos[i].wasChanged();
	}

	
	// grj
	public String seqassignsAsText(int g) {
		return geneTreeInfos[g].seqassignsAsText();
	}

	// grj
	public String genetreeAsText(int g) {
		return geneTreeInfos[g].genetreeAsText();
	}
	
	
   /**
     * Per species coalecent times.
     * <p/>
     * Indexed by sp index, a list of coalescent times of taxa of this sp from all gene trees.
     *
     * @return Per species coalecent times
     */
	// jh
    public double[][] getPopTimesSingle() {
        if (dirty_sg) {
            for (int ns = 0; ns < popTimesSingle.length; ++ns) {
                getAllCoalPoints(ns, popTimesSingle[ns]);
            }
            dirty_sg = false;
        }
        return popTimesSingle;
    }

    // jh
    public double[][] getPopTimesPair() {
        if (dirty_pp) {
            final int nsp = nSpSeqs();
            for (int ns1 = 0; ns1 < nsp - 1; ++ns1) {
                final int z = (ns1 * (2 * nsp - ns1 - 3)) / 2 - 1;

                for (int ns2 = ns1 + 1; ns2 < nsp; ++ns2) {
                    getAllPairCoalPoints(ns1, ns2, popTimesPair[z + ns2]);
                }
            }
        }
        return popTimesPair;
    }

    
    // jh
    private void getAllPairCoalPoints(int ns1, int ns2, double[] popTimes) {

        for (int i = 0; i < geneTreeInfos.length; i++) {
            for (CoalInfo ci : geneTreeInfos[i].getCoalInfo()) {
                if ((ci.sinfo[0].contains(ns1) && ci.sinfo[1].contains(ns2)) ||
                        (ci.sinfo[1].contains(ns1) && ci.sinfo[0].contains(ns2))) {
                    popTimes[i] = ci.ctime;
                    break;
                }
            }
        }
        HeapSort.sort(popTimes);
    }

    
    // jh
    private int allCoalPointsCount(int spseqIndex) {
        int tot = 0;
        for (GeneTreeInfo t : geneTreeInfos) {
            if (t.nLineages(spseqIndex) > 0) {
                tot += t.nLineages(spseqIndex) - 1;
            }
        }
        return tot;
    }

    
    // length of points must be right
    // jh
    void getAllCoalPoints(int spseqIndex, double[] points) {

        int k = 0;
        for (GeneTreeInfo t : geneTreeInfos) {
            final int totCoalEvents = t.nLineages(spseqIndex) - 1;
            int savek = k;
            for (CoalInfo ci : t.getCoalInfo()) {
//               if( ci == null ) {
//                assert ci != null;
//            }
                if (ci.allHas(spseqIndex)) {
                    points[k] = ci.ctime;
                    ++k;
                }
            }
            if (!(totCoalEvents >= 0 && savek + totCoalEvents == k) || (totCoalEvents < 0 && savek == k)) {
                System.err.println(totCoalEvents);
            }
            assert (totCoalEvents >= 0 && savek + totCoalEvents == k) || (totCoalEvents < 0 && savek == k);
        }
        assert k == points.length;
        HeapSort.sort(points);
    }

    
    // jh
    private int allPairCoalPointsCount() {
        return geneTreeInfos.length;
    }

    
    // jh
    public double speciationUpperBound(FixedBitSet sub1, FixedBitSet sub2) {
        //Determined by the last time any pair of sp's in sub1 x sub2 have been seen
        // together in any of the gene trees."""

        double bound = Double.MAX_VALUE;
        for (GeneTreeInfo g : getGeneTrees()) {
            for (CoalInfo ci : g.getCoalInfo()) {
                // if past time of current bound, can't change it anymore
                if (ci.ctime >= bound) {
                    break;
                }
                if ((ci.sinfo[0].intersectCardinality(sub1) > 0 && ci.sinfo[1].intersectCardinality(sub2) > 0)
                        ||
                        (ci.sinfo[0].intersectCardinality(sub2) > 0 && ci.sinfo[1].intersectCardinality(sub1) > 0)) {
                    bound = ci.ctime;
                    break;
                }
            }
        }
        return bound;
    }

    
    // jh
    public void makeCompatible(double rootHeight) {
        for( GeneTreeInfo t : getGeneTrees() ) {

            MutableTree tree = t.tree;

            for (int i = 0; i < tree.getExternalNodeCount(); i++) {
                final NodeRef node = tree.getExternalNode(i);
                final NodeRef p = tree.getParent(node);
                tree.setNodeHeight(p, rootHeight + tree.getNodeHeight(p));
            }
            MutableTree.Utils.correctHeightsForTips(tree);
             // (todo) ugly re-init - can I do something better?
            t.wasChanged();
            t.getCoalInfo();
            t.wasBacked = false;
            //t.wasChanged();
       }
    }


    // jh
    class CoalInfo implements Comparable<CoalInfo> {
        // zero based, 0 is taxa time, i.e. in tree branch units
        final double ctime;
        // sp info for each subtree
        final FixedBitSet[] sinfo;

        CoalInfo(double t, int nc) {
            ctime = t;
            sinfo = new FixedBitSet[nc];
        }

        public int compareTo(CoalInfo o) {
            return o.ctime < ctime ? +1 : (o.ctime > ctime ? -1 : 0);
        }

        /**
         * @param s
         * @return true if all children have at least one taxa from sp 's'
         */
        public boolean allHas(int s) {
            for (FixedBitSet b : sinfo) {
                if (!b.contains(s)) {
                    return false;
                }
            }
            return true;
        }
    }
    
    
    // jh + grj
    /**
     * Collect coalescence information for sub-tree rooted at 'node'.
     *
     * @param tree
     * @param node
     * @param loc  Place node data in loc, sub-tree info before that.
     * @param info array to fill
     * @return location of next available location
     */
    private int collectCoalInfo(Tree tree, NodeRef node, 
    		GeneTreeInfo.SequenceAssignment[] seqassigns, int loc, CoalInfo[] info) {

        info[loc] = new CoalInfo(tree.getNodeHeight(node), tree.getChildCount(node));

        int newLoc = loc - 1;
        for (int i = 0; i < 2; i++) {
            NodeRef child = tree.getChild(node, i);
            info[loc].sinfo[i] = new FixedBitSet(nSpSeqs());

            if (tree.isExternal(child)) {
            	Taxon taxon = tree.getNodeTaxon(child);
                int ti = taxon2index.get(taxon);
				int spseq = spsq[seqassigns[ti].spIndex][seqassigns[ti].seqIndex];
                info[loc].sinfo[i].set(spseq);
                
                assert tree.getNodeHeight(child) == 0;
            } else {
                final int used = collectCoalInfo(tree, child, seqassigns, newLoc, info);
                for (int j = 0; j < info[newLoc].sinfo.length; ++j) {
                    info[loc].sinfo[i].union(info[newLoc].sinfo[j]);
                }
                newLoc = used;
            }
        }
        return newLoc;
    }

    
    // mostly grj
    public class GeneTreeInfo {
        public final TreeModel tree;
    	private SequenceAssignment seqassigns[];
    	private SequenceAssignment oldseqassigns[];
        private final int[] lineagesCount;
        private CoalInfo[] cList;
        private CoalInfo[] savedcList;
        private boolean dirty;
        private boolean wasBacked;
        private final double popFactor;
        
        
    	/* class GeneTreeInfo.SequenceAssignments
    	 * 
    	 * spIndex is an index for an allopolyploid species. For example, it identifies
    	 * a bit in a FixedBitSet (union) in a MulLabTree
    	 * 
    	 * seqIndex identifies a sequence copy for this gene and for each individual.
    	 * seqIndex is 0 or 1 for tetraploids and it is these that get flipped to change
    	 * assignments of sequence copies to legs in AlloppSpeciesNetworkModel (or
    	 * equivalently to tips in a MulLabTree).
    	 * 
    	 * 2011-06-23 spIndex is the same for all gene trees. Maybe
    	 * allow non-rectangular data later.
    	 */
    	private class SequenceAssignment {
     		public int spIndex;
    		public int seqIndex;
    		
    		public SequenceAssignment(int spIndex, int seqIndex) {
    			this.spIndex = spIndex;
    			this.seqIndex = seqIndex;
    		}
    		
    		public String toString() {
    			String s = "" + seqIndex;
    			return s;
    		}
    	} 	

    	
    	
    	// grj
        GeneTreeInfo(TreeModel tree, double popFactor) {
            this.tree = tree;
            this.popFactor = popFactor;
            seqassigns = new SequenceAssignment[taxa.length];
            oldseqassigns = new SequenceAssignment[taxa.length];      
            
            // This uses taxa list for *all* gene trees, not this gene tree.
            for (int s = 0; s < apspecies.length; s++) {
            	for (int i = 0; i < apspecies[s].individuals.length; i++) {
            		int nseqs = apspecies[s].individuals[i].taxa.length;
            		int asgns[] = new int [nseqs];
            		for (int x = 0; x < nseqs; x++) {
            			asgns[x] = x;
            		}
            		MathUtils.permute(asgns);
            		for (int x = 0; x < nseqs; x++) {
            			int t = taxon2index.get(apspecies[s].individuals[i].taxa[x]);
            			seqassigns[t] = new SequenceAssignment(s, asgns[x]);
            			oldseqassigns[t] = new SequenceAssignment(s, asgns[x]);
            		}
            	}
            }
            lineagesCount = new int[nSpSeqs()];
            Arrays.fill(lineagesCount, 0);

            for (int nl = 0; nl < lineagesCount.length; ++nl) {
            	int sp = spseqindex2sp(nl);            	
                for (AlloppSpeciesBindings.Individual indiv : apspecies[sp].individuals) {
                	boolean got = false;
                	for (Taxon t : indiv.taxa) {
                		if (tree.getTaxonIndex(t) >= 0) {
                			got = true;
                		}
                	}
                	for (Taxon t : indiv.taxa) {
                		assert (tree.getTaxonIndex(t) >= 0) == got;
                		}
                	assert got;
                    if (got) {
                        ++lineagesCount[nl];
                	}                	
                }
            }
            
            // this bit jh
            cList = new CoalInfo[tree.getExternalNodeCount() - 1];
            savedcList = new CoalInfo[cList.length];
            wasChanged();
            getCoalInfo();
            wasBacked = false;
        }
        
        
        // grj
        public String seqassignsAsText() {
        	String s = "Sequence assignments" + System.getProperty("line.separator");
        	for (int tx = 0; tx < seqassigns.length; tx++) {
        		s += taxa[tx];
        		s += ":";
        		s += seqassigns[tx].seqIndex;
        		if (tx+1 < seqassigns.length  &&  seqassigns[tx].spIndex != seqassigns[tx+1].spIndex) {
        			s += System.getProperty("line.separator");
        		} else {
        			s += "  ";
        		}
        	}
        	return s;
        }
        
        // grj
         public String genetreeAsText() {
        	 return tree.getNewick();
         }
        
        // grj
        public SequenceAssignment getSeqassigns(int tx) {
        	return seqassigns[tx];
        }
        
        
        // grj
		public void storeSequenceAssignments() {
			for (int i = 0; i < seqassigns.length; i++) {
				oldseqassigns[i].seqIndex = seqassigns[i].seqIndex;
			}
		}
        
        // grj
		public void restoreSequenceAssignments() {
			for (int i = 0; i < seqassigns.length; i++) {
				seqassigns[i].seqIndex = oldseqassigns[i].seqIndex;
			}
		}

        
        // grj
    	public void permuteOneSpeciesOneIndiv() {
    		int sp = MathUtils.nextInt(apspecies.length);
    		int iv = MathUtils.nextInt(apspecies[sp].individuals.length);
    		permuteOneAssignment(sp, iv);
    	}       

        
    	/* grjtodo-oneday
    	 * This is a bit odd. It collects individuals as (sp, iv) indices
    	 * that `belong' to a node in the sense that any taxon (sequence)
    	 * of an individual belongs to the clade of the node.
    	 * I've used a set but not made SpeciesIndivPair's comparable
    	 * so that if both sequences of an individual occurs in clade it appears
    	 * twice. Then permuteOneAssignment() flips everything so that those
    	 * occurring twice get flipped twice and so not changed.
    	 * 
    	 * Result is that individuals with one but not two sequences in
    	 * the clade of the node get flipped. Sometimes all individuals
    	 * are flipped, sometimes none, sometimes just one, the last is the 
    	 * same as permuteOneSpeciesOneIndiv().
    	 * 
    	 * 2011-07-29 it appears to work OK on minimal testing and I
    	 * don't have a good idea for a more rational or efficient version. 
    	 * 
    	 */
    	public void permuteSetOfIndivs() {
    		int num = tree.getInternalNodeCount();
    		int i = MathUtils.nextInt(num);	
    		NodeRef node = tree.getInternalNode(i);
    		Set<SpeciesIndivPair> spivs = new HashSet<SpeciesIndivPair>();
    		collectIndivsOfNode(node, spivs);
    		for (SpeciesIndivPair spiv : spivs) {
    			permuteOneAssignment(spiv.spIndex, spiv.ivIndex);
    		}
    	}        
    	

    	// jh
        int nLineages(int spseqIndex) {
            return lineagesCount[spseqIndex];
        }

        // jh
        public CoalInfo[] getCoalInfo() {
            if (dirty) {
                swap();

                collectCoalInfo(tree, tree.getRoot(), seqassigns, cList.length - 1, cList);
                HeapSort.sort(cList);
                dirty = false;
                wasBacked = true;
            }
            
            /*
            CoalInfo check[] = new CoalInfo[cList.length];
            collectCoalInfo(tree, tree.getRoot(), seqassigns, check.length - 1, check);
            HeapSort.sort(check);
            for (int i=0; i<check.length; ++i) {
            	if (check[i].ctime != cList[i].ctime) {
            		System.err.println("Inconsistent ctimes " + i + " check " + check[i].ctime+ " cList " + cList[i].ctime);
            	}
            	if (!check[i].sinfo[0].equals(cList[i].sinfo[0])) {
            		System.err.println("Inconsistent sinfo[0] " + i );
            	}
            	if (!check[i].sinfo[1].equals(cList[i].sinfo[1])) {
            		System.err.println("Inconsistent sinfo[1] " + i );
            	}
            }*/
            return cList;
        }

        // jh
        private void swap() {
            CoalInfo[] tmp = cList;
            cList = savedcList;
            savedcList = tmp;
        }

        // jh
        void wasChanged() {
            dirty = true;
            wasBacked = false;
        }

        // jh
        boolean restore() {
            if (verbose) System.out.println(" SP binding: restore " + tree.getId() + " (" + wasBacked + ")");
            if (wasBacked) {
//                if( false ) {
//                    swap();
//                    dirty = true;
//                    getCoalInfo();
//                    for(int k = 0; k < cList.length; ++k) {
//                        assert cList[k].ctime == savedcList[k].ctime &&
//                                cList[k].sinfo[0].equals(savedcList[k].sinfo[0]) &&
//                                cList[k].sinfo[1].equals(savedcList[k].sinfo[1]);
//                    }
//                }
                swap();
                wasBacked = false;
                dirty = false;
                return true;
            }
            return false;
        }

        
        // jh
        void accept() {
            if (verbose) System.out.println(" SP binding: accept " + tree.getId());

            wasBacked = false;
        }

        
        // jh
        public double popFactor() {
            return popFactor;
        }
        
        
        // grj
		private void collectIndivsOfNode(NodeRef node, Set<SpeciesIndivPair> spivs) {
			if (tree.isExternal(node)) {
				SpeciesIndivPair x = apspeciesId2speciesindiv(tree.getNodeTaxon(node).getId());
				spivs.add(x);
			} else {
				collectIndivsOfNode(tree.getChild(node, 0), spivs);
				collectIndivsOfNode(tree.getChild(node, 1), spivs);
			}
		}
		
        
        // grj
		private void permuteOneAssignment(int sp, int iv) {
			// grjtodo-tetraonly
			int tx;
			if (apspecies[sp].individuals[iv].taxa.length == 2) {
				tx = taxon2index.get(apspecies[sp].individuals[iv].taxa[0]);
				seqassigns[tx].seqIndex = 1 - seqassigns[tx].seqIndex;
				tx = taxon2index.get(apspecies[sp].individuals[iv].taxa[1]);
				seqassigns[tx].seqIndex = 1 - seqassigns[tx].seqIndex;
			}
		}

    }  // end of GeneTreeInfo

    
    // jh
    public GeneTreeInfo[] getGeneTrees() {
        return geneTreeInfos;
    }

    
    // jh + grj
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (verbose) System.out.println(" SP binding: model changed " + model.getId());

        dirty_sg = true;
        dirty_pp = true;

        for (GeneTreeInfo g : geneTreeInfos) {
            if (g.tree == model) {
                g.wasChanged();
                break;
            }
        }
        fireModelChanged(object, index);
    }

    // jh
    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        assert false;
    }

     // grj
    // jh comment was 'do on a per need basis'. I hope its ok to mix with doing always.
   protected void storeState() {
	for (GeneTreeInfo gti : geneTreeInfos) {
		gti.storeSequenceAssignments();
	}
	if (MulSpeciesTreeModel.DBUGTUNE)
		System.err.println("MulSpeciesBindings.storeState()");
    }

   

   
    // jh + grj
    protected void restoreState() {
    	for (GeneTreeInfo gti : geneTreeInfos) {
    		gti.restoreSequenceAssignments();
    	}
        for (GeneTreeInfo g : geneTreeInfos) {
            if (g.restore()) {
                dirty_sg = true;
                dirty_pp = true;
            }
        }
    	if (MulSpeciesTreeModel.DBUGTUNE)
    		System.err.println("MulSpeciesBindings.restoreState()");
    }

    
    // jh + grj
    protected void acceptState() {
        for (GeneTreeInfo g : geneTreeInfos) {
            g.accept();
        }
    }
    
    
    // grj
	public LogColumn[] getColumns() {
		int ncols = geneTreeInfos.length * taxa.length;
		LogColumn[] columns = new LogColumn[ncols];
		for (int g = 0, i = 0; g < geneTreeInfos.length; g++) {
			for (int tx = 0; tx < taxa.length; tx++, i++) {
				GeneTreeInfo.SequenceAssignment sqa = geneTreeInfos[g].getSeqassigns(tx);
				String header = "Gene" + g + "taxon" + tx; 
				columns[i] = new LogColumn.Default(header, sqa);
			}
				
		}
		
		return columns;
	}

    // jh
    public class SpeciesLimits extends Statistic.Abstract {
        int nDim;
        int c[][];

        SpeciesLimits() {
            super("SpeciationBounds");

            nDim = 0;

            final int nsp = nSpSeqs();

            c = new int[nsp + 1][nsp + 1];
            for(int k = 0; k < nsp + 1; ++k) {
                c[k][0] = 1;
                c[k][k] = 1;
            }
            for(int k = 0; k < nsp + 1; ++k) {
                for(int j = 1; j < k; ++j) {
                    c[k][j] = c[k - 1][j - 1] + c[k - 1][j];
                }
            }

            for(int k = 0; k <= (int) (nsp / 2); ++k) {
                nDim += c[nsp][k];
            }

        }

        // jh
        public int getDimension() {
            return nDim;
        }

        // jh
        private double boundOnRoot() {
            double bound = Double.MAX_VALUE;
            final int nsp = nSpSeqs();
            for(GeneTreeInfo g : getGeneTrees()) {
                for(CoalInfo ci : g.getCoalInfo()) {
                    if( ci.sinfo[0].cardinality() == nsp || ci.sinfo[1].cardinality() == nsp ) {
                        bound = Math.min(bound, ci.ctime);
                        break;
                    }
                }
            }
            return bound;
        }

        // jh
        public double getStatisticValue(int dim) {
            if( dim == 0 ) {
                return boundOnRoot();
            }

            final int nsp = nSpSeqs();
            int r = 0;
            int k;
            for(k = 0; k <= (int) (nsp / 2); ++k) {
                final int i = c[nsp][k];
                if( dim < r + i ) {
                    break;
                }
                r += i;
            }

            // Classic index -> select k of nsp subset

            // number of species in set is k
            int n = dim - r;
            FixedBitSet in = new FixedBitSet(nsp),
                    out = new FixedBitSet(nsp);
            int fr = nsp;
            for(int i = 0; i < nsp; ++i) {
                if( k == 0 ) {
                    out.set(i);
                } else {
                    if( n < c[fr - 1][k - 1] ) {
                        in.set(i);
                        k -= 1;
                    } else {
                        out.set(i);
                        n -= c[fr - 1][k];
                    }
                    fr -= 1;
                }
            }
            return speciationUpperBound(in, out);
        }
    }
    
    
    
    // grj
	private int[] spseqindex2spandseq(int spsqindex) {
		int indexp = -1;
		int indexq = -1;
		for (int p = 0; p < spsq.length; p++) {
			for (int q = 0; q < spsq[p].length; q++) {
				if (spsq[p][q] == spsqindex) {
				assert indexp == -1;
				assert indexq == -1;
				indexp = p;
				indexq = q;
				}
			} 
		}
		assert indexp != -1;
		assert indexq != -1;
		int[] pq = new int[2];
		pq[0] = indexp;
		pq[1] = indexq;
		return pq;
	}	

}