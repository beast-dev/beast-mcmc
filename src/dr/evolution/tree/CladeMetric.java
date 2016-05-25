/*
 * CladeMetric.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree;

import jebl.evolution.graphs.Node;
import jebl.evolution.taxa.Taxon;
import jebl.evolution.treemetrics.RootedTreeMetric;
import jebl.evolution.trees.RootedTree;

import java.util.*;

/**
 * @author Sebastian Hehna
 * @version 1.0
 * 
 */

public class CladeMetric implements RootedTreeMetric {
   public CladeMetric() {
      taxonMap = null;
   }

   public CladeMetric(List<Taxon> taxa) {
      taxonMap = new HashMap<Taxon, Integer>();
      for (int i = 0; i < taxa.size(); i++) {
         taxonMap.put(taxa.get(i), i);
      }
   }

   public double getMetric(RootedTree tree1, RootedTree tree2) {
      // tree1 is the clade for which we are looking
      // tree2 is the current tree

      Map<Taxon, Integer> tm = taxonMap;

      if (tm == null) {
         List<Taxon> taxa = new ArrayList<Taxon>(tree2.getTaxa());

         if (!tree2.getTaxa().equals(taxa))
            tm = new HashMap<Taxon, Integer>();
         for (int i = 0; i < taxa.size(); i++) {
            tm.put(taxa.get(i), i);
         }
      }

      Clade reference = getClade(tm, tree1, tree1.getRootNode(), null);

      List<Clade> clades = new ArrayList<Clade>();
      getClades(tm, tree2, tree2.getRootNode(), clades, null);

      if (contains(clades, reference)){
         return 1.0;
      }
      else {
         return 0.0;
      }
   }

   private void getClades(Map<Taxon, Integer> taxonMap, RootedTree tree,
         Node node, List<Clade> clades, BitSet bits) {

      BitSet bits2 = new BitSet();

      if (tree.isExternal(node)) {

         int index = taxonMap.get(tree.getTaxon(node));
         bits2.set(index);

      } else {

         for (Node child : tree.getChildren(node)) {
            getClades(taxonMap, tree, child, clades, bits2);
         }
         clades.add(new Clade(bits2, tree.getHeight(node)));
      }

      if (bits != null) {
         bits.or(bits2);
      }
   }
   
   private Clade getClade(Map<Taxon, Integer> taxonMap, RootedTree tree,
         Node node, BitSet bits) {

      BitSet bits2 = new BitSet();
      Clade c = null;

      if (tree.isExternal(node)) {

         int index = taxonMap.get(tree.getTaxon(node));
         bits2.set(index);

      } else {

         for (Node child : tree.getChildren(node)) {
            getClade(taxonMap, tree, child, bits2);
         }
         c = new Clade(bits2, tree.getHeight(node));
      }

      if (bits != null) {
         bits.or(bits2);
      }
      
      return c;
      
   }
   
   private boolean contains (List<Clade> list, Clade c){
      for (Clade cl : list){
         if (c.equals(cl)){
            return true;
         }
      }
      return false;
   }

   BitSet tmpBits = new BitSet();

   private final Map<Taxon, Integer> taxonMap;
}
