/**
 * 
 */
package dr.evolution.tree;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import jebl.evolution.taxa.Taxon;
import jebl.evolution.graphs.Node;
import jebl.evolution.treemetrics.RootedTreeMetric;
import jebl.evolution.trees.RootedTree;

/**
 * @author Sebastian Hoehna
 * @version 1.0
 * 
 */
public class BranchScoreMetric implements RootedTreeMetric {

   public BranchScoreMetric() {
      taxonMap = null;
   }

   public BranchScoreMetric(List<Taxon> taxa) {
      taxonMap = new HashMap<Taxon, Integer>();
      for (int i = 0; i < taxa.size(); i++) {
         taxonMap.put(taxa.get(i), i);
      }
   }

   public double getMetric(RootedTree tree1, RootedTree tree2) {

      if (!tree1.getTaxa().equals(tree2.getTaxa())) {
         throw new IllegalArgumentException("Trees contain different taxa");
      }

      Map<Taxon, Integer> tm = taxonMap;

      if (tm == null) {
         List<Taxon> taxa = new ArrayList<Taxon>(tree1.getTaxa());

         if (!tree2.getTaxa().equals(taxa))
            tm = new HashMap<Taxon, Integer>();
         for (int i = 0; i < taxa.size(); i++) {
            tm.put(taxa.get(i), i);
         }
      }

      List<Clade> clades1 = new ArrayList<Clade>();
      getClades(tm, tree1, tree1.getRootNode(), clades1, null);

      List<Clade> clades2 = new ArrayList<Clade>();
      getClades(tm, tree2, tree2.getRootNode(), clades2, null);

      return getDistance(clades1, clades2);
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

   private double getDistance(List<Clade> clades1, List<Clade> clades2) {

      Collections.sort(clades1);
      Collections.sort(clades2);
      double distance = 0.0;
      int indexClade2 = 0;
      Clade clade2 = null;
      Clade parent1, parent2 = null;
      double height1, height2;

      for (Clade clade1 : clades1) {
         
         parent1 = findParent(clade1, clades1);
         height1 = parent1.getHeight() - clade1.getHeight();

         if (indexClade2 < clades2.size()) {
            clade2 = clades2.get(indexClade2);
            parent2 = findParent(clade2, clades2);
         }
         while (clade1.compareTo(clade2) > 0 && indexClade2 < clades2.size()) {
            height2 = parent2.getHeight() - clade2.getHeight();
            distance += height2 * height2;
            indexClade2++;
            if (indexClade2 < clades2.size()) {
               clade2 = clades2.get(indexClade2);
               parent2 = findParent(clade2, clades2);
            }
         }
         if (clade1.compareTo(clade2) == 0) {
            height2 = parent2.getHeight() - clade2.getHeight();
            distance += (height1 - height2) * (height1 - height2);
            indexClade2++;
         } else {
            distance += height1 * height1;
         }
      }

      return Math.sqrt(distance);
   }

   private Clade findParent(Clade clade1, List<Clade> clades) {
      Clade parent = null;
      for (Clade clade2 : clades) {
         if (isParent(clade2, clade1)) {
            if (parent == null || parent.getSize() > clade2.getSize())
            parent = clade2;
         }
      }
      
      if (parent == null){
         //the case that this clade is the whole tree 
         return clade1;
      }

      return parent;
   }

   private boolean isParent(Clade parent, Clade child) {
      if (parent.getSize() <= child.getSize()) {
         return false;
      }

      tmpBits.clear();
      tmpBits.or(parent.getBits());
      tmpBits.xor(child.getBits());

      return tmpBits.cardinality() < parent.getSize();
   }

   BitSet                            tmpBits = new BitSet();

   private final Map<Taxon, Integer> taxonMap;

}
