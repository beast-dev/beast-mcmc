/**
 * 
 */
package dr.inference.markovchain;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import dr.evolution.tree.NodeRef;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.CompoundModel;
import dr.inference.model.Model;

/**
 * @author Sebastian Hoehna
 * 
 */
public class ConvergenceListener implements MarkovChainListener {

   private MarkovChain mc = null;

   private int checkEvery = 0;

   private double threshold = 0.0;

   private BufferedWriter bwDistances = null;

   private long trees = 0;

   private HashMap<BitSet, Double> referenceSplitFrequencies = null;

   private HashMap<BitSet, Double> splitOccurences = null;

   private String outputFilename;

   private HashMap<String, Integer> taxonMap;

   /**
    * 
    */
   public ConvergenceListener(MarkovChain mc, int checkEvery,
         double threshold, String outputFilename,
         String referenceTreeFileName) {
      this.mc = mc;
      this.checkEvery = checkEvery;
      this.threshold = threshold;
      this.outputFilename = outputFilename.substring(0, outputFilename.indexOf('.')) + ".dist";
      
      File f = new File(outputFilename.substring(0, outputFilename.indexOf('.')) + ".dist");
      if (f.exists()){
         f.delete();
      }

      
      taxonMap = new HashMap<String, Integer>();
      splitOccurences = new HashMap<BitSet, Double>();

      parseReferenceFile(referenceTreeFileName);
   }

   private void parseReferenceFile(String log) {
      FileReader fr;
      BufferedReader br;
      String line;
      boolean cladeScetion = false;
      referenceSplitFrequencies = new HashMap<BitSet, Double>();

      cladeScetion = false;
      try {
         fr = new FileReader(log);
         br = new BufferedReader(fr);

         do {
            line = br.readLine();
            if (line != null) {
               line = line.trim();
               if (cladeScetion && !line.equals("")) {
                  String[] taxa = getTaxa(line);
                  BitSet split = new BitSet();
                  for (String t : taxa) {
                     if (t.indexOf(',') != -1) {
                        t = t.substring(0, t.indexOf(','));
                     }
                     if (!taxonMap.containsKey(t)) {
                        taxonMap.put(t, taxonMap.size());
                     }
                     split.set(taxonMap.get(t));
                  }
                  referenceSplitFrequencies.put(split,
                        getPercentage(line));
               } else if (line.contains("%-rule clades")) {
                  cladeScetion = true;
               }
            }
         } while (line != null);

         br.close();
         fr.close();
      } catch (FileNotFoundException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      } catch (IOException e) {
         // TODO Auto-generated catch block
         e.printStackTrace();
      }
   }

   private double getPercentage(String line) {
      String[] tokens = line.split("\\s++");
      return Double.parseDouble(tokens[1]
            .substring(0, tokens[1].length() - 2));
   }

   private String[] getTaxa(String line) {
      line = line.substring(line.indexOf('{') + 1, line.indexOf('}'));
      String[] taxa = line.split("\\s++");
      return taxa;
   }

   private double getMaxDeviation() {
      double[] deviation = getDeviations();

      double max = 0;
      for (double m : deviation) {
         if (m > max) {
            max = m;
         }
      }
      return max;
   }

   private double[] getDeviations() {
      double[] deviations = new double[Math.max(referenceSplitFrequencies
            .size(), splitOccurences.size())];

      Set<BitSet> keys = null;
      if (referenceSplitFrequencies.size() < splitOccurences.size()) {
         keys = splitOccurences.keySet();
      } else {
         keys = referenceSplitFrequencies.keySet();
      }

      Iterator<BitSet> it = keys.iterator();
      for (int i = 0; i < deviations.length; i++) {
         BitSet k = it.next();
         double r = 0.0;
         if (referenceSplitFrequencies.containsKey(k)){
            r = referenceSplitFrequencies.get(k);
         }
         double o = 0.0;
         if (splitOccurences.containsKey(k)){
            o = splitOccurences.get(k);
         }
         deviations[i] = Math.abs(r
               - ((o*100.0) / (double) trees));
      }

      return deviations;
   }

   private void addTree(TreeModel tree) {
      trees++;
      List<BitSet> splits = getSplits(tree);

      for (BitSet s : splits) {
         if (!splitOccurences.containsKey(s)) {
            splitOccurences.put(s, 1.0);
         } else {
            splitOccurences.put(s, splitOccurences.get(s) + 1);
         }
      }
   }

   private List<BitSet> getSplits(TreeModel tree) {
      List<BitSet> splits = new ArrayList<BitSet>();

      NodeRef root = tree.getRoot();
      fillSplits(splits, root, tree);

      return splits;
   }

   private BitSet fillSplits(List<BitSet> splits, NodeRef root, TreeModel tree) {
      BitSet split = new BitSet();
      

      if (!tree.isExternal(root)) {
         split.or(fillSplits(splits, tree.getChild(root, 0), tree));
         split.or(fillSplits(splits, tree.getChild(root, 1), tree));
         splits.add(split);
      }
      else {
         Taxon taxon = tree.getNodeTaxon(root);
         String name = taxon.getId();
         split.set(taxonMap.get(name));
      }

      return split;

   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.markovchain.MarkovChainListener#bestState(int,
    *      dr.inference.model.Model)
    */
   public void bestState(int state, Model bestModel) {
   }

   private TreeModel getTreeModel(Model m) {
      if (m instanceof CompoundModel) {
         Model m1 = getTreeModelFromCompoundModel((CompoundModel) m);
         if (m1 instanceof TreeModel) {
            return (TreeModel) m1;
         }
      } else if (m instanceof TreeModel) {
         return (TreeModel) m;
      }
      return null;
   }

   private TreeModel getTreeModelFromCompoundModel(CompoundModel m) {

      for (int i = 0; i < m.getModelCount(); i++) {
         Object m1 = m.getModel(i);
         if (m1 instanceof CompoundModel) {
            Model m2 = getTreeModelFromCompoundModel((CompoundModel) m1);
            if (m2 instanceof TreeModel) {
               return (TreeModel) m2;
            }
         } else if (m1 instanceof TreeLikelihood) {
            Model m2 = getTreeModelFromTreeLikelihood((TreeLikelihood) m1);
            if (m2 instanceof TreeModel) {
               return (TreeModel) m2;
            }
         } else if (m1 instanceof TreeModel) {
            return (TreeModel) m1;
         }
      }
      return null;
   }

   private TreeModel getTreeModelFromTreeLikelihood(TreeLikelihood tl) {

      for (int i = 0; i < tl.getModelCount(); i++) {
         Object m1 = tl.getModel(i);
         if (m1 instanceof TreeModel) {
            return (TreeModel) m1;
         }
      }
      return null;
   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.markovchain.MarkovChainListener#currentState(int,
    *      dr.inference.model.Model)
    */
   public void currentState(int state, Model currentModel) {

      if (state > 10 * checkEvery && checkEvery > 0
            && (state % checkEvery == 0)) {

         TreeModel tree = getTreeModel(currentModel);

         addTree(tree);

         if (state % (checkEvery * 10) == 0) {
            double distance = getMaxCladeDistance();
            try {
               bwDistances = new BufferedWriter(new FileWriter(
                     outputFilename, true));
               bwDistances.append((state) + ":\t" + distance);
               bwDistances.newLine();
               bwDistances.flush();
            } catch (IOException e) {
               // TODO Auto-generated catch block
               e.printStackTrace();
            }
            if (distance <= threshold) {
               try {
                  bwDistances.close();
               } catch (IOException e) {
                  // TODO Auto-generated catch block
                  e.printStackTrace();
               }
               fireMCStopp();
            }
         }
      }
   }

   private void fireMCStopp() {
      mc.pleaseStop();
   }

   private double getMaxCladeDistance() {

      return getMaxDeviation();
   }

   public void setMarkovChain(MarkovChain mc) {
      this.mc = mc;
   }

   /*
    * (non-Javadoc)
    * 
    * @see dr.inference.markovchain.MarkovChainListener#finished(int)
    */
   public void finished(int chainLength) {
      // TODO Auto-generated method stub

   }

}
