package dr.evomodel.speciation.agedependent.simulation;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.speciation.agedependent.agehazard.AgeHazard;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Frederik M. Andersen
 *
 * Simulates a symmetric or asymmetric time- and age-dependent birth-death tree using a thinning
 * algorithm, and returns the reconstructed tree over the sampled tips.
 *
 * Sampling is serial: sampled-through-time tips are produced by psi events during the simulation
 * (the lineage is removed on sampling), and lineages surviving to the present are each sampled
 * with probability rho. Setting psi = 0 and rho = 1 recovers the ultrametric process, in which
 * every lineage extant at the present becomes a tip.
 *
 * Rates:
 *     lambda(t, a) = birthScale(t) * h_b(a)
 *     mu(t, a)     = deathScale(t) * h_d(a)
 *     psi(t)       = samplingScale(t)
 *
 * where h_b and h_d are supplied as {@link AgeHazard} instances, the same pluggable hazard
 * family the likelihood models use, so a simulation and the likelihood it is validated against
 * can be driven by the identical hazard shape. A hazard covering one epoch is shared across
 * all epochs; one covering numEpochs supplies a distinct shape per epoch.
 *
 * Thinning uses {@link AgeHazard#maxHazard} as the per-epoch envelope, so any hazard providing
 * a correct upper bound over [0, originTime] can be simulated under without further changes.
 */
public class AgeDependentBirthDeathSimulator {
    private final double[] birthScale;
    private final double[] deathScale;
    private final double[] samplingScale;
    private final double extantSamplingProb;
    private final AgeHazard birthHazard;
    private final AgeHazard deathHazard;
    private final double[] epochBounds;
    private final double originTime;
    private final boolean symmetric;
    private final int maxLineages;

    private final double[] maxBirthHaz;
    private final double[] maxDeathHaz;

    /**
     * @param birthScale         piecewise-constant birth scale, one per epoch
     * @param deathScale         piecewise-constant death scale, one per epoch (or length 1)
     * @param samplingScale      piecewise-constant serial sampling rate psi(t), one per epoch (or length 1)
     * @param extantSamplingProb extant sampling probability rho in [0, 1]
     * @param birthHazard        age-hazard h_b(a, epoch); must cover 1 or numEpochs epochs
     * @param deathHazard        age-hazard h_d(a, epoch); must cover 1 or numEpochs epochs
     * @param epochTimes         internal epoch boundaries in backwards time (ascending); excludes origin
     * @param originTime         the origin time (most ancient point)
     * @param symmetric          if true, both daughters get age 0; if false, one inherits parent age
     * @param maxLineages        safety cap to prevent runaway growth
     */
    public AgeDependentBirthDeathSimulator(double[] birthScale,
                                                 double[] deathScale,
                                                 double[] samplingScale,
                                                 double extantSamplingProb,
                                                 AgeHazard birthHazard,
                                                 AgeHazard deathHazard,
                                                 double[] epochTimes,
                                                 double originTime,
                                                 boolean symmetric,
                                                 int maxLineages) {
        this.birthScale = birthScale;
        this.deathScale = deathScale;
        this.samplingScale = samplingScale;
        this.extantSamplingProb = extantSamplingProb;
        int numEpochs = epochTimes.length + 1;
        checkEpochCount(birthHazard, "birthHazard", numEpochs);
        checkEpochCount(deathHazard, "deathHazard", numEpochs);
        this.birthHazard = birthHazard;
        this.deathHazard = deathHazard;
        this.epochBounds = new double[epochTimes.length + 2];
        this.epochBounds[0] = 0.0;
        System.arraycopy(epochTimes, 0, this.epochBounds, 1, epochTimes.length);
        this.epochBounds[epochTimes.length + 1] = originTime;
        this.originTime = originTime;
        this.symmetric = symmetric;
        this.maxLineages = maxLineages;

        // Thinning envelope: one bound per epoch, evaluated over the whole age range a lineage
        // can reach. A hazard covering a single epoch returns the same bound for every k.
        this.maxBirthHaz = new double[numEpochs];
        this.maxDeathHaz = new double[numEpochs];
        for (int k = 0; k < numEpochs; k++) {
            maxBirthHaz[k] = birthHazard.maxHazard(originTime, hazardEpoch(birthHazard, k));
            maxDeathHaz[k] = deathHazard.maxHazard(originTime, hazardEpoch(deathHazard, k));
        }
    }

    /**
     * A hazard covering a single epoch is shared across all of them; otherwise the epoch
     * selects its own shape. Mirrors {@code AgeHazard.idx} at the hazard-object level.
     */
    private static int hazardEpoch(AgeHazard hazard, int epoch) {
        return hazard.getEpochCount() == 1 ? 0 : epoch;
    }

    private static void checkEpochCount(AgeHazard hazard, String name, int numEpochs) {
        int n = hazard.getEpochCount();
        if (n != 1 && n != numEpochs) {
            throw new IllegalArgumentException(name + " must cover 1 (shared) or " + numEpochs
                    + " epochs, got " + n);
        }
    }

    /**
     * Simulate a tree, retrying until the sampled-tip count is in [minTips, maxTips].
     */
    public Tree simulate(int minTips, int maxTips, int maxAttempts) {
        for (int attempt = 0; attempt < maxAttempts; attempt++) {
            Tree tree = simulateOnce();

            if (tree != null) {
                int n = tree.getExternalNodeCount();

                if (n >= minTips && (maxTips <= 0 || n <= maxTips)) {
                    Logger.getLogger("dr.evomodel.speciation").info(
                            "Simulated TABD tree with " + n +
                                    " tips (attempt " + (attempt + 1) + ")");
                    return tree;
                }
            }
        }
        String range = maxTips > 0 ? "[" + minTips + ", " + maxTips + "]" : ">= " + minTips;
        throw new RuntimeException(
                "Failed to simulate TABD tree with " + range +
                        " tips after " + maxAttempts + " attempts");
    }

    public Tree simulate(int minTips, int maxAttempts) {
        return simulate(minTips, 0, maxAttempts);
    }

    /**
     * Single simulation attempt
     */
    private Tree simulateOnce() {
        List<LineAge> lineages = new ArrayList<>();
        List<FlexibleNode> sampledTips = new ArrayList<>();
        FlexibleNode stemNode = new FlexibleNode();
        lineages.add(new LineAge(stemNode, 0.0));

        double currentTime = originTime;
        int epoch = epochBounds.length - 2;
        int taxonCount = 0;

        while (currentTime > 0 && !lineages.isEmpty()) {
            int K = lineages.size();
            if (K > maxLineages) {
                return null;
            }

            double bScale = birthScale[epoch];
            double dScale = (deathScale.length == 1) ? deathScale[0] : deathScale[epoch];
            double sScale = (samplingScale.length == 1) ? samplingScale[0] : samplingScale[epoch];
            double epochBound = epochBounds[epoch];

            double maxRate = K * (bScale * maxBirthHaz[epoch] + dScale * maxDeathHaz[epoch] + sScale);
            if (maxRate <= 0.0) {
                double dt = currentTime - epochBound;
                for (LineAge l : lineages) l.age += dt;
                currentTime = epochBound;
                if (epoch > 0) epoch--;
                continue;
            }

            double prop = -Math.log(MathUtils.nextDouble()) / maxRate;
            if (currentTime - prop <= epochBound) {
                double dt = currentTime - epochBound;
                for (LineAge l : lineages) {
                    l.age += dt;
                }
                currentTime = epochBound;
                if (epoch > 0) {
                    epoch--;
                }
                continue;
            }

            currentTime -= prop;
            for (LineAge l : lineages) {
                l.age += prop;
            }

            int idx = MathUtils.nextInt(K);
            LineAge chosen = lineages.get(idx);

            double lam = bScale * birthHazard.evaluate(chosen.age, hazardEpoch(birthHazard, epoch));
            double mu  = dScale * deathHazard.evaluate(chosen.age, hazardEpoch(deathHazard, epoch));
            double psi = sScale;
            double r = lam + mu + psi;

            if (MathUtils.nextDouble() >= K * r / maxRate) {
                continue;
            }

            // Choose event: birth / death / sampling
            double u = MathUtils.nextDouble() * r;
            if (u < lam) {
                FlexibleNode parent = chosen.node;
                parent.setHeight(currentTime);

                FlexibleNode child1 = new FlexibleNode();
                FlexibleNode child2 = new FlexibleNode();
                parent.addChild(child1);
                parent.addChild(child2);

                if (symmetric) {
                    lineages.set(idx, new LineAge(child1, 0.0));
                    lineages.add(new LineAge(child2, 0.0));
                } else {
                    lineages.set(idx, new LineAge(child1, chosen.age));
                    lineages.add(new LineAge(child2, 0.0));
                }
            } else if (u < lam + mu) {
                FlexibleNode deadNode = chosen.node;
                deadNode.setHeight(currentTime);
                removeAt(lineages, idx);
            } else {
                // Sampling event: the lineage becomes a sampled tip at currentTime > 0.
                FlexibleNode sampled = chosen.node;
                sampled.setHeight(currentTime);
                taxonCount++;
                sampled.setTaxon(new Taxon("taxon" + taxonCount));
                sampledTips.add(sampled);
                removeAt(lineages, idx);
            }
        }

        // Extant sampling: at present, each surviving lineage is sampled with prob rho.
        for (LineAge l : lineages) {
            l.node.setHeight(0.0);
            if (extantSamplingProb >= 1.0 || MathUtils.nextDouble() < extantSamplingProb) {
                taxonCount++;
                l.node.setTaxon(new Taxon("taxon" + taxonCount));
                sampledTips.add(l.node);
            }
        }

        if (sampledTips.isEmpty()) {
            return null;
        }

        FlexibleNode prunedRoot = reconstructTree(stemNode);

        if (prunedRoot == null || prunedRoot.getChildCount() == 0) {
            return null;
        }

        // The simulation only ever assigns node heights, never branch lengths. The single-arg
        // FlexibleTree constructor asserts BOTH are known, so every getBranchLength() would
        // silently return the unset default of 0. Declaring lengths unknown makes FlexibleTree
        // derive them from the heights on first access.
        return new FlexibleTree(prunedRoot, true, false);
    }

    private static void removeAt(List<LineAge> lineages, int idx) {
        int last = lineages.size() - 1;
        if (idx != last) {
            lineages.set(idx, lineages.get(last));
        }
        lineages.remove(last);
    }

    /**
     * Recursively prune unsampled lineages and collapse degree-2 nodes. A leaf with a taxon
     * is a sampled tip (extant or serial); a leaf without a taxon is an unsampled extinction.
     */
    private FlexibleNode reconstructTree(FlexibleNode node) {
        if (node.getChildCount() == 0) {
            return (node.getTaxon() != null) ? node : null;
        }

        List<FlexibleNode> surviving = new ArrayList<>();
        for (int i = 0; i < node.getChildCount(); i++) {
            FlexibleNode pruned = reconstructTree(node.getChild(i));
            if (pruned != null) {
                surviving.add(pruned);
            }
        }

        if (surviving.isEmpty()) {
            return null;
        }

        if (surviving.size() == 1) {
            return surviving.get(0);
        }

        FlexibleNode newNode = new FlexibleNode();
        newNode.setHeight(node.getHeight());
        for (FlexibleNode child : surviving) {
            newNode.addChild(child);
        }
        return newNode;
    }

    /**
     * An active lineage defined by its node and current age, i.e. time since the lineage's own origin
     */
    private static class LineAge {
        FlexibleNode node;
        double age;

        LineAge(FlexibleNode node, double age) {
            this.node = node;
            this.age = age;
        }
    }
}
