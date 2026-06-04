package dr.evomodel.speciation.agedependent.simulation;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Frederik M. Andersen
 *
 * Serially-sampled variant of {@link AgeDependentBirthDeathSimulator}. Adds a per-time serial
 * sampling rate psi(t) (skyline) and an extant sampling probability rho. Sampled-through-time
 * tips are produced by psi events during the simulation (lineage is removed on sampling);
 * surviving lineages at the present are each sampled with probability rho.
 *
 * Rates:
 *     lambda(t, a) = birthScale(t) * h_b(a)
 *     mu(t, a)     = deathScale(t) * h_d(a)
 *     psi(t)       = samplingScale(t)
 *     h(a) = (1 + r * gamma * a) * exp(-gamma * a)
 */
public class AgeDependentBirthDeathSerialSimulator {
    private final double[] birthScale;
    private final double[] deathScale;
    private final double[] samplingScale;
    private final double extantSamplingProb;
    private final double[] birthB;
    private final double[] birthGamma;
    private final double[] deathB;
    private final double[] deathGamma;
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
     * @param birthHazard         [r, gamma] for birth hazard, with b = r*gamma; length 2 (shared) or 2*numEpochs
     * @param deathHazard         [r, gamma] for death hazard, with b = r*gamma; length 2 (shared) or 2*numEpochs
     * @param epochTimes         internal epoch boundaries in backwards time (ascending); excludes origin
     * @param originTime         the origin time (most ancient point)
     * @param symmetric          if true, both daughters get age 0; if false, one inherits parent age
     * @param maxLineages        safety cap to prevent runaway growth
     */
    public AgeDependentBirthDeathSerialSimulator(double[] birthScale,
                                                 double[] deathScale,
                                                 double[] samplingScale,
                                                 double extantSamplingProb,
                                                 double[] birthHazard,
                                                 double[] deathHazard,
                                                 double[] epochTimes,
                                                 double originTime,
                                                 boolean symmetric,
                                                 int maxLineages) {
        this.birthScale = birthScale;
        this.deathScale = deathScale;
        this.samplingScale = samplingScale;
        this.extantSamplingProb = extantSamplingProb;
        int numEpochs = epochTimes.length + 1;
        this.birthGamma = new double[numEpochs];
        this.birthB = new double[numEpochs];
        this.deathGamma = new double[numEpochs];
        this.deathB = new double[numEpochs];
        boolean constBirthShape = birthHazard.length == 2;
        boolean constDeathShape = deathHazard.length == 2;
        for (int k = 0; k < numEpochs; k++) {
            int bOff = constBirthShape ? 0 : 2 * k;
            int dOff = constDeathShape ? 0 : 2 * k;
            this.birthGamma[k] = birthHazard[bOff + 1];
            this.birthB[k]     = birthHazard[bOff] * this.birthGamma[k];
            this.deathGamma[k] = deathHazard[dOff + 1];
            this.deathB[k]     = deathHazard[dOff] * this.deathGamma[k];
        }
        this.epochBounds = new double[epochTimes.length + 2];
        this.epochBounds[0] = 0.0;
        System.arraycopy(epochTimes, 0, this.epochBounds, 1, epochTimes.length);
        this.epochBounds[epochTimes.length + 1] = originTime;
        this.originTime = originTime;
        this.symmetric = symmetric;
        this.maxLineages = maxLineages;

        this.maxBirthHaz = new double[numEpochs];
        this.maxDeathHaz = new double[numEpochs];
        for (int k = 0; k < numEpochs; k++) {
            maxBirthHaz[k] = AgeDependentBirthDeathSimulator.maxHaz(birthB[k], birthGamma[k], originTime);
            maxDeathHaz[k] = AgeDependentBirthDeathSimulator.maxHaz(deathB[k], deathGamma[k], originTime);
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
                            "Simulated serial TABD tree with " + n +
                                    " tips (attempt " + (attempt + 1) + ")");
                    return tree;
                }
            }
        }
        String range = maxTips > 0 ? "[" + minTips + ", " + maxTips + "]" : ">= " + minTips;
        throw new RuntimeException(
                "Failed to simulate serial TABD tree with " + range +
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

            double lam = bScale * AgeDependentBirthDeathSimulator.ageHaz(birthB[epoch], birthGamma[epoch], chosen.age);
            double mu  = dScale * AgeDependentBirthDeathSimulator.ageHaz(deathB[epoch], deathGamma[epoch], chosen.age);
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

        return new FlexibleTree(prunedRoot);
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

    private static class LineAge {
        FlexibleNode node;
        double age;

        LineAge(FlexibleNode node, double age) {
            this.node = node;
            this.age = age;
        }
    }
}
