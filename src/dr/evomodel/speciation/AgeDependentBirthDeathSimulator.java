package dr.evomodel.speciation;

import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.Units;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Frederik M. Andersen
 *
 * Simulates symmetric or asymmetric time- and age-dependent birth-death tree using thinning algorithm
 *
 * Birth and death rates of the form
 *     lambda(t, a) = birthScale(t) * h_b(a)
 *     mu(t, a) = deathScale(t) * h_d(a)
 *     with h(a) = (1 + r * gamma * a) * exp(-gamma * a)
 */
public class AgeDependentBirthDeathSimulator {
    private final double[] birthScale;
    private final double[] deathScale;
    private final double birthB;
    private final double birthGamma;
    private final double deathB;
    private final double deathGamma;
    private final double[] epochBounds;
    private final double originTime;
    private final boolean symmetric;
    private final int maxLineages;

    // Precomputed max hazard values over [0, origin]
    private final double maxBirthHaz;
    private final double maxDeathHaz;

    /**
     * @param birthScale  piecewise-constant birth scale, one per epoch
     * @param deathScale  piecewise-constant death scale, one per epoch (or length 1 for constant)
     * @param birthShape  [r, gamma] for birth hazard, with b = r*gamma
     * @param deathShape  [r, gamma] for death hazard, with b = r*gamma
     * @param epochTimes  internal epoch boundaries in backwards time (ascending); does not include origin
     * @param originTime  the origin time (most ancient point)
     * @param symmetric   if true, both daughters get age 0; if false, one inherits parent age
     * @param maxLineages safety cap to prevent runaway growth
     */
    public AgeDependentBirthDeathSimulator(double[] birthScale,
                                           double[] deathScale,
                                           double[] birthShape,
                                           double[] deathShape,
                                           double[] epochTimes,
                                           double originTime,
                                           boolean symmetric,
                                           int maxLineages) {
        this.birthScale = birthScale;
        this.deathScale = deathScale;
        this.birthGamma = birthShape[1];
        this.birthB = birthShape[0] * this.birthGamma;
        this.deathGamma = deathShape[1];
        this.deathB = deathShape[0] * this.deathGamma;
        this.epochBounds = new double[epochTimes.length + 2];
        this.epochBounds[0] = 0.0;
        System.arraycopy(epochTimes, 0, this.epochBounds, 1, epochTimes.length);
        this.epochBounds[epochTimes.length + 1] = originTime;
        this.originTime = originTime;
        this.symmetric = symmetric;
        this.maxLineages = maxLineages;

        maxBirthHaz = maxHaz(birthB, birthGamma, originTime);
        maxDeathHaz = maxHaz(deathB, deathGamma, originTime);
    }

    /**
     * Linear-exponential hazard function: h(a) = (1 + b*a) * exp(-gamma*a)
     */
    static double ageHaz(double b, double gamma, double age) {
        return (1.0 + b * age) * Math.exp(-gamma * age);
    }

    /**
     * Maximum of h(a) over [0, maxAge]
     */
    static double maxHaz(double b, double gamma, double maxAge) {
        double hMax = Math.max(1.0, ageHaz(b, gamma, maxAge));
        if (b > 0.0 && gamma > 0.0) {
            double aStar = (b - gamma) / (gamma * b);
            if (aStar > 0.0 && aStar < maxAge) {
                hMax = Math.max(hMax, ageHaz(b, gamma, aStar));
            }
        }
        return hMax;
    }

    /**
     * Simulate a tree, retrying until the tip count is in [minTips, maxTips].
     */
    public Tree simulate(int minTips, int maxTips, int maxAttempts) {
        for (int attempt = 0; attempt < maxAttempts; attempt++) {
            Tree tree = simulateOnce();

            if (tree != null) {
                int n = tree.getExternalNodeCount();

                if (n >= minTips && n <= maxTips) {
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

    /**
     * Simulate a tree, retrying until at least minTips lineages survive to the present.
     */
    public Tree simulate(int minTips, int maxAttempts) {
        for (int attempt = 0; attempt < maxAttempts; attempt++) {
            Tree tree = simulateOnce();

            if (tree != null) {
                int n = tree.getExternalNodeCount();

                if (n >= minTips) {
                    Logger.getLogger("dr.evomodel.speciation").info(
                            "Simulated TABD tree with " + n +
                                    " tips (attempt " + (attempt + 1) + ")");
                    return tree;
                }
            }
        }
        String range =  ">= " + minTips;
        throw new RuntimeException(
                "Failed to simulate TABD tree with " + range +
                        " tips after " + maxAttempts + " attempts");
    }

    /**
     * Single simulation attempt
     */
    private Tree simulateOnce() {
        List<LineAge> lineages = new ArrayList<>();
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
            double dScale = deathScale[epoch];
            double epochBound = epochBounds[epoch];

            double maxRate = K * (bScale * maxBirthHaz + dScale * maxDeathHaz);
            if (maxRate <= 0.0) {
                double dt = currentTime - epochBound;
                for (LineAge l : lineages) l.age += dt;
                currentTime = epochBound;
                if (epoch > 0) epoch--;
                continue;
            }


            // Draw exponential(maxRate) time-to-event proposal
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

            // Pick lineage
            int idx = MathUtils.nextInt(K);
            LineAge chosen = lineages.get(idx);

            double lam = bScale * ageHaz(birthB, birthGamma, chosen.age);
            double mu = dScale * ageHaz(deathB, deathGamma, chosen.age);
            double r = lam + mu;

            // Reject
            if (MathUtils.nextDouble() >= K * r / maxRate) {
                continue;
            }

            // Accept and determine whether a speciation or extinction
            if (MathUtils.nextDouble() < lam / r) {
                FlexibleNode parent = chosen.node;
                parent.setHeight(currentTime);

                FlexibleNode child1 = new FlexibleNode();
                FlexibleNode child2 = new FlexibleNode();
                parent.addChild(child1);
                parent.addChild(child2);

                if (symmetric) {
                    // Both daughters start at age 0
                    lineages.set(idx, new LineAge(child1, 0.0));
                    lineages.add(new LineAge(child2, 0.0));
                } else {
                    // child1 inherits parent's age, child2 starts fresh
                    lineages.set(idx, new LineAge(child1, chosen.age));
                    lineages.add(new LineAge(child2, 0.0));
                }

            } else {
                FlexibleNode deadNode = chosen.node;
                deadNode.setHeight(currentTime);
                lineages.remove(idx);
            }
        }

        // Return if all is extinct
        if (lineages.isEmpty()) {
            return null;
        }

        // Set extant tips
        for (LineAge l : lineages) {
            l.node.setHeight(0.0);
            taxonCount++;
            l.node.setTaxon(new Taxon("taxon" + taxonCount));
        }

        // Prune extinct lineages from the tree
        FlexibleNode prunedRoot = reconstructTree(stemNode);

        if (prunedRoot == null || prunedRoot.getChildCount() == 0) {
            return null;
        }

        FlexibleTree tree = new FlexibleTree(prunedRoot);
        return tree;
    }

    /**
     * Recursively prune extinct lineages from the tree and collapse degree 2 nodes
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
            // Collapse: skip this node, return the single surviving child.
            // Height is preserved, so branch lengths remain correct.
            return surviving.get(0);
        }

        // Rebuild node with only surviving children
        FlexibleNode newNode = new FlexibleNode();
        newNode.setHeight(node.getHeight());
        for (FlexibleNode child : surviving) {
            newNode.addChild(child);
        }
        return newNode;
    }

    /**
     * An active lineage defined by its node and current age, i.e. time since last speciation in the lineage
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
