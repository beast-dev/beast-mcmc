package dr.evomodel.speciation;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.coalescent.VDdemographicFunction;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.util.Attributable;
import dr.util.HeapSort;
import dr.xml.*;
import jebl.util.FixedBitSet;

import java.util.*;

/**
 * Species tree which includes demographic function per branch.
 *
 * @author joseph
 *         Date: 24/05/2008
 */
public class SpeciesTreeModel extends AbstractModel implements MutableTree, NodeAttributeProvider {
    public static final String SPECIES_TREE = "speciesTree";

    private static final String SPP_SPLIT_POPULATIONS = "sppSplitPopulations";
    private static final String COALESCENT_POINTS_POPULATIONS = "coalescentPointsPopulations";
    private static final String COALESCENT_POINTS_INDICATORS = "coalescentPointsIndicators";

    private final SimpleTree spTree;
    private SpeciesBindings species;
    Map<NodeRef, NodeProperties> props = new HashMap<NodeRef, NodeProperties>();
    public final Parameter sppSplitPopulations;
    private int[] singleStartPoints;
    private int[] pairStartPoints;

    private Parameter coalPointsPops;
    private Parameter coalPointsIndicator;
    private boolean nodePropsReady;

    private NodeRef[] children;
    private double[] heights;

    // any change of underlying parameters / models
    private boolean anyChange;
    // Tree has been edited in this cycle
    private boolean treeChanged;

    private final String spIndexAttrName = "spi";

    private class NodeProperties {
        final int speciesIndex;
        public VDdemographicFunction demogf;
        FixedBitSet spSet;

        public NodeProperties(int n) {
            speciesIndex = n;
            demogf = null;
            spSet = new FixedBitSet(species.nSpecies());
        }
    }

    SpeciesTreeModel(SpeciesBindings species, Parameter sppSplitPopulations,
                     Parameter coalPointsPops, Parameter coalPointsIndicator,
                     Tree startTree) {
        super(SPECIES_TREE);

        this.species = species;

        this.sppSplitPopulations = sppSplitPopulations;
        this.coalPointsPops = coalPointsPops;
        this.coalPointsIndicator = coalPointsIndicator;

        addParameter(sppSplitPopulations);

        addModel(species);

        if (coalPointsPops != null) {
            assert coalPointsIndicator != null;

            addParameter(coalPointsPops);
            addParameter(coalPointsIndicator);

            final double[][] pts = species.getPopTimesSingle();

            int start = 0;
            singleStartPoints = new int[pts.length];
            for (int i = 0; i < pts.length; i++) {
                singleStartPoints[i] = start;
                start += pts[i].length;
            }

            final double[][] ptp = species.getPopTimesPair();
            pairStartPoints = new int[ptp.length];
            for (int i = 0; i < ptp.length; i++) {
                pairStartPoints[i] = start;
                start += ptp[i].length;
            }
        }

        // build an initial noninformative tree
        spTree = compatibleUninformedSpeciesTree(startTree);

        // some of the code is generic but some parts assume a binary tree.
        assert Tree.Utils.isBinary(spTree);

        final int nNodes = spTree.getNodeCount();
        heights = new double[nNodes];
        children = new NodeRef[2 * nNodes + 1];

        // fixed properties
        for (int k = 0; k < getExternalNodeCount(); ++k) {
            final NodeRef nodeRef = getExternalNode(k);
            int n = (Integer) getNodeAttribute(nodeRef, spIndexAttrName);
            final NodeProperties np = new NodeProperties(n);
            props.put(nodeRef, np);
            np.spSet.set(n);
        }

        for (int k = 0; k < getInternalNodeCount(); ++k) {
            final NodeRef nodeRef = getInternalNode(k);
            props.put(nodeRef, new NodeProperties(-1));
        }

        nodePropsReady = false;
    }


    public boolean isCompatible(SpeciesBindings.GeneTreeInfo geneTreeInfo) {
        // can't set demographics if a tree is not compatible, but we need spSets.
        if (!nodePropsReady) {
            setSPsets(getRoot());
        }
        return isCompatible(getRoot(), geneTreeInfo.getCoalInfo(), 0) >= 0;
    }

    // Not very effecient, should do something better, based on traversing the cList once

    private int isCompatible(NodeRef node, SpeciesBindings.CoalInfo[] cList, int loc) {
        if (!isExternal(node)) {
            int l = -1;
            for (int nc = 0; nc < getChildCount(node); ++nc) {
                int l1 = isCompatible(getChild(node, nc), cList, loc);
                if (l1 < 0) {
                    return -1;
                }
                assert l == -1 || l1 == l;

                l = l1;
            }
            loc = l;

            assert cList[loc].ctime >= getNodeHeight(node);
        }

        if (node == getRoot()) {
            return cList.length;
        }

        // spSet guaranteed to be ready by caller
        final FixedBitSet nodeSps = props.get(node).spSet;

        final double limit = getNodeHeight(getParent(node));

        while (loc < cList.length) {
            final SpeciesBindings.CoalInfo ci = cList[loc];
            if (ci.ctime >= limit) {
                break;
            }
            boolean allIn = true, noneIn = true;

            for (int i = 0; i < 2; ++i) {
                final FixedBitSet s = ci.sinfo[i];

                final int in1 = s.intersectCardinality(nodeSps);
                if (in1 > 0) {
                    noneIn = false;
                }
                if (s.cardinality() != in1) {
                    allIn = false;
                }
            }
            if (!(allIn || noneIn)) {
                return -1;
            }
            ++loc;
        }
        return loc;
    }


    private static double fp(double val, double low, double[][] tt, int[] ii) {
        for (int k = 0; k < ii.length; ++k) {
            int ip = ii[k];
            if (ip == tt[k].length || val <= tt[k][ip]) {
                --ip;
                while (ip >= 0 && val <= tt[k][ip]) {
                    --ip;
                }
                assert ((ip < 0) || (tt[k][ip] < val)) && ((ip + 1 == tt[k].length) || (val <= tt[k][ip + 1]));
                if (ip >= 0) {
                    low = Math.max(low, tt[k][ip]);
                }
            } else {
                ++ip;
                while (ip < tt[k].length && val > tt[k][ip]) {
                    ++ip;
                }
                assert tt[k][ip - 1] < val && ((ip == tt[k].length) || (val <= tt[k][ip]));
                low = Math.max(low, tt[k][ip - 1]);
            }
        }
        return low;
    }

    // Pass arguments of recursive functions in a compact format.

    private class Args {
        final double[][] cps = species.getPopTimesSingle();
        final double[][] cpp = species.getPopTimesPair();
        final int[] iSingle = new int[cps.length];
        final int[] iPair = new int[cpp.length];

        final double[] indicators = ((Parameter.Default) coalPointsIndicator).inspectParameterValues();
        final double[] pops = ((Parameter.Default) coalPointsPops).inspectParameterValues();

        private double findPrev(double val, double low) {
            low = fp(val, low, cps, iSingle);
            low = fp(val, low, cpp, iPair);

            return low;
        }
    }


    static private class Points implements Comparable<Points> {
        final double time;
        final double population;

        Points(double t, double p) {
            time = t;
            population = p;
        }

        public int compareTo(Points points) {
            return time < points.time ? -1 : (time > points.time ? 1 : 0);
        }
    }

    private NodeProperties setSPsets(NodeRef nodeID) {
        final NodeProperties nprop = props.get(nodeID);

        if (!isExternal(nodeID)) {
            nprop.spSet = new FixedBitSet(species.nSpecies());
            for (int nc = 0; nc < getChildCount(nodeID); ++nc) {
                NodeProperties p = setSPsets(getChild(nodeID, nc));
                nprop.spSet.union(p.spSet);
            }
        }
        return nprop;
    }

    //  Assign positions in 'pointsList' for the sub-tree rooted at the ancestor of
    //  nodeID.
    //
    //  pointsList is indexed by node-id. Every element is a list of internal
    //  population points for the branch between nodeID and it's ancestor
    //
    private NodeProperties getDemographicPoints(NodeRef nodeID, Args args, Points[][] pointsList) {

        final NodeProperties nprop = props.get(nodeID);
        final int nSpecies = species.nSpecies();

        if (isExternal(nodeID)) {
            //spset = frozenset((data.speciesSet.argmax(),))
        } else {
            nprop.spSet = new FixedBitSet(nSpecies);
            for (int nc = 0; nc < getChildCount(nodeID); ++nc) {
                NodeProperties p = getDemographicPoints(getChild(nodeID, nc), args, pointsList);
                nprop.spSet.union(p.spSet);
            }
        }

        final double cheight = nodeID != getRoot() ? getNodeHeight(getParent(nodeID)) : Double.MAX_VALUE;

        List<Points> allPoints = new ArrayList<Points>(5);

        for (int x = nprop.spSet.nextOnBit(0); x >= 0; x = nprop.spSet.nextOnBit(x + 1)) {
            final double nodeHeight = spTree.getNodeHeight(nodeID);
            {
                double[] cp = args.cps[x];
                final int upi = singleStartPoints[x];

                int i = args.iSingle[x];

                while (i < cp.length && cp[i] < cheight) {
                    if (args.indicators[upi + i] > 0) {
                        //System.out.println("  popbit s");
                        args.iSingle[x] = i;
                        double prev = args.findPrev(cp[i], nodeHeight);
                        double mid = (prev + cp[i]) / 2.0;
                        assert nodeHeight < mid;
                        allPoints.add(new Points(mid, args.pops[upi + i]));
                    }
                    ++i;
                }
                args.iSingle[x] = i;
            }

            final int kx = (x * (2 * nSpecies - x - 3)) / 2 - 1;
            for (int y = nprop.spSet.nextOnBit(x + 1); y >= 0; y = nprop.spSet.nextOnBit(y + 1)) {

                assert x < y;
                int k = kx + y;

                double[] cp = args.cpp[k];
                int i = args.iPair[k];
                final int upi = pairStartPoints[k];

                while (i < cp.length && cp[i] < cheight) {
                    if (args.indicators[upi + i] > 0) {
                        //System.out.println("  popbit p");
                        args.iPair[k] = i;
                        final double prev = args.findPrev(cp[i], nodeHeight);
                        double mid = (prev + cp[i]) / 2.0;
                        assert nodeHeight < mid;
                        allPoints.add(new Points(mid, args.pops[upi + i]));
                    }
                    ++i;
                }
                args.iPair[k] = i;
            }
        }
        Points[] all = null;

        if (allPoints.size() > 0) {
            all = allPoints.toArray(new Points[allPoints.size()]);
            if (all.length > 1) {
                HeapSort.sort(all);
            }

            // duplications
            int len = all.length;
            int k = 0;
            while (k + 1 < len) {
                double t = all[k].time;
                if (t == all[k + 1].time) {
                    int j = k + 2;
                    double v = all[k].population + all[k + 1].population;
                    while (j < len && t == all[j].time) {
                        v += all[j].population;
                        j += 1;
                    }
                    int removed = (j - k - 1);
                    all[k] = new Points(t, v / (removed + 1));
                    for (int i = k + 1; i < len - removed; ++i) {
                        all[i] = all[i + removed];
                    }
                    //System.arraycopy(all, j, all, k + 1, all.length - j + 1);
                    len -= removed;
                }
                ++k;
            }
            if (len != all.length) {
                Points[] a = new Points[len];
                System.arraycopy(all, 0, a, 0, len);
                all = a;
            }
        }
        pointsList[nodeID.getNumber()] = all;
        return nprop;
    }

    private int setDemographics(NodeRef nodeID, int pStart, int side, double[] pops, Points[][] pointsList) {

        final int nSpecies = species.nSpecies();
        final NodeProperties nprop = props.get(nodeID);
        int pEnd;
        double p0;

        if (isExternal(nodeID)) {
            final int sps = nprop.speciesIndex;
            p0 = pops[sps];
            pEnd = pStart;
        } else {
            assert getChildCount(nodeID) == 2;

            final int iHere = setDemographics(getChild(nodeID, 0), pStart, 0, pops, pointsList);
            pEnd = setDemographics(getChild(nodeID, 1), iHere + 1, 1, pops, pointsList);
            final int i = nSpecies + iHere * 2;
            p0 = pops[i] + pops[i + 1];
        }

        final double t0 = getNodeHeight(nodeID);

        Points[] p = pointsList[nodeID.getNumber()];
        final int plen = p == null ? 0 : p.length;

        final boolean isRoot = nodeID == getRoot();
//        double[] xs = new double[plen + (isRoot ? 1 : 1)];
//        double[] ys = new double[plen + (isRoot ? 2 : 2)];

        double[] xs = new double[plen + 1];
        double[] ys = new double[plen + 2];
        ys[0] = p0;
        for (int i = 0; i < plen; ++i) {
            xs[i] = p[i].time - t0;
            ys[i + 1] = p[i].population;
        }

        if (!isRoot) {
            final int anccIndex = (side == 0) ? pEnd : pStart - 1;
            final double pe = pops[nSpecies + anccIndex * 2 + side];
            final double b = getBranchLength(nodeID);

            xs[xs.length - 1] = b;
            ys[ys.length - 1] = pe;
        } else {
            // extend the last point to most ancient coalescent point. Has no effect on the demographic
            // per se but for use when analyzing the results.

            double h = -1;
            for (SpeciesBindings.GeneTreeInfo t : species.getGeneTrees()) {
                h = Math.max(h, t.tree.getNodeHeight(t.tree.getRoot()));
            }
            xs[xs.length - 1] = h - getNodeHeight(nodeID);
            ys[ys.length - 1] = ys[ys.length - 2];
        }

        nprop.demogf = new VDdemographicFunction(xs, ys, getUnits());
        return pEnd;
    }

    private void setNodeProperties() {
        Args args = new Args();
        Points[][] perBranchPoints = new Points[getNodeCount()][];
        getDemographicPoints(getRoot(), args, perBranchPoints);

        setDemographics(getRoot(), 0, -1, ((Parameter.Default) sppSplitPopulations).inspectParameterValues(), perBranchPoints);
    }

    private Map<NodeRef, NodeProperties> getProps() {
        if (!nodePropsReady) {
            setNodeProperties();
            nodePropsReady = true;
        }
        return props;
    }

    public DemographicFunction getNodeDemographic(NodeRef node) {
        return getProps().get(node).demogf;
    }

    public FixedBitSet spSet(NodeRef node) {
        return getProps().get(node).spSet;
    }

    public int speciesIndex(NodeRef tip) {
        assert isExternal(tip);

        // always ready even if props is dirty
        return props.get(tip).speciesIndex;
    }

    private SimpleTree compatibleUninformedSpeciesTree(Tree startTree) {
        double rootHeight = Double.MAX_VALUE;

        for (SpeciesBindings.GeneTreeInfo t : species.getGeneTrees()) {
            rootHeight = Math.min(rootHeight, t.getCoalInfo()[0].ctime);
        }

        final SpeciesBindings.SPinfo[] spp = species.species;

        if (startTree != null) {
            // Allow start tree to be very basic basic - may be not fully resolved and no
            // branch lengths

            if (startTree.getExternalNodeCount() != spp.length) {
                throw new Error("Start tree error - different number of tips");
            }

            final FlexibleTree tree = new FlexibleTree(startTree);
            tree.resolveTree();
            if (tree.getRootHeight() <= 0.0) {
                tree.setRootHeight(1.0);
            }
            Utils.correctHeightsForTips(tree);
            SimpleTree.Utils.scaleNodeHeights(tree, rootHeight / tree.getRootHeight());

            SimpleTree sTree = new SimpleTree(tree);
            for (int ns = 0; ns < spp.length; ns++) {
                SpeciesBindings.SPinfo sp = spp[ns];
                final int i = sTree.getTaxonIndex(sp.name);
                if (i < 0) {
                    throw new Error(sp.name + " is not present in the start tree");
                }
                final SimpleNode node = sTree.getExternalNode(i);
                node.setAttribute(spIndexAttrName, ns);
            }
            return sTree;
        }


        final double delta = rootHeight / (spp.length + 1);
        double cTime = delta;

        List<SimpleNode> subs = new ArrayList<SimpleNode>(spp.length);

        for (int ns = 0; ns < spp.length; ns++) {
            SpeciesBindings.SPinfo sp = spp[ns];
            final SimpleNode node = new SimpleNode();
            node.setTaxon(new Taxon(sp.name));
            subs.add(node);

            node.setAttribute(spIndexAttrName, ns);
        }

        while (subs.size() > 1) {
            final SimpleNode node = new SimpleNode();
            int i = 0, j = 1;
            node.addChild(subs.get(i));
            node.addChild(subs.get(j));
            node.setHeight(cTime);
            cTime += delta;
            subs.set(j, node);
            subs.remove(i);
        }

        return new SimpleTree(subs.get(0));
    }

    private boolean verbose = false;

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (verbose) System.out.println(" SPtree: model changed " + model.getId());

        nodePropsReady = false;
        anyChange = true;
        // this should happen by default, no?
        fireModelChanged();
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        if (verbose) System.out.println(" SPtree: parameter changed " + parameter.getId());

        nodePropsReady = false;
        anyChange = true;
    }

    protected void storeState() {
        assert !treeChanged;
        assert !anyChange;
    }

    protected void restoreState() {
        if (verbose) System.out.println(" SPtree: restore (" + treeChanged + "," + anyChange + ")");

        if (treeChanged) {
            //
            spTree.beginTreeEdit();

            for (int k = 0; k < getInternalNodeCount(); ++k) {
                final NodeRef node = getInternalNode(k);
                final int index = node.getNumber();
                final double h = heights[index];
                if (getNodeHeight(node) != h) {
                    setNodeHeight(node, h);
                }
                for (int nc = 0; nc < 2; ++nc) {
                    final NodeRef child = getChild(node, nc);

                    final NodeRef child1 = children[2 * index + nc];
                    if (child != child1) {
                        replaceChild(node, child, child1);
                    }
                    assert getParent(child1) == node;
                }
            }
            setRoot(children[children.length - 1]);

            if (verbose) System.out.println("  restored to: " + spTree);

            spTree.endTreeEdit();
        }
        if (treeChanged || anyChange) {
            setNodeProperties();
        }
        treeChanged = false;
        anyChange = false;
    }

    protected void acceptState() {
        if (verbose) System.out.println(" SPtree: accept");

        treeChanged = false;
        anyChange = false;
    }


    public String[] getNodeAttributeLabel() {
        // keep short, repeated endlessly in tree log
        return new String[]{"dmf"};
    }

    public String[] getAttributeForNode(Tree tree, NodeRef node) {
        assert tree == this;

        final VDdemographicFunction df = getProps().get(node).demogf;

        return new String[]{"{" + df.toString() + "}"};
    }

    // boring delagation
    public SimpleTree getSimpleTree() {
        return spTree;
    }

    public Tree getCopy() {
        return spTree.getCopy();
    }

    public Type getUnits() {
        return spTree.getUnits();
    }

    public void setUnits(Type units) {
        spTree.setUnits(units);
    }

    public int getNodeCount() {
        return spTree.getNodeCount();
    }

    public boolean hasNodeHeights() {
        return spTree.hasNodeHeights();
    }

    public double getNodeHeight(NodeRef node) {
        return spTree.getNodeHeight(node);
    }

    public double getNodeRate(NodeRef node) {
        return spTree.getNodeRate(node);
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return spTree.getNodeTaxon(node);
    }

    public int getChildCount(NodeRef node) {
        return spTree.getChildCount(node);
    }

    public boolean isExternal(NodeRef node) {
        return spTree.isExternal(node);
    }

    public boolean isRoot(NodeRef node) {
        return spTree.isRoot(node);
    }

    public NodeRef getChild(NodeRef node, int i) {
        return spTree.getChild(node, i);
    }

    public NodeRef getParent(NodeRef node) {
        return spTree.getParent(node);
    }

    public boolean hasBranchLengths() {
        return spTree.hasBranchLengths();
    }

    public double getBranchLength(NodeRef node) {
        return spTree.getBranchLength(node);
    }

    public void setBranchLength(NodeRef node, double length) {
        spTree.setBranchLength(node, length);
    }

    public NodeRef getExternalNode(int i) {
        return spTree.getExternalNode(i);
    }

    public NodeRef getInternalNode(int i) {
        return spTree.getInternalNode(i);
    }

    public NodeRef getNode(int i) {
        return spTree.getNode(i);
    }

    public int getExternalNodeCount() {
        return spTree.getExternalNodeCount();
    }

    public int getInternalNodeCount() {
        return spTree.getInternalNodeCount();
    }

    public NodeRef getRoot() {
        return spTree.getRoot();
    }

    public void setRoot(NodeRef r) {
        spTree.setRoot(r);
    }

    public void addChild(NodeRef p, NodeRef c) {
        spTree.addChild(p, c);
    }

    public void removeChild(NodeRef p, NodeRef c) {
        spTree.removeChild(p, c);
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        spTree.replaceChild(node, child, newChild);
    }

    public boolean beginTreeEdit() {
        boolean beingEdited = spTree.beginTreeEdit();
        if (!beingEdited) {
            // save tree for restore
            for (int n = 0; n < getInternalNodeCount(); ++n) {
                final NodeRef node = getInternalNode(n);
                final int k = node.getNumber();
                children[2 * k] = getChild(node, 0);
                children[2 * k + 1] = getChild(node, 1);
                heights[k] = getNodeHeight(node);
            }
            children[children.length - 1] = getRoot();

            treeChanged = true;
            nodePropsReady = false;
            //anyChange = true;
        }
        return beingEdited;
    }

    public void endTreeEdit() {
        spTree.endTreeEdit();
        fireModelChanged();
    }

    public void setNodeHeight(NodeRef n, double height) {
        spTree.setNodeHeight(n, height);
    }

    public void setNodeRate(NodeRef n, double rate) {
        spTree.setNodeRate(n, rate);
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        spTree.setNodeAttribute(node, name, value);
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        return spTree.getNodeAttribute(node, name);
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        return spTree.getNodeAttributeNames(node);
    }

    public int getTaxonCount() {
        return spTree.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex) {
        return spTree.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex) {
        return spTree.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id) {
        return spTree.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon) {
        return spTree.getTaxonIndex(taxon);
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        return spTree.getTaxonAttribute(taxonIndex, name);
    }

    public int addTaxon(Taxon taxon) {
        return spTree.addTaxon(taxon);
    }

    public boolean removeTaxon(Taxon taxon) {
        return spTree.removeTaxon(taxon);
    }

    public void setTaxonId(int taxonIndex, String id) {
        spTree.setTaxonId(taxonIndex, id);
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        spTree.setTaxonAttribute(taxonIndex, name, value);
    }

    public String getId() {
        return spTree.getId();
    }

    public void setId(String id) {
        spTree.setId(id);
    }

    public void setAttribute(String name, Object value) {
        spTree.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        return spTree.getAttribute(name);
    }

    public Iterator getAttributeNames() {
        return spTree.getAttributeNames();
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        spTree.addMutableTreeListener(listener);
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        spTree.addMutableTaxonListListener(listener);
    }

    private static Parameter createCoalPointsPopParameter(SpeciesBindings spb, Double value) {
        int dim = 0;
        for (double[] d : spb.getPopTimesSingle()) {
            dim += d.length;
        }
        for (double[] d : spb.getPopTimesPair()) {
            dim += d.length;
        }

        return new Parameter.Default(dim, value);
    }

    private static Parameter createSplitPopulationsParameter(SpeciesBindings spb, double value) {
        final int dim = 3 * spb.nSpecies() - 2;
        return new Parameter.Default(dim, value);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            SpeciesBindings spb = (SpeciesBindings) xo.getChild(SpeciesBindings.class);
            XMLObject cxo = (XMLObject) xo.getChild(SPP_SPLIT_POPULATIONS);

            double value = cxo.getAttribute(Attributable.VALUE, 1.0);
            Parameter sppSplitPopulations =
                    createSplitPopulationsParameter(spb, cxo.getAttribute(Attributable.VALUE, value));
            replaceParameter(cxo, sppSplitPopulations);
            sppSplitPopulations.addBounds(
                    new Parameter.DefaultBounds(Double.MAX_VALUE, 0, sppSplitPopulations.getDimension()));

            cxo = (XMLObject) xo.getChild(COALESCENT_POINTS_POPULATIONS);
            value = cxo.getAttribute(Attributable.VALUE, 1.0);
            Parameter coalPointsPops = createCoalPointsPopParameter(spb, cxo.getAttribute(Attributable.VALUE, value));
            replaceParameter(cxo, coalPointsPops);
            coalPointsPops.addBounds(
                    new Parameter.DefaultBounds(Double.MAX_VALUE, 0, coalPointsPops.getDimension()));

            cxo = (XMLObject) xo.getChild(COALESCENT_POINTS_INDICATORS);
            Parameter coalPointsIndicators = new Parameter.Default(coalPointsPops.getDimension(), 0);
            replaceParameter(cxo, coalPointsIndicators);

            Tree startTree = (Tree) xo.getChild(Tree.class);

            return new SpeciesTreeModel(spb, sppSplitPopulations, coalPointsPops, coalPointsIndicators, startTree);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(SpeciesBindings.class),
                    // A starting tree. Can be very minimal, i.e. no branch lengths and not resolved
                    new ElementRule(Tree.class, true),
                    new ElementRule(SPP_SPLIT_POPULATIONS, new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(Attributable.VALUE, true),
                            new ElementRule(Parameter.class)}),

                    new ElementRule(COALESCENT_POINTS_POPULATIONS, new XMLSyntaxRule[]{
                            AttributeRule.newDoubleRule(Attributable.VALUE, true),
                            new ElementRule(Parameter.class)}, true),

                    new ElementRule(COALESCENT_POINTS_INDICATORS, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)}, true),
            };
        }

        public String getParserDescription() {
            return "Species tree which includes demographic function per branch.";
        }

        public Class getReturnType() {
            return SpeciesTreeModel.class;
        }

        public String getParserName() {
            return SPECIES_TREE;
        }
    };
}
