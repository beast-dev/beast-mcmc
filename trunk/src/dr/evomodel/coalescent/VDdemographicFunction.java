package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.tree.Tree;

import java.util.Arrays;

/**
 * @author Joseph Heled
 */
public class VDdemographicFunction extends DemographicFunction.Abstract {
    private double[] values;
    private double[] times;
    private double[] intervals;
    private double[][] ttimes;
    private double[] alltimes;
    private boolean[] dirtyTrees;
    boolean dirty;

    private final VariableDemographicModel.Type type;

    TreeIntervals[] ti;

    public VDdemographicFunction(Tree[] trees, VariableDemographicModel.Type type,
                                 double[] indicatorParameter, double[] popSizeParameter, boolean logSpace,
                                 boolean mid) {
        super(trees[0].getUnits());
        this.type = type;

        ti = new TreeIntervals[trees.length];
        dirtyTrees = new boolean[trees.length];
        Arrays.fill(dirtyTrees, true);
        ttimes = new double[ti.length][];
        int tot = 0;
        for(int k = 0; k < ti.length; ++k) {
            ttimes[k] = new double[trees[k].getTaxonCount()-1];
            tot += ttimes[k].length;
        }
        alltimes = new double[tot];

        setDirty();
        setup(trees, indicatorParameter, popSizeParameter, logSpace, mid);
    }

    /**
     * Reduce memory footprint of object. After a call to freeze only population/intensity
     * are allowed.
     */
    public void freeze() {
        ttimes = null;
        alltimes = null;
        dirtyTrees = null;
        ti = null;
    }
    
    public VDdemographicFunction(VDdemographicFunction demoFunction) {
        super(demoFunction.getUnits());
        type = demoFunction.type;

        this.ti = demoFunction.ti.clone();
        this.values = demoFunction.values.clone();
        this.times = demoFunction.times.clone();
        this.intervals = demoFunction.intervals.clone();
        this.ttimes = demoFunction.ttimes.clone();
        for(int k = 0; k < ttimes.length; ++k) {
            ttimes[k] = ttimes[k].clone();
        }

        this.alltimes = demoFunction.alltimes.clone();
        this.dirtyTrees = demoFunction.dirtyTrees.clone();
        this.dirty = demoFunction.dirty;
    }

    // Hack so that VDdemo can be used as just a linear piecewise demography (the BEAST one is broken)
    // Alexei fixed PiecewiseLinearPopulation, but did not say yet if it is tested or not.
    public VDdemographicFunction(double[] t, double[] p, Type units) {
        super(units);

// seem safe to remove them for now
//        assert t[0] >= 0;
//        for(int k = 1; k < t.length; ++k) {
//            assert t[k-1] <= t[k];
//        }
//        assert t.length + 1 == p.length;
//
        this.type = VariableDemographicModel.Type.LINEAR;
        final int tot = p.length;
        times = new double[tot+1];
        values = p;
        intervals = new double[tot -1];

        times[0] = 0.0;
        times[tot] = Double.POSITIVE_INFINITY;

//         boolean logSpace = false;
//        values[0] = logSpace ? Math.exp(p[0]) : p[0];
        System.arraycopy(t, 0, times, 1, t.length);

        for(int n = 0; n < intervals.length; ++n) {
           intervals[n] = times[n+1] - times[n];
        }
        dirty = false;
    }

    public int numberOfChanges() {
        return values.length  - 1;
    }
    
    public void treeChanged(int nt) {
        dirtyTrees[nt] = true;
        setDirty();
    }

    public void setDirty() {
        dirty = true;
    }

    private boolean setTreeTimes(int nt, Tree[] trees) {
        if( dirtyTrees[nt] ) {
            /*double[] doubles = null;
            if( ! dirtyTrees[nt] ) {
               doubles = ttimes[nt].clone();

            }*/
            ti[nt] = new TreeIntervals(trees[nt]);

            TreeIntervals nti = ti[nt];
            // make sure we get each coalescent event individually
            nti.setMultifurcationLimit(0);
            // code probably incorrect for serial samples
            final int nLineages = nti.getIntervalCount();
            assert nLineages >= ttimes[nt].length: nLineages + " " + ttimes[nt].length;

            int iCount = 0;
            for(int k = 0; k < ttimes[nt].length; ++k) {
                double timeToCoal = nti.getInterval(iCount);
                while( nti.getIntervalType(iCount) != IntervalType.COALESCENT ) {
                    ++iCount;
                    timeToCoal += nti.getInterval(iCount);
                }

                int linAtStart = nti.getLineageCount(iCount);
                ++iCount;

                assert ! (iCount == nLineages && linAtStart != 2);

                int linAtEnd = (iCount == nLineages) ? 1 : nti.getLineageCount(iCount);

                while( linAtStart <= linAtEnd ) {
                    ++iCount;
                    timeToCoal += nti.getInterval(iCount);

                    linAtStart = linAtEnd;
                    ++iCount;
                    linAtEnd = nti.getLineageCount(iCount);
                }
                ttimes[nt][k] = timeToCoal + (k == 0 ? 0 : ttimes[nt][k-1]);
            }

            /*if( doubles != null ) {
                if( ! Arrays.equals(doubles, ttimes[nt]) ) {
                   System.out.println(Arrays.toString(doubles) + " != " + Arrays.toString(ttimes[nt])
                           + Arrays.toString(dirtyTrees) + " " + dirtyTrees);
                }
            }*/
            dirtyTrees[nt] = false;
            // System.out.print(nt + " " + Arrays.toString(dirtyTrees) + " " + dirtyTrees);
            return true;
        }
        return false;
    }

    void setup(Tree[] trees, double[] indicatorParameter, double[] popSizes, boolean logSpace, boolean mid) {
        // boolean was = dirty;
        if( dirty ) {
            boolean any = false;
            for(int nt = 0; nt < ti.length; ++nt) {
                if( setTreeTimes(nt, trees) ) {
                    any = true;
                }
            }

            final int nd = indicatorParameter.length;

            assert nd == alltimes.length + (type == VariableDemographicModel.Type.STEPWISE ? -1 : 0);

            if( any ) {
                // now we want to merge times together
                int[] inds = new int[ttimes.length];

                for(int k = 0; k < alltimes.length; ++k) {
                    int j = 0;
                    while( inds[j] == ttimes[j].length ) {
                        ++j;
                    }
                    for(int l = j+1; l < inds.length; ++l) {
                        if( inds[l] < ttimes[l].length ) {
                            if( ttimes[l][inds[l]] < ttimes[j][inds[j]] ) {
                                j = l;
                            }
                        }
                    }
                    alltimes[k] = ttimes[j][inds[j]];
                    inds[j] ++;
                }
            }

            // assumes lowest node has time 0. this is probably problematic when we come
            // to deal with multiple trees

            int tot = 1;

            for(int k = 0; k < nd; ++k) {
                if( indicatorParameter[k] > 0 ) {
                    ++tot;
                }
            }

            times = new double[tot+1];
            values = new double[tot];
            intervals = new double[tot -1];

            times[0] = 0.0;
            times[tot] = Double.POSITIVE_INFINITY;

            values[0] = logSpace ? Math.exp(popSizes[0]) : popSizes[0];

            int n = 0;
            for(int k = 0; k < nd && n+1 < tot; ++k) {

                if( indicatorParameter[k] > 0 ) {
                    times[n+1] = mid ? ((alltimes[k] + (k > 0 ? alltimes[k-1] : 0))/2) : alltimes[k];

                    values[n+1] = logSpace ? Math.exp(popSizes[k+1]) : popSizes[k+1];
                    intervals[n] = times[n+1] - times[n];
                    ++n;
                }
            }
            dirty = false;
        }
        //
        /*System.out.println("after setup " + (was ? "(dirty)" : "") + " , alltimes " + Arrays.toString(alltimes)
       + " times " + Arrays.toString(times) + " values " + Arrays.toString(values) +
   " inds " + Arrays.toString(indicatorParameter.getParameterValues())) ;*/
    }

    private int getIntervalIndexStep(final double t) {
        int j = 0;
        // ugly hack,
        // when doubles are added in a different order and compared later, they can be a tiny bit off. With a
        // stepwise model this creates a "one off" situation here, which is unpleasent.
        // use float comarison here to avoid it

        final float tf = (float)t;
        while( tf > (float)times[j+1] ) ++j;
        return j;
    }

    private int getIntervalIndexLin(final double t) {
        int j = 0;
        while( t > times[j+1] ) ++j;
        return j;
    }

    public double getDemographic(double t) {

        double p;
        switch( type ) {
            case STEPWISE:
            {
                final int j = getIntervalIndexStep(t);
                p = values[j];
                break;
            }
            case LINEAR:
            {
                final int j = getIntervalIndexLin(t);
                if( j == values.length - 1 ) {
                    p = values[j];
                    break;
                }

                final double a = (t - times[j]) / (intervals[j]);
                p = a * values[j+1] + (1-a) * values[j];
                break;
            }
            default: throw new IllegalArgumentException("");

        }
        return p;
    }

    public double getIntensity(double t) {
        return getIntegral(0, t);
    }

    public double getInverseIntensity(double x) {
        assert false;
        return 0;
    }

    private double intensityLinInterval(double start, double end, int index) {
        final double dx = end - start;

        final double popStart = values[index];
        final double popDiff = (index < values.length - 1 ) ? values[index + 1] - popStart : 0.0;
        if( popDiff == 0.0 ) {
            return dx / popStart;
        }
        final double time0 = times[index];
        final double interval = intervals[index];

        assert (float)start <= (float)(time0 + interval) && start >= time0 && (float)end <= (float)(time0 + interval) && end >= time0;

//        final double pop0 = popStart + ((start - time0) / interval) * popDiff;
//       final double pop1 = popStart + ((end - time0) / interval) * popDiff;

        // do same as above more effeciently
//        final double r = popDiff / interval;
//        final double x = popStart - time0 * r;
//        final double pop0 = x + start * r;
//        final double pop1 = x + end * r;
        //better numerical stability but not perfect
        double p1minusp0 = ((end-start)/interval) * popDiff;

        double v = interval * (popStart / popDiff);
        final double p1overp0 = (v + (end- time0)) / (v + (start- time0));
        if( p1minusp0 == 0.0 || p1overp0 <= 0 ) {
            // either dx == 0 or is very small (numerical inaccuracy)
            final double pop0 = popStart + ((start - time0) / interval) * popDiff;
            return dx/pop0;
        }

        return dx * Math.log(p1overp0) / p1minusp0;
       // return dx * Math.log(pop1/pop0) / (pop1 - pop0);*/
    }

    private double intensityLinInterval(int index) {
        final double interval = intervals[index];
        final double pop0 = values[index];
        final double pop1 =  values[index+1];
        if( pop0 == pop1 ) {
            return interval / pop0;
        }
        return interval * Math.log(pop1/pop0) / (pop1 - pop0);
    }

    public double getIntegral(double start, double finish) {

        double intensity = 0.0;

        switch( type ) {
            case STEPWISE:
            {
                final int first = getIntervalIndexStep(start);
                final int last = getIntervalIndexStep(finish);

                final double popStart = values[first];
                if( first == last ) {
                    intensity = (finish - start) / popStart;
                } else {
                    intensity = (times[first+1] - start) / popStart;

                    for(int k = first + 1; k < last; ++k) {
                        intensity += intervals[k]/values[k];
                    }
                    intensity += (finish - times[last]) / values[last];
                }
                break;
            }
            case LINEAR:
            {
                final int first = getIntervalIndexLin(start);
                final int last = getIntervalIndexLin(finish);

                if( first == last ) {
                    intensity += intensityLinInterval(start, finish, first);
                } else {
                    // from first to end of interval
                    intensity += intensityLinInterval(start, times[first+1], first);
                    // intervals until (not including) last
                    for(int k = first + 1; k < last; ++k) {
                        intensity += intensityLinInterval(k);
                    }
                    // last interval
                    intensity += intensityLinInterval(times[last], finish, last);
                }
                break;
            }
        }
        return intensity;
    }

    public int getNumArguments() {
        assert false;
        return 0;
    }

    public String getArgumentName(int n) {
        assert false;
        return null;
    }

    public double getArgument(int n) {
        assert false;
        return 0;
    }

    public void setArgument(int n, double value) {
        assert false;
    }

    public double getLowerBound(int n) {
        return 0.0;
    }

    public double getUpperBound(int n) {
        return Double.POSITIVE_INFINITY;
    }

    public DemographicFunction getCopy() {
        return null;
    }

    // not sure why we need this here
    public double value(double x) {
        return 1.0 / getDemographic(x);
    }

    public TreeIntervals getTreeIntervals(int nt) {
        return ti[nt];
    }

    public double[] allTimePoints() {
        return alltimes;
    }

    public String toString() {
        final StringBuilder sb = new StringBuilder(32);

        for(int k = 1; k < times.length-1; ++k) {
            if( k > 1 ) {
                sb.append(",");
            }
            sb.append(times[k]);
        }
        sb.append("|");
        sb.append(values[0]);
        for(int k = 1; k < values.length; ++k) {

            sb.append(",");
            sb.append(values [k]);
        }
        return sb.toString();
    }
}

