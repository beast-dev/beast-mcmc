package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;
import java.util.Arrays;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class VariableDemographicModel extends DemographicModel {
    public static final String MODEL_NAME = "variableDemographic";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String POPULATION_TREE = "populationTree";
    public static final String LOG_SPACE = "logUnits";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    /*private*/ Parameter popSizeParameter;
    /*private*/ Parameter indicatorParameter;
    private Type type;
    private boolean logSpace;
    private TreeModel[] trees;
    private VD demoFunction = null;
    private VD savedDemoFunction = null;

    public Parameter getIndices() {
        return indicatorParameter;
    }

    public enum Type {
        STEPWISE,
        LINEAR,
        EXPONENTIAL
    }

    public VariableDemographicModel(TreeModel[] trees, Parameter popSizeParameter, Parameter indicatorParameter,
                                    Type type, boolean logSpace) {
        super(MODEL_NAME);

        this.popSizeParameter = popSizeParameter;
        this.indicatorParameter = indicatorParameter;
        final int redcueDim = type == Type.STEPWISE ? 1 : 0;
        int events = 0;
        for( Tree t : trees ) {
            events += t.getExternalNodeCount() - redcueDim;
            // we will have to handle this I guess
            assert t.getUnits() == trees[0].getUnits();
        }
        //final int events = tree.getExternalNodeCount() - redcueDim;
        final int paramDim1 = popSizeParameter.getDimension();
        final int paramDim2 = indicatorParameter.getDimension();
        this.type = type;
        this.logSpace = logSpace;
        if( logSpace ) {
            throw new IllegalArgumentException("sorry log space not implemented");
        }

        if (paramDim1 != events) {
            throw new IllegalArgumentException("Dimension of population parameter must be the same as the number of internal nodes in the tree. ("
            + paramDim1 + " != " + events + ")");
        }
        
        if (paramDim2 != events - 1) {
            throw new IllegalArgumentException("Dimension of indicator parameter must one less than the number of internal nodes in the tree. ("
            + paramDim2 + " != " + (events-1) + ")");
        }

        this.trees = trees;

        for( TreeModel t : trees ) {
          addModel(t);
        }

        addParameter(indicatorParameter);
        addParameter(popSizeParameter);
    }

    public int nLoci() {
        return trees.length;
    }
    
    public Tree getTree(int k) {
        return trees[k];
    }

    public VD getDemographicFunction() {
        if( demoFunction == null ) {
            demoFunction = new VD();
        } else {
            demoFunction.setup();
        }
        return demoFunction;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // tree has changed

        if( demoFunction != null ) {
            if( demoFunction == savedDemoFunction ) {
                demoFunction = new VD(demoFunction);
            }
            for(int k = 0; k < trees.length; ++k) {
                if( model == trees[k] ) {
                    demoFunction.treeChanged(k);
                    break;
                }
                assert k+1 < trees.length;
            }
        }
        super.handleModelChangedEvent(model, object, index);
        //demoFunction = null;
        fireModelChanged(this);
    }

    protected void handleParameterChangedEvent(Parameter parameter, int index) {
        super.handleParameterChangedEvent(parameter, index);
        if( demoFunction != null ) {
            if( demoFunction == savedDemoFunction ) {
                demoFunction = new VD(demoFunction);
            }
            demoFunction.setDirty();

           // demoFunction = new VD(demoFunction.getTreeIntervals());
        }
        fireModelChanged(this);
    }

    protected void storeState() {
        savedDemoFunction = demoFunction;
    }

    protected void restoreState() {
        demoFunction = savedDemoFunction;
        savedDemoFunction = null;
    }

    public TreeIntervals[] getTreeIntervals() {
        return getDemographicFunction().getTreeIntervals();
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return VariableDemographicModel.MODEL_NAME;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = (XMLObject) xo.getChild(VariableSkylineLikelihood.POPULATION_SIZES);
            Parameter popParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(VariableSkylineLikelihood.INDICATOR_PARAMETER);
            Parameter indicatorParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = (XMLObject) xo.getChild(POPULATION_TREE);

            TreeModel[] treeModels = new TreeModel[cxo.getChildCount()];
            for( int k = 0; k < treeModels.length; ++k) {
                treeModels[k] = (TreeModel) cxo.getChild(k);
            }
            
            Type type = Type.STEPWISE;

            if (xo.hasAttribute(TYPE)) {
                final String s = xo.getStringAttribute(TYPE);
                if (s.equalsIgnoreCase(STEPWISE)) {
                    type = Type.STEPWISE;
                } else if (s.equalsIgnoreCase(LINEAR)) {
                    type = Type.LINEAR;
                } else if (s.equalsIgnoreCase(EXPONENTIAL)) {
                    type = Type.EXPONENTIAL;
                } else {
                    throw new XMLParseException("Unknown Bayesian Skyline type: " + s);
                }
            }

            boolean logSpace = xo.getBooleanAttribute(LOG_SPACE);

            Logger.getLogger("dr.evomodel").info("Variable demographic: " + type.toString() + " control points");

            return new VariableDemographicModel(treeModels, popParam, indicatorParam, type, logSpace);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element represents the likelihood of the tree given the population size vector.";
        }

        public Class getReturnType() {
            return DemographicModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(VariableSkylineLikelihood.TYPE, true),
                new ElementRule(VariableSkylineLikelihood.POPULATION_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(VariableSkylineLikelihood.INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(POPULATION_TREE, new XMLSyntaxRule[]{
                        new ElementRule(TreeModel.class, 1, Integer.MAX_VALUE)
                }),
                AttributeRule.newBooleanRule(LOG_SPACE)
        };
    };


    class VD implements DemographicFunction {
        private double[] values;
        private double[] times;
        private double[] intervals;
        private double[][] ttimes;
        private double[] alltimes;
        private boolean[] dirtyTrees;
        boolean dirty;

        TreeIntervals[] ti;

        public VD(TreeIntervals[] ti) {
            this.ti = ti;
            setup();
        }

        public VD() {
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
            setup();
        }

        public VD(VD demoFunction) {
            this.ti = demoFunction.ti.clone();
            this.values = demoFunction.values.clone();
            this.times = demoFunction.times.clone();
            this.intervals = demoFunction.intervals.clone();
            this.ttimes = demoFunction.ttimes.clone();
            this.alltimes = demoFunction.alltimes.clone();
            this.dirtyTrees = demoFunction.dirtyTrees.clone();
            this.dirty = demoFunction.dirty;
        }

        public void treeChanged(int nt) {
          dirtyTrees[nt] = true;
          setDirty();
        }

        public void setDirty() {
            dirty = true;
        }

        private boolean setTreeTimes(int nt) {
            if( dirtyTrees[nt] ) {
                ti[nt] = new TreeIntervals(trees[nt]);
                TreeIntervals nti = ti[nt];
                nti.setMultifurcationLimit(1e-9);

                assert ti[nt].getIntervalCount() ==  ttimes[nt].length;

                int iCount = 0;
                for(int k = 0; k < ttimes[nt].length; ++k) {
                    double timeToCoal = nti.getInterval(iCount);
                    int linAtStart = nti.getLineageCount(iCount);
                    ++iCount;

                    assert ! ( iCount == nti.getIntervalCount() && linAtStart != 2);

                    int linAtEnd = iCount == nti.getIntervalCount() ? 1 : nti.getLineageCount(iCount);

                    while( linAtStart <= linAtEnd ) {
                        ++iCount;
                        timeToCoal += nti.getInterval(iCount);

                        linAtStart = linAtEnd;
                        ++iCount;
                        linAtEnd = nti.getLineageCount(iCount);
                    }
                    ttimes[nt][k] = timeToCoal + (k == 0 ? 0 : ttimes[nt][k-1]);
                }
                dirtyTrees[nt] = false;
                return true;
            }
            return false;
        }

        private void setup() {
            if( dirty ) {
                boolean any = false;
                for(int nt = 0; nt < ti.length; ++nt) {
                    if( setTreeTimes(nt) ) {
                        any = true;
                    }
                }

                final int nd = indicatorParameter.getDimension();

                assert nd == alltimes.length-1;

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
                                if( ttimes[l][inds[l]] <  ttimes[j][inds[j]] ) {
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
                    if( indicatorParameter.getParameterValue(k) > 0 ) {
                        ++tot;
                    }
                }

                times = new double[tot+1];
                values = new double[tot];
                intervals = new double[tot -1];

                times[0] = 0.0;
                times[tot] = Double.POSITIVE_INFINITY;

                values[0] = popSizeParameter.getParameterValue(0);


                int n = 0;
                for(int k = 0; k < nd && n+1 < tot; ++k) {

                    if( indicatorParameter.getParameterValue(k) > 0 ) {
                        times[n+1] = alltimes[k];

                        values[n+1] = popSizeParameter.getParameterValue(k+1);
                        intervals[n] = times[n+1] - times[n];
                        ++n;
                    }
                }
                dirty = false;
            }
        }

        private int getIntervalIndex(final double t) {
            int j = 0;
            // ugly hack,
            // when doubles are added in a different order and compared later, they can be a tiny bit off. With a
            // stepwise model this creates a "one off" situation here, which is unpleasent.
            // use float comarison here to avoid it
            final float tf = (float)t;
            while( tf > (float)times[j+1] ) ++j;
            return j;
        }

        public double getDemographic(double t) {
            final int j = getIntervalIndex(t);
            switch( type ) {
                case STEPWISE:
                {
                    return values[j];
                }
                case LINEAR:
                {
                    if( j == values.length - 1 ) return values[j];

                    final double a = (t - times[j]) / (intervals[j]);
                    return values[j] + a * (values[j+1] - values[j]);
                }
            }
            return 0;
        }

        public double getIntensity(double t) {
            final int j = getIntervalIndex(t);
            double intensity = 0.0;
            switch( type ) {
                case STEPWISE:
                {
                    for(int k = 0; k < j; ++k) {
                        intensity += intervals[k]/values[k];
                    }
                    intensity += (t - values[j]) / values[j];
                    break;
                }
            }
            return intensity;
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
            double pop0 = popStart + ((start - time0) / interval) * popDiff;
            double pop1 = popStart + ((end - time0) / interval) * popDiff;
            if( pop0 == pop1 ) {
                // either dx == 0 or very small (numerical inaccuracy)
                return dx/pop0;
            }
            return dx * Math.log(pop1/pop0) / (pop1 - pop0);
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

       // private double populationLin()
        public double getIntegral(double start, double finish) {
            final int first = getIntervalIndex(start);
            final int last = getIntervalIndex(finish);
            double intensity = 0.0;

           final double popStart = values[first];
           switch( type ) {
                case STEPWISE:
                {
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
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public String getArgumentName(int n) {
            assert false;
            return null;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public double getArgument(int n) {
            assert false;
            return 0;  //To change body of implemented methods use File | Settings | File Templates.
        }

        public void setArgument(int n, double value) {
            assert false;
           //To change body of implemented methods use File | Settings | File Templates.
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

        public Type getUnits() {
            return trees[0].getUnits();
        }

        public void setUnits(Type units) {
          assert false;
        }

        public TreeIntervals[] getTreeIntervals() {
            return ti;
        }

    }
}
