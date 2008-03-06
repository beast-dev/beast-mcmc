package dr.evomodel.coalescent;

import dr.evolution.coalescent.MultiLociTreeSet;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * @author Joseph Heled
 * @version $Id$
 */
public class VariableDemographicModel extends DemographicModel implements MultiLociTreeSet {
    static final String MODEL_NAME = "variableDemographic";
    public static final String POPULATION_SIZES = "populationSizes";
    public static final String INDICATOR_PARAMETER = "indicators";
    public static final String POPULATION_TREES = "trees";
    private static String PLOIDY = "ploidy";
    public static String POP_TREE = "ptree";

    public static final String LOG_SPACE = "logUnits";

    public static final String TYPE = "type";
    public static final String STEPWISE = "stepwise";
    public static final String LINEAR = "linear";
    public static final String EXPONENTIAL = "exponential";

    private Parameter popSizeParameter;
    private Parameter indicatorParameter;
    private Type type;
    private boolean logSpace;
    private TreeModel[] trees;
    private VDdemographicFunction demoFunction = null;
    private VDdemographicFunction savedDemoFunction = null;
    private double[] populationFactors;

    public Parameter getIndices() {
        return indicatorParameter;
    }

    public enum Type {
        STEPWISE,
        LINEAR,
        EXPONENTIAL
    }

    public VariableDemographicModel(TreeModel[] trees, double[] popFactors, Parameter popSizeParameter, Parameter indicatorParameter,
                                    Type type, boolean logSpace) {
        super(MODEL_NAME);

        this.popSizeParameter = popSizeParameter;
        this.indicatorParameter = indicatorParameter;

        this.populationFactors = popFactors;

        //final int redcueDim = type == Type.STEPWISE ? 1 : 0;
        int events = 0;
        for( Tree t : trees ) {
            // number of coalescent envents
            events += t.getExternalNodeCount() - 1;
            // we will have to handle this I guess
            assert t.getUnits() == trees[0].getUnits();
        }
        // all trees share time 0, need fixing for serial data

        events +=  type == Type.STEPWISE ? 0 : 1;

        final int popSizes = popSizeParameter.getDimension();
        final int nIndicators = indicatorParameter.getDimension();
        this.type = type;
        this.logSpace = logSpace;
        if( logSpace ) {
            throw new IllegalArgumentException("sorry log space not implemented");
        }

        if (popSizes != events) {
            throw new IllegalArgumentException("Dimension of population parameter (" + popSizes +
                    ") must be the same as the number of internal nodes in the tree. (" + events + ")");
        }
        
        if (nIndicators != events - 1) {
            throw new IllegalArgumentException("Dimension of indicator parameter must one less than the number of internal nodes in the tree. ("
            + nIndicators + " != " + (events-1) + ")");
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

    // must be called before each use demographic integrals/population to set the population scaling
    public TreeIntervals getTreeIntervals(int nt) {
        return getDemographicFunction().getTreeIntervals(nt);
    }

    public double getPopulationFactor(int nt) {
        return populationFactors[nt];
    }

    public void storeTheState() {
        // as a demographic model store/restore is already taken care of 
    }

    public void restoreTheState() {
        // as a demographic model store/restore is already taken care of
    }

    public VDdemographicFunction getDemographicFunction() {
        if( demoFunction == null ) {
            demoFunction = new VDdemographicFunction(trees, type , /*populationFactors,*/
                    indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues());
        } else {
            demoFunction.setup(trees, indicatorParameter.getParameterValues(), popSizeParameter.getParameterValues());
        }
        return demoFunction;
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // tree has changed
        //System.out.println("model changed: " + model);
        if( demoFunction != null ) {
            if( demoFunction == savedDemoFunction ) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            for(int k = 0; k < trees.length; ++k) {
                if( model == trees[k] ) {
                    demoFunction.treeChanged(k);
                    //System.out.println("tree changed: " + k + " " + Arrays.toString(demoFunction.dirtyTrees)
                     //       + " " + demoFunction.dirtyTrees);
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
        //System.out.println("parm changed: " + parameter);
        super.handleParameterChangedEvent(parameter, index);
        if( demoFunction != null ) {
            if( demoFunction == savedDemoFunction ) {
                demoFunction = new VDdemographicFunction(demoFunction);
            }
            demoFunction.setDirty();
        }
        fireModelChanged(this);
    }

    protected void storeState() {
        //System.out.println("store");

        /*for (int u = 0; u < 2; ++u) {
            System.out.println(u + " " + getModel(u));
          System.out.println(Arrays.toString(demoFunction.ttimes[u]));
        }
*/
        savedDemoFunction = demoFunction;
    }

    protected void restoreState() {
        //System.out.println("restore");
        demoFunction = savedDemoFunction;
        savedDemoFunction = null;
        
        /*for (int u = 0; u < 2; ++u) {
            System.out.println(u + " " + getModel(u));
            System.out.println(Arrays.toString(demoFunction.ttimes[u]));
        }*/
        //demoFunction.setDirty(); demoFunction.setup();
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

            cxo = (XMLObject) xo.getChild(POPULATION_TREES);

            final int nc = cxo.getChildCount();
            TreeModel[] treeModels = new TreeModel[nc];
            double[] populationFactor = new double[nc];
            //Parameter[] branchFactors = new Parameter[nc];

            for(int k = 0; k < treeModels.length; ++k) {
                final XMLObject child = (XMLObject) cxo.getChild(k);

               // populationFactor[k] = 1.0;
                //branchFactors[k] = null;

//                if( child.getClass().isAssignableFrom(TreeModel.class) ) {
//                    treeModels[k] = (TreeModel) child;
//                } else {
                    /*final XMLObject bscale = (XMLObject) cxo.getChild("branchScale");

                    if( bscale != null ) {
                      branchFactors[k] = (Parameter) bscale.getChild(Parameter.class);
                    }
                    */
              //  cxo = (XMLObject) child;

                populationFactor[k] = child.hasAttribute(PLOIDY) ? child.getDoubleAttribute(PLOIDY) : 1.0 ;

                treeModels[k] = (TreeModel) child.getChild(TreeModel.class);
               // }

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

            boolean logSpace = xo.hasAttribute(LOG_SPACE) ? xo.getBooleanAttribute(LOG_SPACE) : false;

            Logger.getLogger("dr.evomodel").info("Variable demographic: " + type.toString() + " control points");

            return new VariableDemographicModel(treeModels, populationFactor, popParam, indicatorParam, type, logSpace);
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

        private XMLSyntaxRule[] rules =
        new XMLSyntaxRule[]{
                AttributeRule.newStringRule(VariableSkylineLikelihood.TYPE, true),
                 AttributeRule.newBooleanRule(LOG_SPACE, true),

                new ElementRule(VariableSkylineLikelihood.POPULATION_SIZES, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(VariableSkylineLikelihood.INDICATOR_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(POPULATION_TREES, new XMLSyntaxRule[]{
                        new ElementRule(POP_TREE, new XMLSyntaxRule[] {
                              AttributeRule.newDoubleRule(PLOIDY, true),
                              new ElementRule(TreeModel.class),
                        }, 1, Integer.MAX_VALUE)
                })

//                        new XMLSyntaxRule[]{
//                           new OrRule(new ElementRule(TreeModel.class),
//                                      new ElementRule(POP_TREE,
//                                              new XMLSyntaxRule[]{ new ElementRule(TreeModel.class),
//                                                                   AttributeRule.newDoubleRule(PLOIDY) }))
//                        }, 1, Integer.MAX_VALUE) ,
                       /* new ElementRule(new XMLSyntaxRule[]{new OrRule(new ElementRule(TreeModel.class),
                                new ElementRule(TreeModel.class))} }, 1, Integer.MAX_VALUE)*/


        };
    };

   /*
    class VD implements DemographicFunction {
        private double[] values;
        private double[] times;
        private double[] intervals;
        private double[][] ttimes;
        private double[] alltimes;
        private boolean[] dirtyTrees;
        boolean dirty;

        TreeIntervals[] ti;
        private double popFactor = 1;

//        public VD(TreeIntervals[] ti) {
//            this.ti = ti;
//            setup();
//        }

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
            for(int k = 0; k < ttimes.length; ++k) {
                ttimes[k] = ttimes[k].clone();
            }

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
//                double[] doubles = null;
//                if( ! dirtyTrees[nt] ) {
//                   doubles = ttimes[nt].clone();
//
//                }
                ti[nt] = new TreeIntervals(trees[nt]);

                TreeIntervals nti = ti[nt];
                // make sure we get each coalescent event individually
                nti.setMultifurcationLimit(0);

                final int nLineages = nti.getIntervalCount();
                assert nLineages == ttimes[nt].length: nLineages + " " + ttimes[nt].length;

                int iCount = 0;
                for(int k = 0; k < ttimes[nt].length; ++k) {
                    double timeToCoal = nti.getInterval(iCount);
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

//                if( doubles != null ) {
//                    if( ! Arrays.equals(doubles, ttimes[nt]) ) {
//                       System.out.println(Arrays.toString(doubles) + " != " + Arrays.toString(ttimes[nt])
//                               + Arrays.toString(dirtyTrees) + " " + dirtyTrees);
//                    }
//                }
                dirtyTrees[nt] = false;
               // System.out.print(nt + " " + Arrays.toString(dirtyTrees) + " " + dirtyTrees);
                return true;
            }
            return false;
        }

        private void setup() {
           // boolean was = dirty;
            if( dirty ) {
                boolean any = false;
                for(int nt = 0; nt < ti.length; ++nt) {
                    if( setTreeTimes(nt) ) {
                        any = true;
                    }
                }

                final int nd = indicatorParameter.getDimension();

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
            //
//            System.out.println("after setup " + (was ? "(dirty)" : "") + " , alltimes " + Arrays.toString(alltimes)
//                        + " times " + Arrays.toString(times) + " values " + Arrays.toString(values) +
//                    " inds " + Arrays.toString(indicatorParameter.getParameterValues())) ;
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
                    p = a * values[j+1] + (1-a) * values[j]; // values[j] + a * (values[j+1] - values[j]);
                    break;
                }
                default: throw new IllegalArgumentException("");

            }
            return p * popFactor;
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
            
           // final double pop0 = popStart + ((start - time0) / interval) * popDiff;
           // final double pop1 = popStart + ((end - time0) / interval) * popDiff;

            // do same as above more effeciently
            final double r = popDiff / interval;
            final double x = popStart - time0 * r;
            final double pop0 = x + start * r;
            final double pop1 = x + end * r;
            
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
           return intensity/popFactor;
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

        public TreeIntervals getTreeIntervals(int nt) {
            popFactor = populationFactors[nt];
            return ti[nt];
        }
    }
                */
}
