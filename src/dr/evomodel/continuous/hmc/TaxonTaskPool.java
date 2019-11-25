package dr.evomodel.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;
import java.util.function.BinaryOperator;

/**
 * @author Marc A. Suchard
 * @author Andrew Holbrook
 */

public class TaxonTaskPool<E> {

    class TaxonTaskIndices {

        final int start;
        final int stop;
        final int task;

        TaxonTaskIndices(int start, int stop, int task) {
            this.start = start;
            this.stop = stop;
            this.task = task;
        }

        public String toString() {
            return start + " " + stop;
        }
    }

    private ExecutorService pool = null;
    final private List<TaxonTaskIndices> indices;
    final private int taxonCount;
    final private int threadCount;

    public TaxonTaskPool(int taxonCount, int threadCount) {
        this.indices = setupTasks(taxonCount, Math.abs(threadCount));
        this.taxonCount = taxonCount;
        this.threadCount = threadCount;
    }

    public ExecutorService getPool() { return pool; }

    public List<TaxonTaskIndices> getIndices() { return indices; }

    public int getNumThreads() { return indices.size(); }

    public int getNumTaxon() { return taxonCount; }

    private List<TaxonTaskIndices> setupTasks(int taxonCount, int threadCount) {
        List<TaxonTaskIndices> tasks = new ArrayList<>(threadCount);

        int length = taxonCount / threadCount;
        if (taxonCount % threadCount != 0) ++length;

        int start = 0;

        for (int task = 0; task < threadCount && start < taxonCount; ++task) {
            tasks.add(new TaxonTaskIndices(start, Math.min(start + length, taxonCount), task));
            start += length;
        }

        return tasks;
    }

    private ExecutorService setupParallelServices(int threadCount) {

        final ExecutorService pool;

        if (threadCount > 1) {
            pool = Executors.newFixedThreadPool(threadCount);
        } else if (threadCount < 0) {
            pool = Executors.newCachedThreadPool();
        } else {
            pool = null;
        }

        return pool;
    }

    public interface TaxonCallable {
        void execute(int taxon, int thread);
    }

    public interface RangeCallable<E> {
        E map(int start, int end, int thread);
    }

    public E mapReduce(final RangeCallable<E> map, final BinaryOperator<E> reduce) {

        E result = null;

        if (indices.size() == 1) {

            final TaxonTaskIndices index = indices.get(0);
            result = map.map(index.start, index.stop, 0);

        } else {

            if (pool == null) {
                pool = setupParallelServices(threadCount);
            }

            List<Callable<E>> calls = new ArrayList<>();

            for (final TaxonTaskIndices indexSet : indices) {
                calls.add(() -> map.map(indexSet.start, indexSet.stop, indexSet.task));
            }

            try {

                List<Future<E>> futures = pool.invokeAll(calls);

                result = futures.get(0).get();
                for (int i = 1; i < futures.size(); ++i) {
                    result = reduce.apply(result, futures.get(i).get());
                }

            } catch (InterruptedException | ExecutionException exception) {
                exception.printStackTrace();
            }

        }

        return result;
    }

    public void fork(final TaxonCallable runnable) {
        if (indices.size() == 1) {

            final TaxonTaskIndices index = indices.get(0);
            for (int taxon = index.start; taxon < index.stop; ++taxon) {
                runnable.execute(taxon, 0);
            }

        } else {

            if (pool == null) {
                pool = setupParallelServices(threadCount);
            }

            List<Callable<Object>> calls = new ArrayList<>();

            for (final TaxonTaskIndices indexSet : indices) {

                calls.add(Executors.callable(() -> {
                            for (int taxon = indexSet.start; taxon < indexSet.stop; ++taxon) {
                                runnable.execute(taxon, indexSet.task);
                            }
                        }
                ));
            }

            try {
                pool.invokeAll(calls);
            } catch (InterruptedException exception) {
                exception.printStackTrace();
            }
        }
    }

    private static final String PARSER_NAME = "taxonTaskPool";
    private static final String THREAD_COUNT = "threadCount";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

         @Override
         public Object parseXMLObject(XMLObject xo) throws XMLParseException {

             Tree tree = (Tree) xo.getChild(Tree.class);
             int threadCount = xo.getAttribute(THREAD_COUNT, 1);
             return new TaxonTaskPool(tree.getExternalNodeCount(), threadCount);
         }

         @Override
         public XMLSyntaxRule[] getSyntaxRules() {
             return rules;
         }

         @Override
         public String getParserDescription() {
             return "A thread pool for per-taxon specific operations";
         }

         @Override
         public Class getReturnType() {
             return TaxonTaskPool.class;
         }

         @Override
         public String getParserName() {
             return PARSER_NAME;
         }

         private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
                 AttributeRule.newIntegerRule(THREAD_COUNT, true),
                 new ElementRule(Tree.class),
         };
     };
}
