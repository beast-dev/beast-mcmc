/*
 * ModelRegistry.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.inference.model;

import java.util.*;
import java.util.concurrent.CopyOnWriteArrayList;

/**
 * The sets of likelihoods, models and parameters that have been created (FULL) and that are
 * connected to a Markov chain (CONNECTED). These are used for warnings about unused
 * components, for checkpointing and for releasing resources at the end of a run.
 *
 * Each chain of a multi-chain run has its own registry, bound to the threads that build and
 * run the chain. The static sets such as Parameter.CONNECTED_PARAMETER_SET are views of the
 * registry bound to the current thread, or of the default registry when none is bound (as in
 * a single-chain run), so existing code sees only the objects of its own chain.
 *
 * @author Andrew Rambaut
 */
public final class ModelRegistry {

    public ModelRegistry() {
        ALL.add(this);
    }

    /**
     * @return the registry bound to the current thread or the default one
     */
    public static ModelRegistry current() {
        ModelRegistry registry = CURRENT.get();
        return registry != null ? registry : DEFAULT;
    }

    /**
     * Binds a registry to the current thread. Pass null to go back to the default registry.
     */
    public static void setCurrent(ModelRegistry registry) {
        if (registry == null) {
            CURRENT.remove();
        } else {
            CURRENT.set(registry);
        }
    }

    /**
     * @return every registry that has been created, including the default one
     */
    public static List<ModelRegistry> getAll() {
        return Collections.unmodifiableList(ALL);
    }

    public Set<Likelihood> getFullLikelihoodSet() {
        return fullLikelihoods;
    }

    public Set<Likelihood> getConnectedLikelihoodSet() {
        return connectedLikelihoods;
    }

    public Set<Model> getFullModelSet() {
        return fullModels;
    }

    public Set<Model> getConnectedModelSet() {
        return connectedModels;
    }

    public Set<Parameter> getFullParameterSet() {
        return fullParameters;
    }

    public Set<Parameter> getConnectedParameterSet() {
        return connectedParameters;
    }

    interface SetGetter<T> {
        Set<T> get(ModelRegistry registry);
    }

    /**
     * @return a set that forwards to the given set of the current thread's registry
     */
    static <T> Set<T> view(final SetGetter<T> getter) {
        return new AbstractSet<T>() {
            private Set<T> set() {
                return getter.get(current());
            }

            @Override
            public Iterator<T> iterator() {
                return set().iterator();
            }

            @Override
            public int size() {
                return set().size();
            }

            @Override
            public boolean add(T t) {
                return set().add(t);
            }

            @Override
            public boolean contains(Object o) {
                return set().contains(o);
            }

            @Override
            public boolean remove(Object o) {
                return set().remove(o);
            }

            @Override
            public void clear() {
                set().clear();
            }
        };
    }

    // The set types (and so the iteration orders) are those the static sets always had.
    // They are synchronized because objects can be created on likelihood worker threads.
    private final Set<Likelihood> fullLikelihoods = Collections.synchronizedSet(new HashSet<Likelihood>());
    private final Set<Likelihood> connectedLikelihoods = Collections.synchronizedSet(new HashSet<Likelihood>());
    private final Set<Model> fullModels = Collections.synchronizedSet(new LinkedHashSet<Model>());
    private final Set<Model> connectedModels = Collections.synchronizedSet(new LinkedHashSet<Model>());
    private final Set<Parameter> fullParameters = Collections.synchronizedSet(new LinkedHashSet<Parameter>());
    private final Set<Parameter> connectedParameters = Collections.synchronizedSet(new LinkedHashSet<Parameter>());

    private static final List<ModelRegistry> ALL = new CopyOnWriteArrayList<ModelRegistry>();
    private static final ThreadLocal<ModelRegistry> CURRENT = new ThreadLocal<ModelRegistry>();
    private static final ModelRegistry DEFAULT = new ModelRegistry();
}
