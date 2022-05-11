/*
 * MassPreconditionScheduler.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.inference.operators.hmc;

import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.MCMCOperator;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public interface MassPreconditionScheduler {

    boolean shouldUpdatePreconditioning();

    boolean shouldStoreSecant(double[] lastGradient, double[] lastPosition);

    void forceUpdateCount();

    enum Type {

        NONE("none") {
            @Override
            public MassPreconditionScheduler factory(MassPreconditioningOptions options,
                                                     MCMCOperator operator) {
                return new None();
            }
        },
        DEFAULT("default") {
            @Override
            public MassPreconditionScheduler factory(MassPreconditioningOptions options,
                                                     MCMCOperator operator) {
                return new Default(options, operator);
            }
        };

        private final String name;

        Type(String name) {
            this.name = name;
        }

        public abstract MassPreconditionScheduler factory(MassPreconditioningOptions options,
                                                          MCMCOperator operator);


        public String getName() { return name; }

        public static MassPreconditionScheduler.Type parseFromString(String text) {
            for (MassPreconditionScheduler.Type type : MassPreconditionScheduler.Type.values()) {
                if (type.name.toLowerCase().compareToIgnoreCase(text) == 0) {
                    return type;
                }
            }
            return MassPreconditionScheduler.Type.NONE;
        }
    }

    class None implements MassPreconditionScheduler {

        @Override
        public boolean shouldUpdatePreconditioning() {
            return false;
        }

        @Override
        public boolean shouldStoreSecant(double[] lastGradient, double[] lastPosition) {
            return false;
        }

        @Override
        public void forceUpdateCount() {
        }
    }

    class Default implements MassPreconditionScheduler {

        private MassPreconditioningOptions options;
        private MCMCOperator operator;
        private int totalUpdates = 0;
        private long paramUpdateCount = 0;
        private boolean useOperatorCount = true;

        Default(MassPreconditioningOptions options,
                MCMCOperator operator) {
            this.options = options;
            this.operator = operator;
        }

        @Override
        public boolean shouldUpdatePreconditioning() {

            long count = useOperatorCount ? operator.getCount() : this.paramUpdateCount;
            boolean shouldUpdate = shouldUpdate(count);

            if (shouldUpdate) {
                totalUpdates++;
            }

            return shouldUpdate;
        }

        protected boolean shouldUpdate(long count) {
            return ((options.preconditioningUpdateFrequency() > 0)
                    && (((count % options.preconditioningUpdateFrequency() == 0)))
                    && (options.preconditioningMaxUpdate() == 0 || totalUpdates < options.preconditioningMaxUpdate()));
        }

        @Override
        public boolean shouldStoreSecant(double[] lastGradient, double[] lastPosition) {
            return lastGradient != null && lastPosition != null;
        }

        @Override
        public void forceUpdateCount() {
            paramUpdateCount++;
            useOperatorCount = false;
        }
    }

    class UpdateByProbability extends Default {

        UpdateByProbability(HamiltonianMonteCarloOperator.Options options, AdaptableMCMCOperator operator) {
            super(options, operator);
        }
    }

}
