/*
 * AntigenicGradientWrtParameter.java
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

package dr.evomodel.antigenic;

import dr.inference.model.Parameter;

import java.util.Arrays;

public interface AntigenicGradientWrtParameter {

    boolean requiresLocationGradient();

    boolean requiresObservationGradient();

    int getSize();

    void getGradient(double[] gradient, int offset,
                     double[] locationGradient,
                     double[] observationGradient);

    Parameter getParameter();

    class VirusLocations extends Locations {

        VirusLocations(int viruses, int sera, int mdsDim, Parameter parameter, NewAntigenicLikelihood.Layout layout,
                       int startOffset, int tipSize) {
            super(viruses, sera, mdsDim, parameter, layout);
            this.startOffset = startOffset;
            this.tipSize = tipSize;
        }

        @Override
        public void getGradient(double[] gradient, int offset,
                                double[] locationGradient,
                                double[] observationGradient) {
            if (tipSize == mdsDim) {
                super.getGradient(gradient, offset, locationGradient, observationGradient);
            } else {

                final int locationGradientOffset = getLocationOffset();

                Arrays.fill(gradient, offset, offset + getSize(), 0.0);

                for (int v = 0; v < viruses; ++v) {
                    for (int j = 0; j < mdsDim; ++j) {
                        gradient[offset + v * tipSize + startOffset + j] =
                                locationGradient[locationGradientOffset + v * mdsDim + j];
                    }
                }
            }
        }

        @Override
        int getLocationOffset() {
            return layout.getVirusLocationOffset();
        }

        @Override
        public int getSize() {
            return viruses * tipSize;
        }

        private final int startOffset;
        private final int tipSize;
    }

    class SerumLocations extends Locations {

        SerumLocations(int viruses, int sera, int mdsDim, Parameter parameter, NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, parameter, layout);
        }

        @Override
        int getLocationOffset() {
            return layout.getSerumLocationOffset();
        }

        @Override
        public int getSize() {
            return sera * mdsDim;
        }
    }

    abstract class Locations extends Base {

        final NewAntigenicLikelihood.Layout layout;

        Locations(int viruses, int sera, int mdsDim,
                  Parameter parameter,
                  NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, parameter);
            this.layout = layout;
        }

        abstract int getLocationOffset();

        abstract public int getSize();

        @Override
        public boolean requiresLocationGradient() {
            return true;
        }

        @Override
        public boolean requiresObservationGradient() {
            return false;
        }

        @Override
        public void getGradient(double[] gradient, int offset,
                                double[] locationGradient,
                                double[] observationGradient) {
            System.arraycopy(locationGradient, getLocationOffset(), gradient, offset, getSize());
        }
    }

    class Drift extends Locations {

        private final Parameter virusTime;
        private final Parameter serumTime;

        Drift(int viruses, int sera, int mdsDim,
              Parameter locationDrift, Parameter virusTime, Parameter serumTime,
              NewAntigenicLikelihood.Layout layout) {
            super(viruses, sera, mdsDim, locationDrift, layout);
            this.virusTime = virusTime;
            this.serumTime = serumTime;
        }

        @Override
        int getLocationOffset() {
            throw new RuntimeException("Should not be called");
        }

        @Override
        public int getSize() { return 1; }

        @Override
        public void getGradient(double[] gradient, int offset,
                                double[] locationGradient, double[] observationGradient) {

            double derivative = 0;

            int virusOffset = layout.getVirusLocationOffset();
            for (int i = 0; i < viruses; ++i) {
                derivative += locationGradient[virusOffset + i * mdsDim] * virusTime.getParameterValue(i);
            }

            int serumOffset = layout.getSerumLocationOffset();
            for (int i = 0; i < sera; ++i) {
                derivative += locationGradient[serumOffset + i * mdsDim] * serumTime.getParameterValue(i);
            }

            gradient[offset] = derivative;
        }
    }

    abstract class Base implements AntigenicGradientWrtParameter {

        Base(int viruses, int sera, int mdsDim, Parameter parameter) {
            this.viruses = viruses;
            this.sera = sera;
            this.mdsDim = mdsDim;
            this.parameter = parameter;
        }

        @Override
        public Parameter getParameter() {
            return parameter;
        }

        final int viruses;
        final int sera;
        final int mdsDim;
        final Parameter parameter;
    }
}
