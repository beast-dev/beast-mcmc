/*
 * MascotPreparedInput.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

import dr.evomodel.coalescent.StructuredCoalescentSchedule;
import dr.evomodel.coalescent.StructuredCoalescentTipData;

/**
 * Validated fixed-tree MASCOT input: biological schedule plus fixed tip payload.
 */
public final class MascotPreparedInput {

    final StructuredCoalescentSchedule schedule;
    final StructuredCoalescentTipData tipData;
    final int nodeCount;
    final int maxLineageId;

    private MascotPreparedInput(StructuredCoalescentSchedule schedule, StructuredCoalescentTipData tipData) {
        this.schedule = schedule;
        this.tipData = tipData;
        this.nodeCount = tipData.nodeCount;
        this.maxLineageId = tipData.nodeCount - 1;
    }

    public static MascotPreparedInput prepare(StructuredCoalescentSchedule schedule,
                                              StructuredCoalescentTipData tipData) {
        if (schedule == null) {
            throw new IllegalArgumentException("schedule must not be null");
        }
        if (tipData == null) {
            throw new IllegalArgumentException("tipData must not be null");
        }
        tipData.validateCompatibleWith(schedule);
        return new MascotPreparedInput(schedule, tipData);
    }
}
