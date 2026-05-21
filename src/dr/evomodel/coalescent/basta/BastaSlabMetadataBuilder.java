package dr.evomodel.coalescent.basta;

import java.util.Arrays;

public final class BastaSlabMetadataBuilder {
    public static final int VERSION = 2;

    public static final int FIXED_HEADER_INTS = 11;

    private static final int INITIAL_HISTORY_CAP = 1024;
    private static final int INITIAL_COAL_CAP = 64;
    private static final int INITIAL_SNAP_CAP = 1024;


    private int[] nodeTail;
    private int[] nodeKLen;
    private int[] nodeBornBuf;
    private int[] nodeProducerCoal;
    private boolean[] nodeActive;


    private int[] histBuf;
    private int[] histPrev;
    private int histLen;

    private int coalCountRec;

    private int[] coalIntervalRec;
    private int[] coalDestBufRec;
    private int[] coalLeftAccBufRec;
    private int[] coalRightAccBufRec;
    private int[] coalParentNodeRec;
    private int[] coalLeftChildCoal;
    private int[] coalRightChildCoal;
    private int[] coalLeftBotBufRec;
    private int[] coalRightBotBufRec;
    private int[] coalLeftKbRec;
    private int[] coalRightKbRec;
    private int[] coalLeftSnapStart;
    private int[] coalLeftSnapLen;

    private int[] coalRightSnapStart;
    private int[] coalRightSnapLen;

    private int[] coalSnapArena;
    private int coalSnapArenaLen;

    private long[] forwardBufBitset;
    private int forwardBufMaxIdx;

    private int[] opReduce;
    private int operationCountRecorded;
    private int builderTipCount;
    private int builderTipPartialsBufferCount;
    private int builderIndexOffsetPat;
    private int nodeCount;

    private int intervalCountSlot;

    private int[] payloadBuf;

    private int lastPackedLength = -1;

    private int[] scratchCoalDepth;

    private int[] scratchStackCoal;
    private int[] scratchStackDepth;


    private int[] scratchBranchCursor;
    private int[] scratchCoalCursor;
    private int[] scratchIntCursor;



    public BastaSlabMetadataBuilder() {
    }

    public void beginTraversal(int nodeCount) {
        beginTraversal(nodeCount, Integer.MAX_VALUE, Integer.MAX_VALUE, 1);
    }

    public void beginTraversal(int nodeCount,
                               int tipCount,
                               int tipPartialsBufferCount,
                               int kIndexOffsetPat) {
        this.nodeCount = nodeCount;
        this.builderTipCount = tipCount;
        this.builderTipPartialsBufferCount = tipPartialsBufferCount;
        this.builderIndexOffsetPat = kIndexOffsetPat;

        if (nodeTail == null || nodeTail.length < nodeCount) {
            nodeTail = new int[nodeCount];
            nodeKLen = new int[nodeCount];
            nodeBornBuf = new int[nodeCount];
            nodeProducerCoal = new int[nodeCount];
            nodeActive = new boolean[nodeCount];
        }
        Arrays.fill(nodeTail, 0, nodeCount, -1);
        Arrays.fill(nodeKLen, 0, nodeCount, 0);
        Arrays.fill(nodeBornBuf, 0, nodeCount, -1);
        Arrays.fill(nodeProducerCoal, 0, nodeCount, -1);
        Arrays.fill(nodeActive, 0, nodeCount, false);

        if (histBuf == null) {
            histBuf = new int[INITIAL_HISTORY_CAP];
            histPrev = new int[INITIAL_HISTORY_CAP];
        }
        histLen = 0;

        if (coalIntervalRec == null) {
            allocCoalArrays(INITIAL_COAL_CAP);
        }
        coalCountRec = 0;

        if (coalSnapArena == null) {
            coalSnapArena = new int[INITIAL_SNAP_CAP];
        }
        coalSnapArenaLen = 0;

        if (forwardBufBitset == null) {
            forwardBufBitset = new long[64];
        } else {
            Arrays.fill(forwardBufBitset, 0L);
        }
        forwardBufMaxIdx = -1;

        if (opReduce == null) {
            opReduce = new int[5 * INITIAL_HISTORY_CAP];
        }
        operationCountRecorded = 0;

        intervalCountSlot = 0;
    }

    private void allocCoalArrays(int cap) {
        coalIntervalRec = new int[cap];
        coalDestBufRec = new int[cap];
        coalLeftAccBufRec = new int[cap];
        coalRightAccBufRec = new int[cap];
        coalParentNodeRec = new int[cap];
        coalLeftChildCoal = new int[cap];
        coalRightChildCoal = new int[cap];
        coalLeftBotBufRec = new int[cap];
        coalRightBotBufRec = new int[cap];
        coalLeftKbRec = new int[cap];
        coalRightKbRec = new int[cap];
        coalLeftSnapStart = new int[cap];
        coalLeftSnapLen = new int[cap];
        coalRightSnapStart = new int[cap];
        coalRightSnapLen = new int[cap];
    }

    private void growCoalArrays() {
        int oldCap = coalIntervalRec.length;
        int newCap = oldCap * 2;
        coalIntervalRec = Arrays.copyOf(coalIntervalRec, newCap);
        coalDestBufRec = Arrays.copyOf(coalDestBufRec, newCap);
        coalLeftAccBufRec = Arrays.copyOf(coalLeftAccBufRec, newCap);
        coalRightAccBufRec = Arrays.copyOf(coalRightAccBufRec, newCap);
        coalParentNodeRec = Arrays.copyOf(coalParentNodeRec, newCap);
        coalLeftChildCoal = Arrays.copyOf(coalLeftChildCoal, newCap);
        coalRightChildCoal = Arrays.copyOf(coalRightChildCoal, newCap);
        coalLeftBotBufRec = Arrays.copyOf(coalLeftBotBufRec, newCap);
        coalRightBotBufRec = Arrays.copyOf(coalRightBotBufRec, newCap);
        coalLeftKbRec = Arrays.copyOf(coalLeftKbRec, newCap);
        coalRightKbRec = Arrays.copyOf(coalRightKbRec, newCap);
        coalLeftSnapStart = Arrays.copyOf(coalLeftSnapStart, newCap);
        coalLeftSnapLen = Arrays.copyOf(coalLeftSnapLen, newCap);
        coalRightSnapStart = Arrays.copyOf(coalRightSnapStart, newCap);
        coalRightSnapLen = Arrays.copyOf(coalRightSnapLen, newCap);
    }


    private int transformBuffer(int raw) {
        if (raw < 0) {
            return raw;
        }
        return builderIndexOffsetPat * raw;
    }

    public void onSample(int tipNodeNumber) {
        if (tipNodeNumber < 0 || tipNodeNumber >= nodeCount) {
            throw new IllegalArgumentException("Tip node out of range: " + tipNodeNumber + " (nodeCount=" + nodeCount + ")");
        }
        if (nodeActive[tipNodeNumber]) {
            throw new IllegalStateException("Tip already active: node " + tipNodeNumber);
        }
        final int tipBuf = transformBuffer(tipNodeNumber);
        nodeBornBuf[tipNodeNumber] = tipBuf;
        nodeTail[tipNodeNumber] = -1;
        nodeKLen[tipNodeNumber] = 0;
        nodeProducerCoal[tipNodeNumber] = -1;
        nodeActive[tipNodeNumber] = true;
        markBufferUsed(tipBuf);
    }

    public void onPropagate(int nodeNumber, int outputBuffer, int inputBuffer, int subInterval) {
        if (nodeNumber < 0 || nodeNumber >= nodeCount || !nodeActive[nodeNumber]) {
            throw new IllegalStateException("Propagate on inactive node: " + nodeNumber);
        }
        final int xfOutput = transformBuffer(outputBuffer);
        final int xfInput = transformBuffer(inputBuffer);

        ensureHistoryCapacity(histLen + 1);
        final int slot = histLen++;
        histBuf[slot] = xfOutput;
        histPrev[slot] = nodeTail[nodeNumber];
        nodeTail[nodeNumber] = slot;
        nodeKLen[nodeNumber]++;
        markBufferUsed(xfOutput);
        markBufferUsed(xfInput);
        if (subInterval + 1 > intervalCountSlot) {
            intervalCountSlot = subInterval + 1;
        }

        ensureReduceCapacity(operationCountRecorded + 1);
        final int rbase = 5 * operationCountRecorded;
        opReduce[rbase] = xfInput;
        opReduce[rbase + 1] = -1;
        opReduce[rbase + 2] = xfOutput;
        opReduce[rbase + 3] = -1;
        opReduce[rbase + 4] = subInterval;
        operationCountRecorded++;
    }

    public void onCoalescent(int parentNodeNumber,
                             int leftChildNodeNumber,
                             int rightChildNodeNumber,
                             int outputBuffer,
                             int leftInputBuffer,
                             int rightInputBuffer,
                             int leftAccBuffer,
                             int rightAccBuffer,
                             int subInterval) {
        if (!nodeActive[leftChildNodeNumber] || !nodeActive[rightChildNodeNumber]) {
            throw new IllegalStateException("Coalescent on inactive child(ren): L=" + leftChildNodeNumber + " R=" + rightChildNodeNumber);
        }
        if (parentNodeNumber < 0 || parentNodeNumber >= nodeCount) {
            throw new IllegalArgumentException("Parent node out of range: " + parentNodeNumber);
        }

        if (coalCountRec == coalIntervalRec.length) {
            growCoalArrays();
        }
        final int ci = coalCountRec++;


        final int xfOutput = transformBuffer(outputBuffer);
        final int xfLeftIn = transformBuffer(leftInputBuffer);
        final int xfRightIn = transformBuffer(rightInputBuffer);
        final int xfLeftAcc = transformBuffer(leftAccBuffer);
        final int xfRightAcc = transformBuffer(rightAccBuffer);

        coalIntervalRec[ci] = subInterval;
        coalDestBufRec[ci] = xfOutput;
        coalLeftAccBufRec[ci] = xfLeftAcc;
        coalRightAccBufRec[ci] = xfRightAcc;
        coalParentNodeRec[ci] = parentNodeNumber;
        coalLeftChildCoal[ci] = nodeProducerCoal[leftChildNodeNumber];
        coalRightChildCoal[ci] = nodeProducerCoal[rightChildNodeNumber];

        final int leftHistLen = nodeKLen[leftChildNodeNumber];
        final int leftSnapBegin = coalSnapArenaLen;
        ensureSnapCapacity(coalSnapArenaLen + leftHistLen);
        int leftU1 = nodeBornBuf[leftChildNodeNumber];
        {
            int idx = nodeTail[leftChildNodeNumber];
            if (idx >= 0) leftU1 = histBuf[idx];
            for (int k = 0; k < leftHistLen; k++) {
                coalSnapArena[coalSnapArenaLen++] = histBuf[idx];
                idx = histPrev[idx];
            }
            if (idx != -1) {
                throw new IllegalStateException(
                        "Left history walk overran nodeKLen at coal ci=" + ci
                                + " leftChild=" + leftChildNodeNumber
                                + " expected " + leftHistLen + " entries");
            }
        }
        coalLeftSnapStart[ci] = leftSnapBegin;
        coalLeftSnapLen[ci] = leftHistLen;
        coalLeftKbRec[ci] = leftHistLen + 1;
        coalLeftBotBufRec[ci] = nodeBornBuf[leftChildNodeNumber];

        // Validation: u_1 must equal the input buffer the coalescent op declares.
        if (leftU1 != xfLeftIn) {
            throw new IllegalStateException(
                    "Left u_1 mismatch at coal ci=" + ci
                            + " parent=" + parentNodeNumber
                            + " leftChild=" + leftChildNodeNumber
                            + " expected u_1=" + leftU1
                            + " actual xfLeftIn=" + xfLeftIn
                            + " (raw leftInputBuffer=" + leftInputBuffer + ")");
        }

        // Same for RIGHT child.
        final int rightHistLen = nodeKLen[rightChildNodeNumber];
        final int rightSnapBegin = coalSnapArenaLen;
        ensureSnapCapacity(coalSnapArenaLen + rightHistLen);
        int rightU1 = nodeBornBuf[rightChildNodeNumber];
        {
            int idx = nodeTail[rightChildNodeNumber];
            if (idx >= 0) rightU1 = histBuf[idx];
            for (int k = 0; k < rightHistLen; k++) {
                coalSnapArena[coalSnapArenaLen++] = histBuf[idx];
                idx = histPrev[idx];
            }
            if (idx != -1) {
                throw new IllegalStateException(
                        "Right history walk overran nodeKLen at coal ci=" + ci
                                + " rightChild=" + rightChildNodeNumber
                                + " expected " + rightHistLen + " entries");
            }
        }
        coalRightSnapStart[ci] = rightSnapBegin;
        coalRightSnapLen[ci] = rightHistLen;
        coalRightKbRec[ci] = rightHistLen + 1;
        coalRightBotBufRec[ci] = nodeBornBuf[rightChildNodeNumber];

        if (rightU1 != xfRightIn) {
            throw new IllegalStateException(
                    "Right u_1 mismatch at coal ci=" + ci
                            + " parent=" + parentNodeNumber
                            + " rightChild=" + rightChildNodeNumber
                            + " expected u_1=" + rightU1
                            + " actual xfRightIn=" + xfRightIn
                            + " (raw rightInputBuffer=" + rightInputBuffer + ")");
        }

        markBufferUsed(xfOutput);
        markBufferUsed(xfLeftIn);
        markBufferUsed(xfRightIn);
        markBufferUsed(xfLeftAcc);
        markBufferUsed(xfRightAcc);

        nodeActive[leftChildNodeNumber] = false;
        nodeActive[rightChildNodeNumber] = false;


        nodeBornBuf[parentNodeNumber] = xfOutput;
        nodeTail[parentNodeNumber] = -1;
        nodeKLen[parentNodeNumber] = 0;
        nodeProducerCoal[parentNodeNumber] = ci;
        nodeActive[parentNodeNumber] = true;

        if (subInterval + 1 > intervalCountSlot) {
            intervalCountSlot = subInterval + 1;
        }

        ensureReduceCapacity(operationCountRecorded + 1);
        final int rbase = 5 * operationCountRecorded;
        opReduce[rbase] = xfLeftIn;
        opReduce[rbase + 1] = xfRightIn;
        opReduce[rbase + 2] = xfLeftAcc;
        opReduce[rbase + 3] = xfRightAcc;
        opReduce[rbase + 4] = subInterval;
        operationCountRecorded++;
    }

    private void ensureHistoryCapacity(int needed) {
        if (needed > histBuf.length) {
            int newCap = histBuf.length;
            while (newCap < needed) newCap *= 2;
            histBuf = Arrays.copyOf(histBuf, newCap);
            histPrev = Arrays.copyOf(histPrev, newCap);
        }
    }

    private void ensureSnapCapacity(int needed) {
        if (needed > coalSnapArena.length) {
            int newCap = coalSnapArena.length;
            while (newCap < needed) newCap *= 2;
            coalSnapArena = Arrays.copyOf(coalSnapArena, newCap);
        }
    }

    private void ensureReduceCapacity(int opsNeeded) {
        final int needed = 5 * opsNeeded;
        if (needed > opReduce.length) {
            int newCap = opReduce.length;
            while (newCap < needed) newCap *= 2;
            opReduce = Arrays.copyOf(opReduce, newCap);
        }
    }

    private void markBufferUsed(int buf) {
        if (buf < 0) {
            return;
        }
        final int word = buf >>> 6;
        if (word >= forwardBufBitset.length) {
            int newLen = forwardBufBitset.length;
            while (newLen <= word) newLen *= 2;
            forwardBufBitset = Arrays.copyOf(forwardBufBitset, newLen);
        }
        forwardBufBitset[word] |= 1L << (buf & 63);
        if (buf > forwardBufMaxIdx) {
            forwardBufMaxIdx = buf;
        }
    }

    public int[] buildAndPackDirect(int slabOpsPerBlock) {
        if (slabOpsPerBlock <= 0) {
            throw new IllegalArgumentException("slabOpsPerBlock must be > 0, got " + slabOpsPerBlock);
        }
        if (coalCountRec == 0) {
            return emptyPackedDirect(slabOpsPerBlock);
        }

        final int coalCount = coalCountRec;
        final int branchCount = 2 * coalCount;
        final int rootCoalIdx = coalCount - 1;

        if (scratchCoalDepth == null || scratchCoalDepth.length < coalCount) {
            scratchCoalDepth = new int[coalCount + (coalCount >> 2)];
        }

        java.util.Arrays.fill(scratchCoalDepth, 0, coalCount, -1);
        scratchCoalDepth[rootCoalIdx] = 0;
        int maxDepth = 0;
        for (int ci = coalCount - 1; ci >= 0; ci--) {
            final int d = scratchCoalDepth[ci];
            if (d < 0) continue;
            final int leftCC = coalLeftChildCoal[ci];
            if (leftCC >= 0) {
                scratchCoalDepth[leftCC] = d + 1;
                if (d + 1 > maxDepth) maxDepth = d + 1;
            }
            final int rightCC = coalRightChildCoal[ci];
            if (rightCC >= 0) {
                scratchCoalDepth[rightCC] = d + 1;
                if (d + 1 > maxDepth) maxDepth = d + 1;
            }
        }
        final int slabDepthCap = maxDepth + 2;

        int branchOpTotal = 0;
        int totalSlabBlocks = 0;
        for (int ci = 0; ci < coalCount; ci++) {
            final int lKb = coalLeftKbRec[ci];
            final int rKb = coalRightKbRec[ci];
            branchOpTotal += lKb + rKb;
            totalSlabBlocks += (lKb + slabOpsPerBlock - 1) / slabOpsPerBlock;
            totalSlabBlocks += (rKb + slabOpsPerBlock - 1) / slabOpsPerBlock;
        }
        final int flatTapeLen = branchOpTotal + branchCount;

        int forwardBufCount = 0;
        for (long w : forwardBufBitset) forwardBufCount += Long.bitCount(w);

        final int intervalCount = intervalCountSlot;

        final int operationCount = operationCountRecorded;
        final int reduceLen = 5 * operationCount;

        final int total = FIXED_HEADER_INTS
                + 8 * branchCount
                + (branchCount + 1)
                + (branchCount + 1)
                + flatTapeLen
                + 4 * coalCount
                + 4 * branchOpTotal
                + 3 * slabDepthCap
                + (intervalCount + 1)
                + branchOpTotal
                + forwardBufCount
                + 4 * totalSlabBlocks
                + reduceLen;

        if (payloadBuf == null || payloadBuf.length < total) {
            final int grow = (payloadBuf == null) ? total : total + (total >> 2);
            payloadBuf = new int[grow];
        }
        final int[] payload = payloadBuf;

        int p = FIXED_HEADER_INTS;
        final int offBrKb = p; p += branchCount;
        final int offBrKTop = p; p += branchCount;
        final int offBrTopBuf = p; p += branchCount;
        final int offBrBotBuf = p; p += branchCount;
        final int offBrDepth = p; p += branchCount;
        final int offBrOpFirst = p; p += branchCount;
        final int offBrTimeStart = p; p += branchCount + 1;
        final int offBrFirstBlock = p; p += branchCount;
        final int offBrBufOff = p; p += branchCount + 1;
        final int offFlatBufTape = p; p += flatTapeLen;
        final int offCoalDestBufs = p; p += coalCount;
        final int offCoalLeftAccBufs = p; p += coalCount;
        final int offCoalRightAccBufs = p; p += coalCount;
        final int offCoalIntervals = p; p += coalCount;
        final int offOpInBufOff = p; p += branchOpTotal;
        final int offOpKIn = p; p += branchOpTotal;
        final int offOpKAcc = p; p += branchOpTotal;
        final int offOpHasAcc = p; p += branchOpTotal;
        final int offBranchSlabStart = p; p += slabDepthCap;
        final int offCoalSlabStart = p; p += slabDepthCap;
        final int offSlabBlockStart = p; p += slabDepthCap;
        final int offBranchSlabList = p; p += branchCount;
        final int offIntStartCSR = p; p += intervalCount + 1;
        final int offIntListCSR = p; p += branchOpTotal;
        final int offForwardBufList = p; p += forwardBufCount;
        final int offSlabBlockBI = p; p += totalSlabBlocks;
        final int offSlabBlockCS = p; p += totalSlabBlocks;
        final int offSlabBlockCL = p; p += totalSlabBlocks;
        final int offSlabBlockCI = p; p += totalSlabBlocks;
        final int offOpReduce = p; p += reduceLen;
        assert p == total : "offset arithmetic disagrees with size calculation: p=" + p + " total=" + total;

        java.util.Arrays.fill(payload, offBranchSlabStart, offBranchSlabStart + slabDepthCap, 0);
        java.util.Arrays.fill(payload, offCoalSlabStart,offCoalSlabStart + slabDepthCap, 0);
        java.util.Arrays.fill(payload, offIntStartCSR,offIntStartCSR + intervalCount + 1, 0);

        java.util.Arrays.fill(payload, offOpHasAcc, offOpHasAcc + branchOpTotal, 1);

        payload[0] = VERSION;
        payload[1] = branchCount;
        payload[2] = coalCount;
        payload[3] = branchOpTotal;
        payload[4] = maxDepth;
        payload[5] = forwardBufCount;
        payload[6] = totalSlabBlocks;
        payload[7] = intervalCount;
        payload[8] = slabOpsPerBlock;
        payload[9] = slabDepthCap;
        payload[10] = operationCount;

        if (scratchStackCoal == null || scratchStackCoal.length < coalCount) {
            scratchStackCoal  = new int[coalCount + (coalCount >> 2)];
            scratchStackDepth = new int[coalCount + (coalCount >> 2)];
        }

        final int[] lKbRec  = coalLeftKbRec;
        final int[] rKbRec  = coalRightKbRec;
        final int[] lAccBuf = coalLeftAccBufRec;
        final int[] rAccBuf = coalRightAccBufRec;
        final int[] lBotBuf = coalLeftBotBufRec;
        final int[] rBotBuf = coalRightBotBufRec;
        final int[] lSnapS  = coalLeftSnapStart;
        final int[] rSnapS  = coalRightSnapStart;
        final int[] lSnapL  = coalLeftSnapLen;
        final int[] rSnapL  = coalRightSnapLen;
        final int[] lChildC = coalLeftChildCoal;
        final int[] rChildC = coalRightChildCoal;
        final int[] kTopRec = coalIntervalRec;
        final int[] snapAr  = coalSnapArena;

        int stackTop = 0;
        scratchStackCoal[stackTop]  = rootCoalIdx;
        scratchStackDepth[stackTop] = 0;
        stackTop++;
        payload[offCoalSlabStart + 1]++;

        int branchIdx = 0;
        int tapeOff   = 0;
        int opOff     = 0;
        int timeOff   = 0;

        while (stackTop > 0) {
            stackTop--;
            final int ci    = scratchStackCoal[stackTop];
            final int depth = scratchStackDepth[stackTop];
            final int kTop  = kTopRec[ci];
            final int depthHistOff = offCoalSlabStart + depth + 2;

            {
                final int Kb = lKbRec[ci];
                final int topBuf = lAccBuf[ci];
                final int botBuf = lBotBuf[ci];
                final int snapStart = lSnapS[ci];
                final int snapLen = lSnapL[ci];

                payload[offBrKb + branchIdx] = Kb;
                payload[offBrKTop + branchIdx] = kTop;
                payload[offBrTopBuf + branchIdx] = topBuf;
                payload[offBrBotBuf + branchIdx] = botBuf;
                payload[offBrDepth + branchIdx] = depth;
                payload[offBrOpFirst + branchIdx] = opOff;
                payload[offBrTimeStart + branchIdx] = timeOff;
                payload[offBrBufOff + branchIdx] = tapeOff;
                payload[offBranchSlabStart + depth + 1]++;

                final int tapeBase = offFlatBufTape + tapeOff;
                payload[tapeBase] = topBuf;
                if (snapLen > 0) {
                    System.arraycopy(snapAr, snapStart, payload, tapeBase + 1, snapLen);
                }
                payload[tapeBase + 1 + snapLen] = botBuf;

                final int opBase = offOpInBufOff + opOff;
                if (snapLen > 0) {
                    System.arraycopy(snapAr, snapStart, payload, opBase, snapLen);
                }
                payload[opBase + snapLen] = botBuf;

                final int lastN = Kb - 1;
                final int kKinBase = offOpKIn + opOff;
                final int kAccBase = offOpKAcc + opOff;
                for (int n = 0; n < lastN; n++) {
                    final int kk = kTop - n;
                    payload[kKinBase + n] = kk;
                    payload[kAccBase + n] = kk - 1;
                    if (kk >= 0 && kk < intervalCount) {
                        payload[offIntStartCSR + kk + 1]++;
                    }
                }

                {
                    final int kk = kTop - lastN;
                    payload[kKinBase + lastN] = kk;
                    payload[kAccBase + lastN] = 0;
                    payload[offOpHasAcc + opOff + lastN] = 0;
                    if (kk >= 0 && kk < intervalCount) {
                        payload[offIntStartCSR + kk + 1]++;
                    }
                }

                tapeOff += Kb + 1;
                opOff   += Kb;
                timeOff += Kb;
                branchIdx++;
                final int leftCC = lChildC[ci];
                if (leftCC >= 0) {
                    scratchStackCoal[stackTop]  = leftCC;
                    scratchStackDepth[stackTop] = depth + 1;
                    stackTop++;
                    payload[depthHistOff]++;
                }
            }

            {
                final int Kb = rKbRec[ci];
                final int topBuf = rAccBuf[ci];
                final int botBuf = rBotBuf[ci];
                final int snapStart = rSnapS[ci];
                final int snapLen = rSnapL[ci];

                payload[offBrKb + branchIdx] = Kb;
                payload[offBrKTop + branchIdx] = kTop;
                payload[offBrTopBuf + branchIdx] = topBuf;
                payload[offBrBotBuf + branchIdx] = botBuf;
                payload[offBrDepth + branchIdx] = depth;
                payload[offBrOpFirst + branchIdx] = opOff;
                payload[offBrTimeStart + branchIdx] = timeOff;
                payload[offBrBufOff + branchIdx] = tapeOff;
                payload[offBranchSlabStart + depth + 1]++;

                final int tapeBase = offFlatBufTape + tapeOff;
                payload[tapeBase] = topBuf;
                if (snapLen > 0) {
                    System.arraycopy(snapAr, snapStart, payload, tapeBase + 1, snapLen);
                }
                payload[tapeBase + 1 + snapLen] = botBuf;

                final int opBase = offOpInBufOff + opOff;
                if (snapLen > 0) {
                    System.arraycopy(snapAr, snapStart, payload, opBase, snapLen);
                }
                payload[opBase + snapLen] = botBuf;

                final int lastN = Kb - 1;
                final int kKinBase = offOpKIn + opOff;
                final int kAccBase = offOpKAcc + opOff;
                for (int n = 0; n < lastN; n++) {
                    final int kk = kTop - n;
                    payload[kKinBase + n] = kk;
                    payload[kAccBase + n] = kk - 1;
                    if (kk >= 0 && kk < intervalCount) {
                        payload[offIntStartCSR + kk + 1]++;
                    }
                }
                {
                    final int kk = kTop - lastN;
                    payload[kKinBase + lastN] = kk;
                    payload[kAccBase + lastN] = 0;
                    payload[offOpHasAcc + opOff + lastN] = 0;
                    if (kk >= 0 && kk < intervalCount) {
                        payload[offIntStartCSR + kk + 1]++;
                    }
                }

                tapeOff += Kb + 1;
                opOff   += Kb;
                timeOff += Kb;
                branchIdx++;
                final int rightCC = rChildC[ci];
                if (rightCC >= 0) {
                    scratchStackCoal[stackTop]  = rightCC;
                    scratchStackDepth[stackTop] = depth + 1;
                    stackTop++;
                    payload[depthHistOff]++;
                }
            }
        }

        payload[offBrBufOff + branchCount]    = tapeOff;
        payload[offBrTimeStart + branchCount] = timeOff;

        if (branchIdx != branchCount) {
            throw new IllegalStateException("DFS enumerated " + branchIdx + " branches, expected " + branchCount);
        }

        for (int d = 1; d < slabDepthCap; d++) {
            payload[offBranchSlabStart + d] += payload[offBranchSlabStart + d - 1];
        }
        if (scratchBranchCursor == null || scratchBranchCursor.length < slabDepthCap) {
            scratchBranchCursor = new int[slabDepthCap + (slabDepthCap >> 2)];
        }
        System.arraycopy(payload, offBranchSlabStart, scratchBranchCursor, 0, slabDepthCap);
        for (int bi = 0; bi < branchCount; bi++) {
            final int d = payload[offBrDepth + bi];
            payload[offBranchSlabList + scratchBranchCursor[d]++] = bi;
        }

        for (int d = 1; d < slabDepthCap; d++) {
            payload[offCoalSlabStart + d] += payload[offCoalSlabStart + d - 1];
        }
        if (scratchCoalCursor == null || scratchCoalCursor.length < slabDepthCap) {
            scratchCoalCursor = new int[slabDepthCap + (slabDepthCap >> 2)];
        }
        System.arraycopy(payload, offCoalSlabStart, scratchCoalCursor, 0, slabDepthCap);
        for (int ci = 0; ci < coalCount; ci++) {
            final int d = scratchCoalDepth[ci];
            if (d < 0) continue;
            final int pos = scratchCoalCursor[d]++;
            payload[offCoalDestBufs     + pos] = coalDestBufRec[ci];
            payload[offCoalLeftAccBufs  + pos] = coalLeftAccBufRec[ci];
            payload[offCoalRightAccBufs + pos] = coalRightAccBufRec[ci];
            payload[offCoalIntervals    + pos] = coalIntervalRec[ci];
        }


        {
            int g = 0;
            for (int d = 0; d <= maxDepth; d++) {
                payload[offSlabBlockStart + d] = g;
                final int brLo = payload[offBranchSlabStart + d];
                final int brHi = payload[offBranchSlabStart + d + 1];
                for (int j = brLo; j < brHi; j++) {
                    final int b  = payload[offBranchSlabList + j];
                    final int Kb = payload[offBrKb + b];
                    payload[offBrFirstBlock + b] = g;
                    final int numChunks = (Kb + slabOpsPerBlock - 1) / slabOpsPerBlock;
                    for (int c = 0; c < numChunks; c++) {
                        final int chunkStart = c * slabOpsPerBlock;
                        payload[offSlabBlockBI + g] = b;
                        payload[offSlabBlockCS + g] = chunkStart;
                        payload[offSlabBlockCL + g] = Math.min(slabOpsPerBlock, Kb - chunkStart);
                        payload[offSlabBlockCI + g] = c;
                        g++;
                    }
                }
            }
            payload[offSlabBlockStart + maxDepth + 1] = g;
            if (g != totalSlabBlocks) {
                throw new IllegalStateException("slab-block count mismatch: " + g + " vs " + totalSlabBlocks);
            }
        }


        if (intervalCount > 0 && branchOpTotal > 0) {
            for (int k = 1; k <= intervalCount; k++) {
                payload[offIntStartCSR + k] += payload[offIntStartCSR + k - 1];
            }
            if (scratchIntCursor == null || scratchIntCursor.length < intervalCount) {
                scratchIntCursor = new int[intervalCount + (intervalCount >> 2)];
            }
            System.arraycopy(payload, offIntStartCSR, scratchIntCursor, 0, intervalCount);
            for (int op = 0; op < branchOpTotal; op++) {
                final int kk = payload[offOpKIn + op];
                if (kk >= 0 && kk < intervalCount) {
                    payload[offIntListCSR + scratchIntCursor[kk]++] = op;
                }
            }
        }


        {
            final int divisor = Math.max(1, builderIndexOffsetPat);
            int idx = 0;
            for (int word = 0; word < forwardBufBitset.length; word++) {
                long bits = forwardBufBitset[word];
                while (bits != 0L) {
                    final int bit = Long.numberOfTrailingZeros(bits);
                    final int offset = (word << 6) | bit;
                    payload[offForwardBufList + idx++] = offset / divisor;
                    bits &= bits - 1L;
                }
            }
            if (idx != forwardBufCount) {
                throw new IllegalStateException("forward-buf count mismatch: " + idx + " vs " + forwardBufCount);
            }
        }

        if (reduceLen > 0) {
            System.arraycopy(opReduce, 0, payload, offOpReduce, reduceLen);
        }

        lastPackedLength = total;
        return payloadBuf;
    }


    public int getLastPackedLength() {
        return lastPackedLength;
    }


    private int[] emptyPackedDirect(int slabOpsPerBlock) {
        final int total = FIXED_HEADER_INTS + 1 + 1 + 2 + 2 + 2 + 1;
        if (payloadBuf == null || payloadBuf.length < total) {
            payloadBuf = new int[total];
        }
        java.util.Arrays.fill(payloadBuf, 0, total, 0);
        payloadBuf[0] = VERSION;
        payloadBuf[8] = slabOpsPerBlock;
        payloadBuf[9] = 2;
        payloadBuf[10] = 0;
        lastPackedLength = total;
        return payloadBuf;
    }
}
