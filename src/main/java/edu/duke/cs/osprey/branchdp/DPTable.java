/*
** This file is part of OSPREY 3.0
**
** OSPREY Protein Redesign Software Version 3.0
** Copyright (C) 2001-2018 Bruce Donald Lab, Duke University
**
** OSPREY is free software: you can redistribute it and/or modify
** it under the terms of the GNU General Public License version 2
** as published by the Free Software Foundation.
**
** You should have received a copy of the GNU General Public License
** along with OSPREY.  If not, see <http://www.gnu.org/licenses/>.
*/

package edu.duke.cs.osprey.branchdp;

import java.nio.DoubleBuffer;
import java.util.Arrays;

public interface DPTable {

    long size();

    double lower(long idx);

    double upper(long idx);

    void set(long idx, double lower, double upper);

    void fill(double lower, double upper);

    void copyFrom(long start, DoubleBuffer lowerSrc, DoubleBuffer upperSrc, int count);

    boolean isDenseArrayBacked();

    double[] lowerArrayUnsafe();

    double[] upperArrayUnsafe();

    long estimatedBytes();
}

final class DenseDPTable implements DPTable {

    private final double[] lower;
    private final double[] upper;

    DenseDPTable(long size) {
        if (size > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Dense DP table size exceeds Java array limits: " + size);
        }
        this.lower = new double[(int) size];
        this.upper = new double[(int) size];
    }

    @Override
    public long size() {
        return lower.length;
    }

    @Override
    public double lower(long idx) {
        return lower[(int) idx];
    }

    @Override
    public double upper(long idx) {
        return upper[(int) idx];
    }

    @Override
    public void set(long idx, double lower, double upper) {
        int i = (int) idx;
        this.lower[i] = lower;
        this.upper[i] = upper;
    }

    @Override
    public void fill(double lower, double upper) {
        Arrays.fill(this.lower, lower);
        Arrays.fill(this.upper, upper);
    }

    @Override
    public void copyFrom(long start, DoubleBuffer lowerSrc, DoubleBuffer upperSrc, int count) {
        if (count < 0 || start < 0 || start + (long)count > lower.length) {
            throw new IndexOutOfBoundsException("DP table copy [" + start + ","
                    + (start + (long)count) + ") outside [0," + lower.length + ")");
        }
        int offset = (int)start;
        lowerSrc.rewind();
        upperSrc.rewind();
        lowerSrc.get(lower, offset, count);
        upperSrc.get(upper, offset, count);
    }

    @Override
    public boolean isDenseArrayBacked() {
        return true;
    }

    @Override
    public double[] lowerArrayUnsafe() {
        return lower;
    }

    @Override
    public double[] upperArrayUnsafe() {
        return upper;
    }

    @Override
    public long estimatedBytes() {
        return 2L * (long) lower.length * (long) Double.BYTES;
    }
}

final class ShardedDPTable implements DPTable {

    private final long size;
    private final int shardSize;
    private final double[][] lowerShards;
    private final double[][] upperShards;

    ShardedDPTable(long size, int shardSize) {
        if (size < 0) {
            throw new IllegalArgumentException("Negative DP table size: " + size);
        }
        if (shardSize <= 0) {
            throw new IllegalArgumentException("Shard size must be positive: " + shardSize);
        }
        long shardCountLong = (size + shardSize - 1L) / shardSize;
        if (shardCountLong > Integer.MAX_VALUE) {
            throw new IllegalArgumentException("Too many DP shards: " + shardCountLong);
        }
        this.size = size;
        this.shardSize = shardSize;
        this.lowerShards = new double[(int) shardCountLong][];
        this.upperShards = new double[(int) shardCountLong][];

        for (int s = 0; s < lowerShards.length; s++) {
            long shardStart = (long) s * (long) shardSize;
            int len = (int) Math.min((long) shardSize, size - shardStart);
            lowerShards[s] = new double[len];
            upperShards[s] = new double[len];
        }
    }

    @Override
    public long size() {
        return size;
    }

    @Override
    public double lower(long idx) {
        int shard = shardIndex(idx);
        int offset = shardOffset(idx);
        return lowerShards[shard][offset];
    }

    @Override
    public double upper(long idx) {
        int shard = shardIndex(idx);
        int offset = shardOffset(idx);
        return upperShards[shard][offset];
    }

    @Override
    public void set(long idx, double lower, double upper) {
        int shard = shardIndex(idx);
        int offset = shardOffset(idx);
        lowerShards[shard][offset] = lower;
        upperShards[shard][offset] = upper;
    }

    @Override
    public void fill(double lower, double upper) {
        for (double[] shard : lowerShards) {
            Arrays.fill(shard, lower);
        }
        for (double[] shard : upperShards) {
            Arrays.fill(shard, upper);
        }
    }

    @Override
    public void copyFrom(long start, DoubleBuffer lowerSrc, DoubleBuffer upperSrc, int count) {
        if (count < 0 || start < 0 || start + (long)count > size) {
            throw new IndexOutOfBoundsException("DP table copy [" + start + ","
                    + (start + (long)count) + ") outside [0," + size + ")");
        }
        lowerSrc.rewind();
        upperSrc.rewind();
        long idx = start;
        int remaining = count;
        while (remaining > 0) {
            int shard = shardIndex(idx);
            int offset = shardOffset(idx);
            int len = Math.min(remaining, lowerShards[shard].length - offset);
            lowerSrc.get(lowerShards[shard], offset, len);
            upperSrc.get(upperShards[shard], offset, len);
            idx += len;
            remaining -= len;
        }
    }

    @Override
    public boolean isDenseArrayBacked() {
        return false;
    }

    @Override
    public double[] lowerArrayUnsafe() {
        throw new IllegalStateException("Sharded DP table is not backed by one lower[] array");
    }

    @Override
    public double[] upperArrayUnsafe() {
        throw new IllegalStateException("Sharded DP table is not backed by one upper[] array");
    }

    @Override
    public long estimatedBytes() {
        if (size > Long.MAX_VALUE / (2L * (long) Double.BYTES)) {
            return Long.MAX_VALUE;
        }
        return 2L * size * (long) Double.BYTES;
    }

    private int shardIndex(long idx) {
        checkIndex(idx);
        return (int) (idx / (long) shardSize);
    }

    private int shardOffset(long idx) {
        return (int) (idx % (long) shardSize);
    }

    private void checkIndex(long idx) {
        if (idx < 0 || idx >= size) {
            throw new IndexOutOfBoundsException("DP table index " + idx + " outside [0," + size + ")");
        }
    }
}
