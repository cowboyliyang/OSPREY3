package edu.duke.cs.osprey.ematrix;

import edu.duke.cs.osprey.confspace.SimpleConfSpace;

import java.util.*;

/**
 * Branch Decomposition for protein structures
 *
 * Decomposes the interaction graph of a protein into a binary tree structure,
 * where each node represents a subset of positions and edges represent separators.
 *
 * This enables:
 * 1. Subtree DOF caching (reuse minimization results for common subtrees)
 * 2. Efficient divide-and-conquer minimization
 * 3. Reduced redundant computation in A* search
 *
 * Algorithm: Greedy balanced partitioning
 * - Splits positions into roughly equal halves
 * - Minimizes separator size (positions connecting the two halves)
 * - Recursively builds the tree
 */
public class BranchDecomposition {

    public final TreeNode root;
    public final int branchWidth;  // Maximum separator size across all nodes
    private final SimpleConfSpace confSpace;

    /**
     * Tree node in the branch decomposition
     */
    public static class TreeNode {
        public final Set<Integer> positions;      // Positions in this subtree
        public final Set<Integer> separator;      // Separator positions (connecting left/right)
        public final TreeNode leftChild;
        public final TreeNode rightChild;
        public final boolean isLeaf;

        // For leaf nodes
        public TreeNode(Set<Integer> positions) {
            this.positions = new HashSet<>(positions);
            this.separator = new HashSet<>();
            this.leftChild = null;
            this.rightChild = null;
            this.isLeaf = true;
        }

        // For internal nodes
        public TreeNode(Set<Integer> positions, Set<Integer> separator,
                       TreeNode leftChild, TreeNode rightChild) {
            this.positions = new HashSet<>(positions);
            this.separator = new HashSet<>(separator);
            this.leftChild = leftChild;
            this.rightChild = rightChild;
            this.isLeaf = false;
        }

        public int getSeparatorSize() {
            return separator.size();
        }

        @Override
        public String toString() {
            if (isLeaf) {
                return "Leaf" + positions;
            } else {
                return "Node" + positions + " sep=" + separator;
            }
        }
    }

    /**
     * Build branch decomposition for a conformation space
     */
    public BranchDecomposition(SimpleConfSpace confSpace) {
        this.confSpace = confSpace;

        // Get all positions
        Set<Integer> allPositions = new HashSet<>();
        for (int i = 0; i < confSpace.positions.size(); i++) {
            allPositions.add(i);
        }

        // Build interaction graph (simplified - assume all positions interact)
        // In a real implementation, would analyze energy matrix for actual interactions
        Map<Integer, Set<Integer>> graph = buildInteractionGraph(allPositions);

        // Build tree
        this.root = buildTree(allPositions, graph);

        // Compute branch width
        this.branchWidth = computeBranchWidth(root);
    }

    /**
     * Build interaction graph (simplified version)
     * Returns adjacency list: position -> set of neighboring positions
     */
    private Map<Integer, Set<Integer>> buildInteractionGraph(Set<Integer> positions) {
        Map<Integer, Set<Integer>> graph = new HashMap<>();

        // Simplified: assume positions interact with their neighbors
        // In reality, should analyze pairwise energy matrix
        List<Integer> sortedPos = new ArrayList<>(positions);
        Collections.sort(sortedPos);

        for (int i = 0; i < sortedPos.size(); i++) {
            int pos = sortedPos.get(i);
            Set<Integer> neighbors = new HashSet<>();

            // Connect to adjacent positions (simplified interaction model)
            if (i > 0) {
                neighbors.add(sortedPos.get(i - 1));
            }
            if (i < sortedPos.size() - 1) {
                neighbors.add(sortedPos.get(i + 1));
            }

            graph.put(pos, neighbors);
        }

        return graph;
    }

    /**
     * Recursively build the tree using greedy balanced partitioning
     */
    private TreeNode buildTree(Set<Integer> positions, Map<Integer, Set<Integer>> graph) {
        // Base case: 3 or fewer positions - stop splitting
        if (positions.size() <= 3) {
            return new TreeNode(positions);
        }

        // Recursive case: partition into two roughly equal parts
        Partition partition = greedyPartition(positions, graph);

        // Build subtrees
        TreeNode leftNode = buildTree(partition.left, graph);
        TreeNode rightNode = buildTree(partition.right, graph);

        return new TreeNode(positions, partition.separator, leftNode, rightNode);
    }

    /**
     * Greedy balanced partitioning
     * Splits positions into two sets to minimize separator size
     */
    private static class Partition {
        Set<Integer> left;
        Set<Integer> right;
        Set<Integer> separator;

        Partition(Set<Integer> left, Set<Integer> right, Set<Integer> separator) {
            this.left = left;
            this.right = right;
            this.separator = separator;
        }
    }

    private Partition greedyPartition(Set<Integer> positions, Map<Integer, Set<Integer>> graph) {
        List<Integer> posList = new ArrayList<>(positions);
        Collections.sort(posList);

        // Simple strategy: split in half
        int mid = posList.size() / 2;

        Set<Integer> left = new HashSet<>();
        Set<Integer> right = new HashSet<>();

        for (int i = 0; i < posList.size(); i++) {
            if (i < mid) {
                left.add(posList.get(i));
            } else {
                right.add(posList.get(i));
            }
        }

        // Find separator (positions with edges crossing the partition)
        Set<Integer> separator = findSeparator(left, right, graph);

        return new Partition(left, right, separator);
    }

    /**
     * Find separator: positions with edges connecting left and right partitions
     */
    private Set<Integer> findSeparator(Set<Integer> left, Set<Integer> right,
                                       Map<Integer, Set<Integer>> graph) {
        Set<Integer> separator = new HashSet<>();

        // Check positions in left partition
        for (int pos : left) {
            Set<Integer> neighbors = graph.getOrDefault(pos, new HashSet<>());
            for (int neighbor : neighbors) {
                if (right.contains(neighbor)) {
                    separator.add(pos);
                    break;
                }
            }
        }

        // Check positions in right partition
        for (int pos : right) {
            Set<Integer> neighbors = graph.getOrDefault(pos, new HashSet<>());
            for (int neighbor : neighbors) {
                if (left.contains(neighbor)) {
                    separator.add(pos);
                    break;
                }
            }
        }

        return separator;
    }

    /**
     * Compute branch width (maximum separator size)
     */
    private int computeBranchWidth(TreeNode node) {
        if (node.isLeaf) {
            return 0;
        }

        int maxWidth = node.getSeparatorSize();

        if (node.leftChild != null) {
            maxWidth = Math.max(maxWidth, computeBranchWidth(node.leftChild));
        }

        if (node.rightChild != null) {
            maxWidth = Math.max(maxWidth, computeBranchWidth(node.rightChild));
        }

        return maxWidth;
    }

    /**
     * Check if this decomposition is suitable for DOF caching
     * (branch width should be small)
     */
    public boolean isSuitable(int maxBranchWidth) {
        return this.branchWidth <= maxBranchWidth;
    }

    /**
     * Print tree structure for debugging
     */
    public void printTree() {
        printTree(root, 0);
    }

    private void printTree(TreeNode node, int depth) {
        String indent = "  ".repeat(depth);
        System.out.println(indent + node);

        if (!node.isLeaf) {
            if (node.leftChild != null) {
                printTree(node.leftChild, depth + 1);
            }
            if (node.rightChild != null) {
                printTree(node.rightChild, depth + 1);
            }
        }
    }

    /**
     * Get the number of positions in this branch decomposition
     */
    public int getPositionCount() {
        return confSpace.positions.size();
    }

    /**
     * Get statistics about this decomposition
     */
    public String getStats() {
        int[] counts = new int[3];  // [totalNodes, leafNodes, internalNodes]
        collectStats(root, counts);

        return String.format(
            "BranchDecomposition[nodes=%d, leaves=%d, internal=%d, branchWidth=%d]",
            counts[0], counts[1], counts[2], branchWidth
        );
    }

    private void collectStats(TreeNode node, int[] counts) {
        counts[0]++;  // total nodes

        if (node.isLeaf) {
            counts[1]++;  // leaf nodes
        } else {
            counts[2]++;  // internal nodes
            if (node.leftChild != null) {
                collectStats(node.leftChild, counts);
            }
            if (node.rightChild != null) {
                collectStats(node.rightChild, counts);
            }
        }
    }
}
