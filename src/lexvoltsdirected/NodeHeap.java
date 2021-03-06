package lexvoltsdirected;

/**
 * Created by sushant on 1/31/15.
 */
public class NodeHeap {

    private double[] keys;
    private int[] heap;
    private int[] index; // where in the heap the ith node is
    private int nitems;

    private int usageCount;
    private int[] lastSeen;
    private int maxSize;

    // create a heap with max maxitems in it
    public NodeHeap(int maxitems) {
        this.maxSize = maxitems;
        keys = new double[maxitems];
        heap = new int[maxitems];
        index = new int[maxitems];
        lastSeen = new int[maxitems];
        nitems = 0;

        // use -1 to indicate not in the heap
        for (int i = 0; i < maxitems; i++) {
            index[i] = -1;
            keys[i] = Double.POSITIVE_INFINITY;
            lastSeen[i] = -1;
        }
        usageCount = 0;
    }

    public void resetHeap() {
        nitems = 0;
        usageCount++;
        if (usageCount < 0) {
            // Time to reset all counters
            System.out.println("Heap count overflow. Linear time reset. Don't worry, not an error!");
            for (int i = 0; i < maxSize; i++)
                lastSeen[i] = -1;
            usageCount = 0;
        }
    }

    // Returns a boolean indicating if this vertex has been finished in this iteration
    public boolean finishedThisTime(int node) {
        return (lastSeen[node] == usageCount);
    }

    // will check if was already in the heap
    // if it was, allow a decreasekey
    // The returned boolean value indicates whether there was a reduction in the node's key. If it is added for the first time, we return true.
    public boolean add(int node, double key) {
        if (index[node] >= 0) {
            if (key < keys[node]) {
                changeKey(node, key);
                return true;
            } else
                return false;
        } else {
            keys[node] = key;
            heap[nitems] = node;
            index[node] = nitems;

            nitems++;
            percUp(node);
            return true;
        }
    }

    private void percUp(int node) {

        int pos = index[node];
        boolean moved = true;

        while ((pos > 0) && (moved)) {
            double key = keys[node];

            int parentPos = (pos - 1) / 2;
            int parentNode = heap[parentPos];
            double parentKey = keys[parentNode];

            moved = false;

            if (parentKey > key) {
                heap[parentPos] = node;
                heap[pos] = parentNode;
                index[node] = parentPos;
                index[parentNode] = pos;
                pos = parentPos;
                moved = true;
            }

        }

    }

    public int popMin() {
        if (nitems == 0)
            throw new Error("Heap is empty.");

        nitems--;
        int minNode = heap[0];

        index[minNode] = -1;
        keys[minNode] = Double.POSITIVE_INFINITY;
        lastSeen[minNode] = usageCount;

        if (nitems > 0) {

            // get the last node, put it at the root, and push it down
            int pos = nitems;
            int node = heap[pos];

            heap[0] = node;
            index[node] = 0;

            siftDown(node);
        }


        return minNode;
    }

    private void siftDown(int node) {

        int pos = index[node];
        double key = keys[node];

        int leftPos = pos * 2 + 1;

        boolean moved = true;

        while ((leftPos < nitems) && moved) {

            moved = false;

            int rightPos = pos * 2 + 2;

            int childPos;
            int childNode;
            double childKey;
            if (rightPos == nitems) {
                // is no right child, so swap with left
                childPos = leftPos;
                childNode = heap[childPos];
                childKey = keys[childNode];

            } else {
                int leftNode = heap[leftPos];
                double leftKey = keys[leftNode];

                int rightNode = heap[rightPos];
                double rightKey = keys[rightNode];


                if (leftKey < rightKey) {
                    childPos = leftPos;
                    childNode = leftNode;
                    childKey = leftKey;
                } else {

                    childPos = rightPos;
                    childNode = rightNode;
                    childKey = rightKey;
                }
            }


            if (childKey < key) {
                heap[childPos] = node;
                heap[pos] = childNode;
                index[node] = childPos;
                index[childNode] = pos;

                pos = childPos;
                leftPos = pos * 2 + 1;
                moved = true;
            }

        }

    }

    // this was just written for testing
    public double[] heapSort(double[] toSort) {
        int n = toSort.length;
        for (int i = 0; i < n; i++)
            add(i, toSort[i]);

        double out[] = new double[n];

        for (int i = 0; i < n; i++)
            out[i] = keys[popMin()];

        return out;
    }

    public void changeKey(int node, double key) {
        double oldKey = keys[node];

        keys[node] = key;

        if (key < oldKey)
            percUp(node);
        else
            siftDown(node);
    }

    public boolean hasMore() {
        return (nitems > 0);
    }

    public double minKey() {
        assert (nitems > 0);
        return keys[heap[0]];
    }
}

