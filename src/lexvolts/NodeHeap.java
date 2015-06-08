package lexvolts;

/**
 * Created by spielman on 8/31/14.
 * A simple implementation of a min-heap, with real number keys and small integer values (nodes)
 */
public class NodeHeap {

    public double[] keys;
    public int[] heap;
    public int[] index; // where in the heap the ith node is
    public int nitems;

    /* **********************
        FOR COMPAT WITH NODEQUEUE
     */
    public double[] values;

    // create a heap with max maxitems in it
    public NodeHeap(int maxitems) {
        keys = new double[maxitems];
        heap = new int[maxitems];
        index = new int[maxitems];
        nitems = 0;


        values = keys;

        // use -1 to indicate not in the heap

        for (int i = 0; i < maxitems; i++) {
            index[i] = -1;
        }

    }

    // will check if was already in the heap
    // if it was, allow a decreasekey
    public boolean add(int node, double key) {
        if (index[node] >= 0) {
            if (key < keys[node]) {
                changeKey(node, key);
                return true;
            } else
                return false;
        }
        else {
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

            int parentPos = (pos-1) / 2;
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

        int leftPos = pos*2 + 1;

        boolean moved = true;

        while ((leftPos < nitems) && moved) {

            moved = false;

            int rightPos = pos*2 + 2;

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
                leftPos = pos*2 + 1;
                moved = true;
            }

        }

    }

    // this was just written for testing
    public double[] heapSort(double[] toSort) {
        int n = toSort.length;
        for (int i = 0; i < n; i++)
            add(i,toSort[i]);

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

    // This procedure assumes that the key is decreasing
    // WARNING: This procedure performs no checks
    public void decreaseKey(int node, double key) {
        keys[node] = key;
        percUp(node);
    }

    public boolean hasMore() {
        return (nitems > 0);
    }


    /* *******************************************
    these routines are for compat with NodeQueue
    ********************************************* */
    public int popMinNode() {
        return popMin();
    }

    public int peekMinNode() {
        return heap[0];
    }

    public void changeValue(int node, double value) {
        changeKey(node, value);
    }

    public boolean inQueue(int node) {
        return (index[node] >= 0);

    }
}

