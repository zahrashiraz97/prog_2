#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <string>


using namespace std;

int MAXN = 100000000;
long num_comparisons = 0; // Counter to keep track of the number of comparisons
long num_edge = 0; //Counter to keep track of the number of edges in the MST

void increment_comparison() // Function to increment the number of comparisons
{
    num_comparisons++;
}
struct Edge {
    long from, to, weight;
    Edge() {}
    Edge(int from, int to, int weight) : from(from), to(to), weight(weight) {}
    bool operator<(const Edge &rhs) const { //for comparing two edges based on their weights.
        return weight < rhs.weight;
    }
};

struct DisjointSet { // Disjoint set data structure
    int *parent = new int [MAXN];// Dynamically allocate an array for the parent pointers of the nodes.
    void init(int n) {
        for (int i = 0; i <= n; i++) {
            parent[i] = i;
        }
    }
    int find(int x) {// Find the parent of a node x.
        if (parent[x] == x) {
            return x;
        }
        return parent[x] = find(parent[x]);
    }
    void unite(int x, int y) {// Union operation for two nodes x and y.
        x = find(x);
        y = find(y);
        if (x != y) {
            parent[x] = y;
        }
    }
    bool same(int x, int y) {// Check if two nodes x and y belong to the same set.
        return find(x) == find(y);
    }
};

struct BinaryHeap {// Binary heap data structure
    vector<Edge> heap;
    void push(const Edge &e) {// Insert a new element into the heap.
        heap.push_back(e);
        int i = heap.size() - 1;
        while (i > 0) {
            int p = (i - 1) / 2;// Compute the parent index.
            increment_comparison();
            if (heap[p] < heap[i]) {// If the parent is less than the child, stop the loop.
                break;
            }
            swap(heap[p], heap[i]);// Swap the parent and child nodes.
            
            i = p;
        }
    }
    Edge removeMin() {// Remove the minimum element from the heap.
        Edge ret = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        int i = 0;
        while (i * 2 + 1 < heap.size()) {// Traverse the heap
            int left = i * 2 + 1, right = i * 2 + 2;
            if (right < heap.size()) {
                increment_comparison();
                if (heap[right] < heap[left])// compare the left and right children and choose the smaller one.
                left = right;
            }
            increment_comparison();
            if (heap[i] < heap[left]) {
                break;
            }
            swap(heap[i], heap[left]);
            // comparisons++;
            i = left;
        }
        return ret;// Return the minimum element that was stored at the root of the heap.
    }
    bool empty() {
        return heap.empty();
    }
};

vector<Edge> kruskal(int n, vector<Edge> edges) {
    // sort(edges.begin(), edges.end());


    DisjointSet dsu;
    dsu.init(n);

    BinaryHeap heap;
    for (auto &e : edges) {
        heap.push(e);
    }

    vector<Edge> ret;
   
    while (!heap.empty() && ret.size() < n - 1) {
        Edge e = heap.removeMin();
        // comparisons++;
        if (!dsu.same(e.from, e.to)) {
            dsu.unite(e.from, e.to);
            ret.push_back(e);
        }
    }
    return ret;
}

long get_total_weight(vector<Edge> &edges) {
    long total_weight = 0;
    for (auto &e : edges) {
        total_weight += e.weight;
    }
    return total_weight;
}


int main() {

    char c;
    long n = 0;
    long m = 0;
    int flag = 0;
    int edge_counter = 0;
     vector<Edge> edges;
     int i = 0;
    while((c = cin.peek())) {
        if (c == 'c')
        {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');// Ignore comment lines
        }
        else if (c == 'n')
        {
            cin.ignore(numeric_limits<streamsize>::max(), '\n'); //Ignore lines starting with 'n'
        }
        else if ((c == 'g') && (flag != 1))// // Start getting input
        {
            cin>>c;
            cin>>n;
            cin>>m;
            flag = 1;
            edge_counter = m;
            edges.resize(m);
        }
        else if ((c == 'e') && (flag == 1))
        {
            if (edge_counter > 0)
            {
                long source;
                long target;
                long weight;
                cin>>c;
                cin>>source;
                cin>>target;
                cin>>weight;
                edges[i] = Edge(source, target, weight);
                i++;
                edge_counter--;
                if (edge_counter == 0)// All edges have been read in, so stop processing
                {
                    break;
                }
            }
            else
            {
                cin.ignore(numeric_limits<streamsize>::max(), '\n');// Ignore extra lines after reading all edges

            }

        }
        else
        {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');// Ignore any other input lines
        }

    }


    auto start_time = chrono::high_resolution_clock::now();

    

    vector<Edge> mst = kruskal(n, edges);
    cout << "g " << n << " " << mst.size() << endl;

    for (auto &e : mst) {

        // cout << "e " << e.from << " " << e.to << " " << e.weight << endl;
        num_edge ++ ;


    }

    if (num_edge!= n -1 ) {
        cerr<<"Error! Disconnected graph provided"<<endl;
    }

    auto end_time = chrono::high_resolution_clock::now();

    // output statistics to standard error

    cerr << "weight\t" << get_total_weight(mst) << endl; //show total weight
    cerr << "runtime\t" << chrono::duration_cast<chrono::microseconds>(end_time - start_time).count()<<endl;
     //double seconds = static_cast<double>((end_time - start_time).count()) / 1000000 ;
    // cerr<<"seconds\t"<<seconds<<endl;
    cerr<<"comparisons\t"<<num_comparisons<<endl;// show number of comparisons
    // cerr<<"edges\t"<<num_edge<<endl;





    return 0;

}