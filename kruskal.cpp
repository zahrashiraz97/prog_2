#include <iostream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <sstream>
#include <string>


using namespace std;

const int MAXN = 100005;
long num_comparisons = 0;

void increment_comparison()
{
    num_comparisons++;
}
struct Edge {
    long from, to, weight;
    Edge() {}
    Edge(int from, int to, int weight) : from(from), to(to), weight(weight) {}
    bool operator<(const Edge &rhs) const {
        return weight < rhs.weight;
    }
};

struct DisjointSet {
    int parent[MAXN];
    void init(int n) {
        for (int i = 0; i <= n; i++) {
            parent[i] = i;
        }
    }
    int find(int x) {
        if (parent[x] == x) {
            return x;
        }
        return parent[x] = find(parent[x]);
    }
    void unite(int x, int y) {
        x = find(x);
        y = find(y);
        if (x != y) {
            parent[x] = y;
        }
    }
    bool same(int x, int y) {
        return find(x) == find(y);
    }
};

struct BinaryHeap {
    vector<Edge> heap;
    void push(const Edge &e) {
        heap.push_back(e);
        int i = heap.size() - 1;
        while (i > 0) {
            int p = (i - 1) / 2;
            increment_comparison();
            if (heap[p] < heap[i]) {
                break;
            }
            swap(heap[p], heap[i]);
            
            i = p;
        }
    }
    Edge removeMin() {
        Edge ret = heap[0];
        heap[0] = heap.back();
        heap.pop_back();
        int i = 0;
        while (i * 2 + 1 < heap.size()) {
            int left = i * 2 + 1, right = i * 2 + 2;
            if (right < heap.size()) {
                increment_comparison();
                if (heap[right] < heap[left])
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
        return ret;
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
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        else if (c == 'n')
        {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }
        else if ((c == 'g') && (flag != 1))
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
                if (edge_counter == 0)
                {
                    break;
                }
            }
            else
            {
                cin.ignore(numeric_limits<streamsize>::max(), '\n');

            }

        }
        else
        {
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
        }

    }


    auto start_time = chrono::high_resolution_clock::now();

    vector<Edge> mst = kruskal(n, edges);

    auto end_time = chrono::high_resolution_clock::now();



    // output results

    cout << "g " << n << " " << mst.size() << endl;

    for (auto &e : mst) {

        cout << "e " << e.from << " " << e.to << " " << e.weight << endl;

    }



    // output statistics to standard error

    cerr << "weight\t" << get_total_weight(mst) << endl;
    cerr << "runtime\t" << chrono::duration_cast<chrono::seconds>(end_time - start_time).count()<<endl;
    cerr<<"comparisons\t"<<num_comparisons<<endl;




    return 0;

}




