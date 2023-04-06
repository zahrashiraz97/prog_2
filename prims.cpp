#include <iostream>
#include <vector>
#include <string>
#include <chrono>

using namespace std;
using namespace std::chrono;

typedef struct ge //graph edge
{
	long source; //vertex on one end
	long target; //vertex on the other end
	long weight;
} ge;

typedef struct gv //graph vertex
{
	long vertex; //number representing the vertex
	vector<ge> edges; //https://www.geeksforgeeks.org/cpp-vector-of-pointers/
} gv;

vector<gv> vertices;
vector<ge> mst_edges;

// lPair ==>  long Pair
typedef pair<long, long> lPair;

long num_comparisons = 0;

void increment_comparison()
{
	num_comparisons++;
}

template<typename T>
class min_heap {
  private:
  vector<T> data;
 
  public:
  // Implementing Priority Queue using inbuilt available vector in C++
  min_heap() {}
 
  // Element Inserting function
  void min_heap_enqueue(T item) {
    // item Insertion
    data.push_back(item);
    int ci = data.size() - 1;
 
    // Re-structure heap(min Heap) so that after
    // addition min element will lie on top of mh
    while (ci > 0) {
      int pi = (ci - 1) / 2;
      if (data[ci] >= data[pi])
        break;
      T tmp = data[ci];
      data[ci] = data[pi];
      data[pi] = tmp;
      ci = pi;
    }
  }
 
  T min_heap_dequeue() {
    // deleting top element of mh
    int li = data.size() - 1;
    T frontItem = data[0];
    data[0] = data[li];
    data.pop_back();
 
    --li;
    int pi = 0;
 
    // Re-structure heap(min Heap) so that after
    // deletion min element will lie on top of mh
    while (true) {
      int min_index = pi;
      int ci = pi * 2 + 1; //left child
      if (ci > li)
        break;
      if (ci <= li && data[ci] < data[min_index])
      {
        min_index = ci;
      }
      int rc = ci + 1; //right child
      if (rc <= li && data[rc] < data[min_index])
      {

        min_index = rc;
      }
      
      if (data[pi] < data[ci])
        break;
      if (min_index != pi)
      {
        T tmp = data[pi];
        data[pi] = data[min_index];
        data[min_index] = tmp;
        pi = min_index;
      }

    }
    return frontItem;
  }
 
  // Function which returns peek element
  T min_heap_check_top() {
    T frontItem = data[0];
    return frontItem;
  }

  int Count() {
    return data.size();
  }
};

// for each line of the form e v_1 v_2 wt do
// either …
// create a new edge object with two endpoints and weight
// add a reference/pointer to it to the ArrayList/vector
// … or, if you did 2(b)
// put the data into the edge array at position index and increment index
// put a reference to the the edge object on the list of incident edges of v_1 and v_2

void add_new_edge(long source, long target, long weight)
{
	ge edge;
	edge.source = source;
	edge.target = target;
	edge.weight = weight;
	// //cout<<"source = "<<source<<" target = "<<target<<" weight = "<<weight<<endl;
	int flag = 0; //source and target not updated
	vector<gv>::iterator ptr;
	bool flag_src = true;
	bool flag_target = true;
	for(ptr=vertices.begin(); ptr<vertices.end(); ptr++)
	{
		if (ptr->vertex == source)
		{
			(ptr->edges).push_back(edge);
			// vector<ge>::iterator ptr2;
			// for(ptr2=(ptr->edges).begin(); ptr2<(ptr->edges).end(); ptr2++)
			// {
			// 	//cout<<"Vertex (source) already in vector "<<ptr->vertex<<"edge "<<ptr2->source<<" "<<ptr2->target<<" "<<ptr2->weight<<endl;
			// }
			flag++;
			flag_src = false;
		}
		if (ptr->vertex == target)
		{
			(ptr->edges).push_back(edge);
			// vector<ge>::iterator ptr2;
			// for(ptr2=(ptr->edges).begin(); ptr2<(ptr->edges).end(); ptr2++)
			// {
			// 	//cout<<"Vertex (target) already in vector "<<ptr->vertex<<"edge "<<ptr2->source<<" "<<ptr2->target<<" "<<ptr2->weight<<endl;
			// }
			flag++;
			flag_target = false;
		}
		if (flag == 2)
		{
			break;
		}
	}
	// //cout<<"flag is "<<flag<<endl;
	if (flag == 0) //add both source and target
	{
		gv s;
		s.vertex = source;
		(s.edges).push_back(edge);
		vertices.push_back(s);
		// vector<ge>::iterator ptr;
		// for(ptr=(s.edges).begin(); ptr<(s.edges).end(); ptr++)
		// {
		// 	//cout<<"Adding Vertex (source) "<<s.vertex<<"edge "<<ptr->source<<" "<<ptr->target<<" "<<ptr->weight<<endl;
		// }

		gv t;
		t.vertex = target;
		(t.edges).push_back(edge);
		vertices.push_back(t);
		// vector<ge>::iterator ptr;
		// for(ptr=(t.edges).begin(); ptr<(t.edges).end(); ptr++)
		// {
		// 	//cout<<"Adding Vertex (target)"<<t.vertex<<"edge "<<ptr->source<<" "<<ptr->target<<" "<<ptr->weight<<endl;
		// }

	}
	else if (flag_target == true) //target to be added
	{
		gv t;
		t.vertex = target;
		(t.edges).push_back(edge);
		vertices.push_back(t);
		// vector<ge>::iterator ptr;
		// for(ptr=(t.edges).begin(); ptr<(t.edges).end(); ptr++)
		// {
		// 	//cout<<"Adding Vertex (target)"<<t.vertex<<"edge "<<ptr->source<<" "<<ptr->target<<" "<<ptr->weight<<endl;
		// }
	}
	else if (flag_src == true) //source to be added
	{
		gv s;
		s.vertex = source;
		(s.edges).push_back(edge);
		vertices.push_back(s);
		// vector<ge>::iterator ptr;
		// for(ptr=(s.edges).begin(); ptr<(s.edges).end(); ptr++)
		// {
		// 	//cout<<"Adding Vertex (source) "<<s.vertex<<"edge "<<ptr->source<<" "<<ptr->target<<" "<<ptr->weight<<endl;
		// }
	}
	else //nothing to be done
	{

	}

	// //cout<<"New updated made to vertices vector: "<<endl;
	// vector<gv>::iterator ptr2;
	// for(ptr2=vertices.begin(); ptr2<vertices.end(); ptr2++)
	// {
	// 	gv v = *ptr2;
	// 	//cout<<"Vertex is "<<v.vertex<<endl;
	// 	vector<ge>::iterator ptr;
	// 	for(ptr=(v.edges).begin(); ptr<(v.edges).end(); ptr++)
	// 	{
	// 		//cout<<"Vertex "<<v.vertex<<"edge "<<ptr->source<<" "<<ptr->target<<" "<<ptr->weight<<endl;
	// 	}
	// }


	// cout << endl;

}

# define INF 0x3f3f3f3f
 

// Prints shortest paths from src to all other vertices
long prim_mst()
{
	long num_vertices = vertices.size();
    min_heap< lPair> mh;
 
    long src = vertices[0].vertex; // Taking the first vertex as source
 
    // Create a vector for keys and initialize all
    // keys as infinite (INF)
    vector<long> key(num_vertices, INF);
 
    // To store parent array which in turn store MST
    vector<long> parent(num_vertices, -1);
 
    // To keep track of vertices included in MST
    vector<bool> inMST(num_vertices, false);
 
    // Insert source itself in priority queue and initialize
    // its key as 0.
    mh.min_heap_enqueue(make_pair(0, src));
    key[src] = 0;
 
    /* Looping till priority queue becomes empty */
    while (mh.Count() > 0)
    {
        // The first vertex in pair is the minimum key
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted key (key must be first item
        // in pair)
        int u = mh.min_heap_check_top().second;
        mh.min_heap_dequeue();
         
          //Different key values for same vertex may exist in the priority queue.
          //The one with the least key value is always processed first.
          //Therefore, ignore the rest.
          if(inMST[u] == true){
            continue;
        }
       
        inMST[u] = true;  // Include vertex in MST
        if (key[u] != 0) //do not include the source vertex "edge"
        {
        	ge e;
	        e.source = u;
	        e.target = parent[u];
	        e.weight = key[u];
	        mst_edges.push_back(e);
    	}

        //get all adjacent vertices of a vertex
        	vector<gv>::iterator ptr;
			for(ptr=vertices.begin(); ptr<vertices.end(); ptr++)
			{
				if (ptr->vertex == u)
				{
					vector<ge>::iterator ptr2;
					for(ptr2=(ptr->edges).begin(); ptr2<(ptr->edges).end(); ptr2++)
					{
						long v;
						if (ptr2->source == u)
						{
							v = ptr2->target;
						}
						else if (ptr2->target == u)
						{
							v = ptr2->source;
						}
            			long weight = ptr2->weight;
            			//cout<<"Considering vertex "<<v<<" and weight "<<weight<<endl;
						if (inMST[v] == false)
			            {
			           		increment_comparison();


			            	if (key[v] > weight)
			            	{
			                // Updating key of v
			                key[v] = weight;
			                //cout<<"New key["<<v<<"] = "<<key[v]<<endl;
			                mh.min_heap_enqueue(make_pair(key[v], v));
			                parent[v] = u;			            		
			            	}

			            }

					}


					break;
				}

			}
    }
 
    // Print edges of MST using parent array
    // for (int i = 1; i < num_vertices; ++i)
    //     printf("%ld - %d\n", parent[i], i);

    vector<ge>::iterator mst;
    long mst_weight = 0;
	for(mst=mst_edges.begin(); mst<mst_edges.end(); mst++)
	{
		ge e = *mst;
		cout<<"e "<<e.source<<" "<<e.target<<" "<<e.weight<<endl;
		mst_weight += e.weight;
	}

	return mst_weight;

}

 

int main()
{

	char c;
	long num_vertices = 0;
	long num_edges = 0;
	int flag = 0;
	int edge_counter = 0;
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
	    	cin>>num_vertices;
	    	cin>>num_edges;
	    	flag = 1;
	    	edge_counter = num_edges;
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
	    		add_new_edge(source, target, weight);
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
	auto start = high_resolution_clock::now();
	long mst_weight = prim_mst();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cerr<<"weight\t"<<mst_weight<<endl;
	cerr<<"runtime\t"<<duration.count()<<endl;
	cerr<<"comparisons\t"<<num_comparisons<<endl;
	return 0;
}






// Algorithm PrimJarnikMST(G) – edge heap with lazy deletion
// Q = new heap-based priority queue
// – entries are edges of G with keys = to their weights
// 	s = a vertex of G
// 	put all edges incident on s on Q
// 	mark s as in tree
// while Q is not empty do
//  		e = Q.removeMin()
// 		if both endpoints of e are in tree continue endif
// 		let u = endpoint of e not in tree
// 		add e to A and mark u as in tree
// 		for all edges f incident on u do
// 			let z = other endpoint of f (so f = uz)
// 			if z not in tree then
// 				put f on Q
// 			endif
// end for
// end while
// end PrimJarnikMST(G) – edge heap with lazy deletion
