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
	struct gv *source_ptr;
	struct gv *target_ptr;

} ge;

typedef struct gv //graph vertex
{
	long vertex; //number representing the vertex
	vector<ge*> edges;
} gv;

vector<gv*> vertices;
vector<ge*> g_edges;
vector<ge*> mst_edges;

// ePair ==>  pair of weight of edge and ptr to edge object
typedef pair<long, ge*> ePair;

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
       increment_comparison();

      if (data[ci] >= data[pi])
        {
        	break;
        }
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
      if (ci <= li)
      {
      	increment_comparison();
      	if (data[ci] < data[min_index])
      	{
      		 min_index = ci;

      	}
      }
      int rc = ci + 1; //right child
      if (rc <= li)
      {

      	increment_comparison();
      	if  (data[rc] < data[min_index])
      	{
      		        min_index = rc;
      	}

      }
      
      increment_comparison();
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


void add_new_edge(long source, long target, long weight)
{
	ge *new_edge = new ge; // creating new edge 
	new_edge->source = source;
	new_edge->target = target;
	new_edge->weight = weight;

	ge *edge = new_edge;
	int flag = 0; //source and target not updated
	vector<gv*>::iterator v_ptr;
	bool flag_src = true;
	bool flag_target = true;
	for(v_ptr=vertices.begin(); v_ptr<vertices.end(); v_ptr++)
	{
		gv *ptr = *v_ptr;
		if (ptr->vertex == source)
		{
			new_edge->source_ptr = ptr;
			flag++;
			flag_src = false;
		}
		if (ptr->vertex == target)
		{
			new_edge->target_ptr = ptr;
			flag++;
			flag_target = false;
		}
		if (flag == 2)
		{
			break;
		}
	}
	if (flag == 0) //add both source and target
	{
		gv *s = new gv; // creating new vertex
		s->vertex = source;
		new_edge->source_ptr = s;
		vertices.push_back(s);

		gv *t = new gv;
		t->vertex = target;
		new_edge->target_ptr = t;
		vertices.push_back(t);

	}
	else if (flag_target == true) //target to be added
	{
		gv *t = new gv;
		t->vertex = target;
		new_edge->target_ptr = t;
		vertices.push_back(t);
	}
	else if (flag_src == true) //source to be added
	{

		gv *s = new gv;
		s->vertex = source;
		new_edge->source_ptr = s;
		vertices.push_back(s);

	}
	else //nothing to be done
	{

	}
	gv *s_ptr = new_edge->source_ptr;
	gv *t_ptr = new_edge->target_ptr;

	((s_ptr)->edges).push_back(new_edge);
	((t_ptr)->edges).push_back(new_edge);
 	g_edges.push_back(new_edge);

}


long prim_jarnik_mst()
{
	long num_vertices = vertices.size();
    min_heap<ePair> mh;

 
    // To keep track of vertices included in MST
    vector<bool> inMST(num_vertices, false);

      vector<gv*>::iterator ptr;
			for(ptr=vertices.begin(); ptr<vertices.end(); ptr++)
			{
				gv *v_ptr = *ptr;
				gv v = *v_ptr;

				vector<ge*>::iterator ptr2;
				for(ptr2=(v.edges).begin(); ptr2<(v.edges).end(); ptr2++)
				{
					ge *ptr_to_edge = *ptr2;
					ge e = *ptr_to_edge;
					mh.min_heap_enqueue(make_pair(e.weight, ptr_to_edge));
				}		
				inMST[v.vertex] = true;
				break; //only for the first vertex
			}

			int edges_check = 0;
 
    /* Looping till priority queue becomes empty */
    while (mh.Count() > 0)
    {
     	ge *ptr_to_edge = mh.min_heap_check_top().second; 
    	ge e = *ptr_to_edge;

    	mh.min_heap_dequeue();

    	if (inMST[e.source] && inMST[e.target])
    	{
    		continue;
    	}
    	long v;
    	gv *other_vertex = NULL;
    	if (inMST[e.source])
    	{
    		v = e.target;
    		other_vertex = e.target_ptr;
    	}
    	else if (inMST[e.target])
    	{
    		v = e.source;
    		other_vertex = e.source_ptr;
    	}
    	mst_edges.push_back(ptr_to_edge);
    	inMST[v] = true;

			vector<ge*>::iterator ptr2;
			for(ptr2=(other_vertex->edges).begin(); ptr2<(other_vertex->edges).end(); ptr2++)
			{
				long z;
				ge *ptr_to_edge = *ptr2;
				ge e = *ptr_to_edge;
				ge *edge_ptr = (ge *) &ptr2;
				if (e.source == v)
				{
					z = e.target;
				}
				else if (e.target == v)
				{
					z = e.source;
				}
				if (inMST[z] == false)
				{
					mh.min_heap_enqueue(make_pair(e.weight, ptr_to_edge));
				}
			}		

		}


    vector<ge*>::iterator mst;
    long mst_weight = 0;
    cout<<"g "<<num_vertices<<" "<<num_vertices-1<<endl;
	for(mst=mst_edges.begin(); mst<mst_edges.end(); mst++)
	{
		ge *ptr_to_edge = *mst;
		ge e = *ptr_to_edge;
		cout<<"e "<<e.source<<" "<<e.target<<" "<<e.weight<<endl;
		mst_weight += e.weight;
		edges_check++;
	}
	if (edges_check != num_vertices-1)
	{
		cout<<"Error! Disconnected graph provided"<<endl;
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
	long mst_weight = prim_jarnik_mst();
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cerr<<"weight\t"<<mst_weight<<endl;
	cerr<<"runtime\t"<<duration.count()<<endl;
	cerr<<"comparisons\t"<<num_comparisons<<endl;
	return 0;
}

