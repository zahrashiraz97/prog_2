/**
 * @file delauney-triangulation.c
 * @brief Creates a Delaunay triangulation on a random set of points; edges
 * are directed so as to form a dag.
 *
 * @author Tang Yu
 * @date Summer, 2010
 *
 * $Id: delauney-triangulation.c 305 2014-07-23 14:05:04Z mfms $
 */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <stdbool.h>
#include <assert.h>
#include <limits.h>
#include "random.h"

#define Vector(p1, p2, u, v) (u = p2->x - p1->x, v = p2->y - p1->y)
#define Cross_product_2v(u1, v1, u2, v2) (u1 * v2 - v1 * u2)
#define Cross_product_3p(p1, p2, p3) ((p2->x - p1->x) * (p3->y - p1->y) - (p2->y - p1->y) * (p3->x - p1->x))
#define Dot_product_2v(u1, v1, u2, v2) (u1 * u2 + v1 * v2)

#define Org(e) ((e)->org)
#define Dest(e) ((e)->dest)
#define Onext(e) ((e)->onext)
#define Oprev(e) ((e)->oprev)
#define Dnext(e) ((e)->dnext)
#define Dprev(e) ((e)->dprev)

#define Other_point(e, p) ((e)->org == p ? (e)->dest : (e)->org)
#define Next(e, p) ((e)->org == p ? (e)->onext : (e)->dnext)
#define Prev(e, p) ((e)->org == p ? (e)->oprev : (e)->dprev)

#define Visited(p) ((p)->f)

#define Identical_refs(e1, e2) (e1 == e2)

/* guarantee that no edge has cost more than UINT_MAX */
#define RANGE (UINT_MAX / 2)

typedef enum {right, left} side;

typedef double Real;

/** @todo names of typedef'd structs should differ from original structs */
typedef struct point point;
typedef struct edge edge;
typedef struct EDGE EDGE;

struct point
{
    Real x, y;
    edge *entry_pt;
};

struct edge
{
    point *org, *dest;
    edge *onext, *oprev, *dnext, *dprev;
};


struct EDGE
{
    int u, v;
    double w;
};

static point *p_array;

EDGE * E;
EDGE * mst;
double * dist;
int * p, * d;

/**
 * true if a vertex has a predecessor -- used to ensure that the outer face is
 * also a triangle; doesn't completely solve the problem, though; consider a
 * simple square divided into two triangles; however, it's probably good
 * enough for use with crossing minimization as it guarantees a single source.
 */
bool * has_predecessor;

int n, en, k;
unsigned long seed = 0;

double Distance(point * p1, point * p2) {
  double x_diff = p1->x - p2->x;
  double y_diff = p1->y - p2->y;
  return sqrt(x_diff * x_diff + y_diff * y_diff);
}

int ecmp(const void * a, const void * b)
{
    return ( ((EDGE *)a)->w - ((EDGE *)b)->w );
}

int pcmp(const void * a, const void * b)
{
    point * pa = (point *)a;
    point * pb = (point *)b;
    if ((pa->x < pb->x) || (pa->x == pb->x && pa->y < pb->y)) return -1;
    else if ((pa->x == pb->x) && (pa->y == pb->y)) return 0;
    else return 1;
}

void output_graph(const char * output_file_name)
{
  FILE * output_file;
  if ( output_file_name != NULL ) {
    /**
     * @todo print useful message if file can't be written
     */
    output_file = fopen(output_file_name, "w");
    assert( output_file != NULL );
  }
  else {
    output_file = stdout;
  }
  fprintf( output_file, "c Delaunay triangulation: %d vertices, %d edges, seed = %lu",
           n, en, seed );
  time_t raw_time = time(NULL);
  struct tm * current_time = gmtime(&raw_time);
  fprintf( output_file, ", %s", asctime(current_time));
    fprintf(output_file, "g %d %d\n", n, en);
    // print position information for vertices in format
    //      n id x y
    // the id's are 1-based rather than 0-based, as they are for the edge
    // endpoints
    for ( int v = 0; v < n; v++ ) {
      fprintf( output_file, "n %d %.0f %.0lf\n", // ell f - not sure what
                                                 // difference it makes
               v + 1, p_array[v].x, p_array[v].y ); 
    }
    for (int i = 0; i < en; i++)
    {
        fprintf(output_file, "e %d %d %d\n",
                E[i].u + 1, E[i].v + 1, (int) (E[i].w) + 1);
    }
    fclose(output_file);
}


static edge *e_array;
static edge **free_list_e;
static int n_free_e;

void alloc_memory(int n)
{
    E = (EDGE *)calloc(3 * n + 1, sizeof(EDGE));
    dist = (double *)calloc( n + 1, sizeof(double));
    has_predecessor = (bool *)calloc( n + 1, sizeof(bool));
    for ( int i = 0; i < n; i++ ) has_predecessor[i] = false;
    //mst = (EDGE *)calloc(3 * n + 1, sizeof(EDGE));
    //p = (int *)calloc( n + 1, sizeof(int));
    //d = (int *)calloc( n + 1, sizeof(int));
    edge *e;
    int i;
    p_array = (point *)calloc(n, sizeof(point));
    n_free_e = 3 * n;
    e_array = e = (edge *)calloc(n_free_e, sizeof(edge));
    free_list_e = (edge **)calloc(n_free_e, sizeof(edge *));
    for(i = 0; i < n_free_e; i++, e++)
        free_list_e[i] = e;
}

void free_memory(void)
{
    free(p_array);
    free(e_array);
    free(free_list_e);
    free(has_predecessor);
}

edge *get_edge(void)
{
    if(n_free_e == 0)
        printf("Out of memory for edges\n");
    return (free_list_e[--n_free_e]);
}

void free_edge(edge *e)
{
    free_list_e[n_free_e++] = e;
}

void splice(edge *a, edge *b, point *v)
{
    edge *next;
    if(Org(a) == v)
    {
        next = Onext(a);
        Onext(a) = b;
    }
    else
    {
        next = Dnext(a);
        Dnext(a) = b;
    }
    if(Org(next) == v)
        Oprev(next) = b;
    else
        Dprev(next) = b;
    if(Org(b) == v)
    {
        Onext(b) = next;
        Oprev(b) = a;
    }
    else
    {
        Dnext(b) = next;
        Dprev(b) = a;
    }
}

edge *make_edge(point *u, point *v)
{
    edge *e;
    e = get_edge();
    e->onext = e->oprev = e->dnext = e->dprev = e;
    /* ensure that the edge goes is directed so that the graph forms a dag */
    if ( pcmp( u, v ) <= 0 ) {
      e->org = u;
      e->dest = v;
    }
    else {
      e->org = v;
      e->dest = u;
    }
    if(u->entry_pt == NULL)
        u->entry_pt = e;
    if(v->entry_pt == NULL)
        v->entry_pt = e;
    return e;
}

edge *join(edge *a, point *u, edge *b, point *v, side s)
{
    edge *e;
    e = make_edge(u, v);
    if(s == left)
    {
        if (Org(a) == u)
            splice(Oprev(a), e, u);
        else
            splice(Dprev(a), e, u);
        splice(b, e, v);
    }
    else
    {
        splice(a, e, u);
        if(Org(b) == v)
            splice(Oprev(b), e, v);
        else
            splice(Dprev(b), e, v);
    }
    return e;
}

void delete_edge(edge *e)
{
    point *u, *v;
    u = Org(e);
    v = Dest(e);
    if(u->entry_pt == e)
        u->entry_pt = e->onext;
    if(v->entry_pt == e)
        v->entry_pt = e->dnext;
    if(Org(e->onext) == u)
        e->onext->oprev = e->oprev;
    else
        e->onext->dprev = e->oprev;
    if(Org(e->oprev) == u)
        e->oprev->onext = e->onext;
    else
        e->oprev->dnext = e->onext;
    if(Org(e->dnext) == v)
        e->dnext->oprev = e->dprev;
    else
        e->dnext->dprev = e->dprev;
    if(Org(e->dprev) == v)
        e->dprev->onext = e->dnext;
    else
        e->dprev->dnext = e->dnext;
    free_edge(e);
}

/**
 *  @todo very confusing use of max (use construction in geometric.c!)
 */
double max(double a, double b)
{
    a = a > 0 ? a : -a;
    b = b > 0 ? b : -b;
    return a > b ? a : b;
}

static void create_edges(int n)
{
    edge *e_start, *e;
    point *u, *v;
    int i;
    for(i = 0; i < n; i++)
    {
        u = &p_array[i];
        e_start = e = u->entry_pt;
        do
        {
            v = Other_point(e, u);
            if(u < v)
            {
              /**
               * @todo This is pretty obscure. The number of a vertex is its
               * position in p_array. The variables u and v are points, but
               * E[en].u and E[en].v are vertex numbers! This code probably
               * also ensures that edges go in the right direction -- since
               * p_array is sorted -- but it's not obvious. Also, here the
               * vertex numbers are 0-based; there's a '+ 1' when they are
               * printed.
               */
                E[en].u = u - p_array, E[en].v = v - p_array;
                E[en].w = Distance(u, v);
                has_predecessor[E[en].v] = true;
                en++;
            }
            e = Next(e, u);
        }while(!Identical_refs(e, e_start));
    }
    // create extra edges to make sure that the outer face is a triangle by
    // connecting vertices that don't have predecessors to the vertex
    // corresponding to the first point.
/*     int first_point_x = p_array[0].x; */
/*     int first_point_y = p_array[0].y; */
/*     for ( int v = 1; v < n; v++ ) { */
/*       if ( ! has_predecessor[v] ) { */
/*         E[en].u = 0; E[en].v = v; */
/*         int v_x = p_array[v].x; */
/*         int v_y = p_array[v].y; */
/*         E[en++].w = max(first_point_x - v_x, first_point_y - v_y); */
/*       } */
/*     } */
}

void merge_sort(point *p[], point *p_temp[], int l, int r)
{
    int i, j, k, m;
    if(r - l > 0)
    {
        m = (r + l) / 2;
        merge_sort(p, p_temp, l, m);
        merge_sort(p, p_temp, m + 1, r);
        for(i = m + 1; i > l; i--)
            p_temp[i - 1] = p[i - 1];
        for(j = m; j < r; j++)
            p_temp[r + m - j] = p[j + 1];
        for(k = l; k <= r; k++)
            if(p_temp[i]->x < p_temp[j]->x)
            {
                p[k] = p_temp[i];
                i = i + 1;
            }
            else if(p_temp[i]->x == p_temp[j]->x && p_temp[i]->y < p_temp[j]->y)
            {
                p[k] = p_temp[i];
                i = i + 1;
            }
            else
            {
                p[k] = p_temp[j];
                j = j - 1;
            }
    }
}

static void lower_tangent(edge *r_cw_l, point *s, edge *l_ccw_r, point *u, edge **l_lower, point **org_l_lower, edge **r_lower, point **org_r_lower)
{
    edge *l, *r;
    point *o_l, *o_r, *d_l, *d_r;
    bool finished;
    l = r_cw_l;
    r = l_ccw_r;
    o_l = s;
    d_l = Other_point(l, s);
    o_r = u;
    d_r = Other_point(r, u);
    finished = false;
    while(!finished)
        if(Cross_product_3p(o_l, d_l, o_r) > 0.0)
        {
            l = Prev(l, d_l);
            o_l = d_l;
            d_l = Other_point(l, o_l);
        }
        else if(Cross_product_3p(o_r, d_r, o_l) < 0.0)
        {
            r = Next(r, d_r);
            o_r = d_r;
            d_r = Other_point(r, o_r);
        }
        else
            finished = true;
        *l_lower = l;
        *r_lower = r;
        *org_l_lower = o_l;
        *org_r_lower = o_r;
}

static void merge(edge *r_cw_l, point *s, edge *l_ccw_r, point *u, edge **l_tangent)
{
    edge *base, *l_cand, *r_cand;
    point *org_base, *dest_base;
    Real u_l_c_o_b, v_l_c_o_b, u_l_c_d_b, v_l_c_d_b;
    Real u_r_c_o_b, v_r_c_o_b, u_r_c_d_b, v_r_c_d_b;
    Real c_p_l_cand, c_p_r_cand;
    Real d_p_l_cand, d_p_r_cand;
    bool above_l_cand, above_r_cand, above_next, above_prev;
    point *dest_l_cand, *dest_r_cand;
    Real cot_l_cand, cot_r_cand;
    edge *l_lower, *r_lower;
    point *org_r_lower, *org_l_lower;
    lower_tangent(r_cw_l, s, l_ccw_r, u, &l_lower, &org_l_lower, &r_lower, &org_r_lower);
    base = join(l_lower, org_l_lower, r_lower, org_r_lower, right);
    org_base = org_l_lower;
    dest_base = org_r_lower;
    *l_tangent = base;
    do
    {
        l_cand = Next(base, org_base);
        r_cand = Prev(base, dest_base);
        dest_l_cand = Other_point(l_cand, org_base);
        dest_r_cand = Other_point(r_cand, dest_base);
        Vector(dest_l_cand, org_base, u_l_c_o_b, v_l_c_o_b);
        Vector(dest_l_cand, dest_base, u_l_c_d_b, v_l_c_d_b);
        Vector(dest_r_cand, org_base, u_r_c_o_b, v_r_c_o_b);
        Vector(dest_r_cand, dest_base, u_r_c_d_b, v_r_c_d_b);
        c_p_l_cand = Cross_product_2v(u_l_c_o_b, v_l_c_o_b, u_l_c_d_b, v_l_c_d_b);
        c_p_r_cand = Cross_product_2v(u_r_c_o_b, v_r_c_o_b, u_r_c_d_b, v_r_c_d_b);
        above_l_cand = c_p_l_cand > 0.0;
        above_r_cand = c_p_r_cand > 0.0;
        if(!above_l_cand && !above_r_cand)
            break;
        if(above_l_cand)
        {
            Real u_n_o_b, v_n_o_b, u_n_d_b, v_n_d_b;
            Real c_p_next, d_p_next, cot_next;
            edge *next;
            point *dest_next;
            d_p_l_cand = Dot_product_2v(u_l_c_o_b, v_l_c_o_b, u_l_c_d_b, v_l_c_d_b);
            cot_l_cand = d_p_l_cand / c_p_l_cand;
            do
            {
                next = Next(l_cand, org_base);
                dest_next = Other_point(next, org_base);
                Vector(dest_next, org_base, u_n_o_b, v_n_o_b);
                Vector(dest_next, dest_base, u_n_d_b, v_n_d_b);
                c_p_next = Cross_product_2v(u_n_o_b, v_n_o_b, u_n_d_b, v_n_d_b);
                above_next = c_p_next > 0.0;
                if(!above_next)
                    break;
                d_p_next = Dot_product_2v(u_n_o_b, v_n_o_b, u_n_d_b, v_n_d_b);
                cot_next = d_p_next / c_p_next;
                if(cot_next > cot_l_cand)
                    break;
                delete_edge(l_cand);
                l_cand = next;
                cot_l_cand = cot_next;
            }while(true);
        }
        if(above_r_cand)
        {
            Real u_p_o_b, v_p_o_b, u_p_d_b, v_p_d_b;
            Real c_p_prev, d_p_prev, cot_prev;
            edge *prev;
            point *dest_prev;
            d_p_r_cand = Dot_product_2v(u_r_c_o_b, v_r_c_o_b, u_r_c_d_b, v_r_c_d_b);
            cot_r_cand = d_p_r_cand / c_p_r_cand;
            do
            {
                prev = Prev(r_cand, dest_base);
                dest_prev = Other_point(prev, dest_base);
                Vector(dest_prev, org_base, u_p_o_b, v_p_o_b);
                Vector(dest_prev, dest_base, u_p_d_b, v_p_d_b);
                c_p_prev = Cross_product_2v(u_p_o_b, v_p_o_b, u_p_d_b, v_p_d_b);
                above_prev = c_p_prev > 0.0;
                if(!above_prev)
                    break;
                d_p_prev = Dot_product_2v(u_p_o_b, v_p_o_b, u_p_d_b, v_p_d_b);
                cot_prev = d_p_prev / c_p_prev;
                if(cot_prev > cot_r_cand)
                    break;
                delete_edge(r_cand);
                r_cand = prev;
                cot_r_cand = cot_prev;
            }while(true);
        }
        dest_l_cand = Other_point(l_cand, org_base);
        dest_r_cand = Other_point(r_cand, dest_base);
        if(!above_l_cand || (above_l_cand && above_r_cand && cot_r_cand < cot_l_cand))
        {
            base = join(base, org_base, r_cand, dest_r_cand, right);
            dest_base = dest_r_cand;
        }
        else
        {
            base = join(l_cand, dest_l_cand, base, dest_base, right);
            org_base = dest_l_cand;
        }
    }while(true);
}

void divide(point *p_sorted[], int l, int r, edge **l_ccw, edge **r_cw)
{
    int n;
    int split;
    edge *l_ccw_l, *r_cw_l, *l_ccw_r, *r_cw_r, *l_tangent;
    edge *a, *b, *c;
    Real c_p;
    n = r - l + 1;
    if(n == 2)
    {
        *l_ccw = *r_cw = make_edge(p_sorted[l], p_sorted[r]);
    }
    else if(n == 3)
    {
        a = make_edge(p_sorted[l], p_sorted[l + 1]);
        b = make_edge(p_sorted[l + 1], p_sorted[r]);
        splice(a, b, p_sorted[l + 1]);
        c_p = Cross_product_3p(p_sorted[l], p_sorted[l + 1], p_sorted[r]);
        if(c_p > 0.0)
        {
            c = join(a, p_sorted[l], b, p_sorted[r], right);
            *l_ccw = a;
            *r_cw = b;
        }
        else if(c_p < 0.0)
        {
            c = join(a, p_sorted[l], b, p_sorted[r], left);
            *l_ccw = c;
            *r_cw = c;
        }
        else
        {
            *l_ccw = a;
            *r_cw = b;
        }
    }
    else if(n > 3)
    {
        split = (l + r) / 2;
        divide(p_sorted, l, split, &l_ccw_l, &r_cw_l);
        divide(p_sorted, split+1, r, &l_ccw_r, &r_cw_r);
        merge(r_cw_l, p_sorted[split], l_ccw_r, p_sorted[split + 1], &l_tangent);
        if(Org(l_tangent) == p_sorted[l])
            l_ccw_l = l_tangent;
        if(Dest(l_tangent) == p_sorted[r])
            r_cw_r = l_tangent;
        *l_ccw = l_ccw_l;
        *r_cw = r_cw_r;
    }
}

int main ( int argc, char * argv[] )
{
  /**
   * @todo check that argc is 1 or 2 and print a usage message otherwise
   */
  const char * output_file_name;
  if ( argc < 3 || argc > 4 ) {
    printf( "Usage: %s number_of_vertices seed [output_file]\n", argv[0] );
    printf( "       the graph is printed to standard output if output_file is missing\n");
    exit( 1 );
  }
  n = strtod( argv[1], NULL );
  seed = strtoul( argv[2], NULL, 10 );
  if ( argc == 4 ) {
    output_file_name = argv[3];
  }
  else {
    output_file_name = NULL;
  }
  init_genrand( seed );
    edge *l_cw, *r_ccw;
    int i;
    point **p_sorted, **p_temp;
    alloc_memory(n);
    for (i = 0; i < n; i++)
    {
        p_array[i].x = genrand_int32() % RANGE;
        p_array[i].y = genrand_int32() % RANGE;
    }
    for(i = 0; i < n; i++)
        p_array[i].entry_pt = NULL;
    p_sorted = (point **)malloc((unsigned)n * sizeof(point *));
    p_temp = (point **)malloc((unsigned)n * sizeof(point *));
    for(i = 0; i < n; i++)
        p_sorted[i] = p_array + i;
    merge_sort(p_sorted, p_temp, 0, n - 1);
    free((char *)p_temp);
    divide(p_sorted, 0, n - 1, &l_cw, &r_ccw);
    free((char *)p_sorted);
    en = 0;
    create_edges(n);
    output_graph(output_file_name);
    free_memory();
    return EXIT_SUCCESS;
}

/*  [Last modified: 2021 05 26 at 15:20:41 GMT] */
