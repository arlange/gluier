//
// File: g.h
//

#ifndef _G_H
#define _G_H

#include <vector>
#include <list>
#include <iostream>
#include <string>
#include "set_operations2.h"

using namespace std;

class g {

 public:
  g( int vsize );
  g( const g &otherG );
  ~g();
  struct edge;

 public:
  int order() const;
  int num_edges();
  int num_tris();
  vset get_adj( int v );
  void add_edge( int u, int v );
  void remove_edge( int u, int v);
  void remove_circ_edge( int e );
  bool is_edge( int u, int v ) const;
  int min_degree();
  int degree( int v );
  vector<int> neighbors( int v );
  vector<int> max_clique( bool print = false, int k = -1 );
  vector<int> max_independent_set( bool print = false, int k = -1);
  bool has_clique( int k, bool is = false );
  void remove_vs( vector<int> vs, int k );
  void make_complement();
  void make_complete();
  int make_rand_er( float sigma );
  bool glue_graphs( g * x, g * y, vector<int> cones, int p );
  bool v_extend( int dx, g * y, vset cone );
  bool glue_graphs( int d, int k, g * y, vector<int>cones );
  bool has_c( int c = 4 );
  void print( ostream * o = &cout );
  string to_g6();
  void read_g6( string g6 );
  void print_g6( ostream *o = &cout );
  string to_canon_g6();
  int get_p3s( int * tab, int p );
  int get_p3s2( bool * tab, int p );
  void get_closures( uint64_t * tab, int p );
  void get_closures2( uint32_t * tab, int p );
  void get_independences( int * tab, int p, int max_is = -1);
  void get_independences2( int * tab, int p, int max_is = -1);
  void get_independences3( int * tab, int p, int max_is = -1);
  void get_independences4( uint8_t * tab, int p, int max_is = -1);
  void get_is( int max_is );
  int get_num_tri_edges( int e );

  // irredundant set stuff
  int max_irs( int max_ir_possible );
  int max_irs( int max_ir_possible, vset of_these );	
  int max_irs( int max_ir_possible, vector<vset>& max_irs );
  int get_irs( int * tab, vector<vset>& max_irs, int max_ir_possible );
  int get_independent_sets( int * tab, int max_is );
  vector<vset> vec_is_sets;
  vset private_neighbors( int v, vset S );
  bool not_isolate( int v, vset S );
  int get_p5s( bool * tab, int p );
  void recursive_p5( vset cur_cl, int cur_v );

 private:
  void set_up();
  void max_clique_backtrack( int l, int k );
  void all_cliques_backtrack( uint64_t cur_cl, int l, int k, int v );
  void all_cliques( int k, int p );
  void all_independent_sets( int k, int p );
  void get_tri_stats();
  void recursive_is( uint64_t cur_cl, int k, int cur_v );
  void recursive_p3( uint64_t cur_cl, int cur_v );
  int recursive_clos( uint32_t s );

  // irredundant set stuff
  void all_ir_bt(vset cur_ir, int l, int v);
  void all_irs(int max_ir);
  void recursive_ir( uint64_t cur_cl, int k, int cur_v );


 private:
  int n, arraySize, numEdges, oldN;
  int * vdegree;
  vector<vset> gA;
  vector<vset> gB;
  vset gV;
  vector<uint64_t> * cliques;
  struct set_list;
  int * tri_stats;
  bool trisCalced;

  // used for max-clique algorithm
  vector<vset> gC;
  vector< int > X, optClique;
  int optSize;
  uint8_t * cur_tab;
  bool * bcur_tab;
  uint32_t * cur_utab;
  int * clique_num;



};

#endif
