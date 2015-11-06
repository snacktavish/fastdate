/*
    Copyright (C) 2015 Tomas Flouri

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: Tomas Flouri <Tomas.Flouri@h-its.org>,
    Heidelberg Institute for Theoretical Studies,
    Schloss-Wolfsbrunnenweg 35, D-69118 Heidelberg, Germany
*/

#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <pthread.h>
#include <getopt.h>
#include <regex.h>
#include <x86intrin.h>
#include <search.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <pthread.h>

/* constants */

#define PROG_NAME "fastdate"
#define PROG_VERSION "v0.0.0"

#ifdef __APPLE__
#define PROG_ARCH "macosx_x86_64"
#else
#define PROG_ARCH "linux_x86_64"
#endif

/* structures and data types */

typedef unsigned int UINT32;
typedef unsigned short WORD;
typedef unsigned char BYTE;

/* parameters of exponential distribution node prior */
typedef struct exp_params_s
{
  double mean;
  double offset;
} exp_params_t;

/* parameters of uniform distribution node prior */
typedef struct uni_params_s
{
  double min_age;
  double max_age;
} uni_params_t;

/* parameters of lognormal distribution node prior */
typedef struct ln_params_s
{
  double mean;
  double stdev;
  double offset;
} ln_params_t;

/* parameters of normal distribution node prior */
typedef struct normal_params_s
{
  double mean;
  double variance;
  double offset;
} normal_params_t;

typedef struct prior_s
{
  int dist;
  int lineno;
  char * taxa;
  void * params;
} prior_t;

typedef struct list_s
{
  prior_t * prior;
  struct list_s * next;
} list_t;

typedef struct tree_noderec
{
  char * label;
  double length;
  struct tree_noderec * left;
  struct tree_noderec * right;
  struct tree_noderec * parent;
  long height;
  long leaves;

  /* grid related data */
  long entries;
  double * matrix;
  long * matrix_left;
  long * matrix_right;

  /* age specific */
  long interval_line;

  /* prior specific */
  int prior;
  int prior_lineno;
  void * prior_params;

  /* sampling specific */
  long sampled_gridline;

} tree_node_t;

/* definitions */

#define OUTPUT_ULTRAMETRIC      0
#define OUTPUT_DATED            1

#define NODEPRIOR_NONE    0
#define NODEPRIOR_EXP     1
#define NODEPRIOR_LN      2
#define NODEPRIOR_UNI     3
#define NODEPRIOR_NORMAL  4

#define PARAM_LAMBDA   1
#define PARAM_MU       2
#define PARAM_PSI      4
#define PARAM_RHO      8

#define DEFAULT_LAMBDA 10.0
#define DEFAULT_MU     0.001
#define DEFAULT_PSI    0.5
#define DEFAULT_RHO    0.5

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* options */

extern int opt_quiet;
extern int opt_outform;
extern int opt_help;
extern int opt_version;
extern int opt_method_relative;
extern int opt_method_nodeprior;
extern int opt_method_tipdates;
extern int opt_showtree;
extern char * opt_treefile;
extern char * opt_outfile;
extern char * opt_priorfile;
extern long opt_seed;
extern long opt_sample;
extern long opt_grid_intervals;
extern long opt_threads;
extern double opt_max_age;
extern double opt_lambda;
extern double opt_mu;
extern double opt_rho;
extern double opt_psi;
extern double opt_rate_mean;
extern double opt_rate_var;
extern double opt_cred_interval;

extern unsigned int opt_parameters_bitv;

/* matrices */

extern const char map_nt[256];
extern const char map_bin[256];
extern const char map_nt_2bit[256];

/* common data */

extern long mmx_present;
extern long sse_present;
extern long sse2_present;
extern long sse3_present;
extern long ssse3_present;
extern long sse41_present;
extern long sse42_present;
extern long popcnt_present;
extern long avx_present;
extern long avx2_present;

/* functions in util.c */

void fatal(const char * format, ...) __attribute__ ((noreturn));
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned long progress);
void progress_done(void);
void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
void encode_sequence(char * s, const char * map);
long getusec(void);
void show_rusage(void);

/* functions in fastdate.c */

void args_init(int argc, char ** argv);
void cmd_help(void);
void getentirecommandline(int argc, char * argv[]);
void fillheader(void);
void show_header(void);
void cmd_method_relative(void);
void set_seed(void);
void cmd_method_nodeprior(void);
void cmd_method_tipdates(void);

/* functions in tree.c */

tree_node_t * yy_create_tree(void);
void yy_dealloc_tree(tree_node_t * tree);
void show_ascii_tree(tree_node_t * tree);
long set_node_heights(tree_node_t * root);
void write_newick_tree(tree_node_t * node);
int tree_traverse(tree_node_t * root, tree_node_t ** outbuffer);
int rtree_query_tipnodes(tree_node_t * root, tree_node_t ** node_list);

/* functions in dp.c */

void dp(tree_node_t * tree);
double dp_evaluate(tree_node_t * tree);

/* functions in gamma.c */

void gamma_dist_init(void);
double gamma_dist_logpdf(double x);

/* functions in bd.c */

void bd_init(long fossils_count, long extinct_leaves_count);
double bd_relative_prod(double t);
double bd_relative_root(long leaves, double t);
double bd_tipdates_root(long leaves, double t);
double bd_tipdates_prod_inner(double t);
double bd_tipdates_prod_tip(double t);


/* functions in newick.y */

tree_node_t * yy_parse_tree(const char * filename);

/* functions in arch.c */

unsigned long arch_get_memused(void);
unsigned long arch_get_memtotal(void);

/* functions in exp.c */

double exp_dist_pdf(double lambda, double x);
double exp_dist_logpdf(double lambda, double x);

/* functions in ln.c */

double ln_dist_pdf(double mean, double variance, double x);
double ln_dist_logpdf(double mean, double variance, double x);

/* functions in normal.c */

double normal_dist_pdf(double mean, double variance, double x);
double normal_dist_logpdf(double mean, double variance, double x);

/* functions in uni.c */

double uni_dist_pdf(double a, double b, double x);
double uni_dist_logpdf(double a, double b, double x);

/* functions in nodeprior.c */

void set_node_priors(tree_node_t * root,
                     long * fossils_count,
                     long * extinct_leaves_count);

/* functions in lca.c */

tree_node_t * lca_compute(tree_node_t * root,
                          tree_node_t ** tip_nodes,
                          unsigned int count);

/* functions in optimize.c */

double opt_parameters(tree_node_t * tree, int which, double factr, double pgtol);

/* functions in sample.c */

void sample(tree_node_t * root);

/* functions in parse_prior.y */

list_t * yy_parse_nodeprior(const char * filename);
