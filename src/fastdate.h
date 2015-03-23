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
#include <x86intrin.h>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include <locale.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>

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

typedef struct tree_noderec
{
  char * label;
  double length;
  struct tree_noderec * left;
  struct tree_noderec * right;
  int height;
  int leaves;

  /* grid related data */
  int entries;
  double * matrix;
  int * matrix_left;
  int * matrix_right;

  /* age specific */
  int interval_line;

  void * data;
} tree_node_t;

/* definitions */

#define OUTPUT_ULTRAMETRIC      0
#define OUTPUT_DATED            1

/* macros */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* options */

extern int opt_quiet;
extern int opt_exhaustive_bd;
extern int opt_outformat;
extern char * opt_treefile;
extern char * opt_outfile;
extern long opt_grid_intervals;
extern long opt_help;
extern long opt_version;
extern long opt_divtimes;
extern long opt_show_tree;
extern long opt_threads;
extern double opt_birth_rate;
extern double opt_death_rate;
extern double opt_edgerate_mean;
extern double opt_edgerate_var;

/* matrices */

extern const char map_nt[256];
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

void fatal(const char * format, ...);
void progress_init(const char * prompt, unsigned long size);
void progress_update(unsigned int progress);
void progress_done();
void * xmalloc(size_t size);
void * xrealloc(void *ptr, size_t size);
char * xstrchrnul(char *s, int c);
char * xstrdup(const char * s);
char * xstrndup(const char * s, size_t len);
void encode_sequence(char * s, const char * map);
long getusec(void);
void show_rusage();

/* functions in fastdate.c */

void args_init(int argc, char ** argv);
void cmd_help();
void getentirecommandline(int argc, char * argv[]);
void fillheader();
void show_header();
void cmd_divtimes();

/* functions in tree.c */

tree_node_t * yy_create_tree();
void yy_dealloc_tree(tree_node_t * tree);
void show_ascii_tree(tree_node_t * tree);
int set_node_heights(tree_node_t * root);
void write_newick_tree(tree_node_t * node);

/* functions in dp.c */

void dp(tree_node_t * tree);

/* functions in gamma.c */

void gamma_dist_init(double mean, double variance);
double gamma_dist_logpdf(double x);

/* functions in bd.c */

void bd_init(double birth_rate, double death_rate);
double bd_prob(int leaves, double t);

/* functions in newick.y */

tree_node_t * yy_parse_tree(const char * filename);

/* functions in arch.c */

unsigned long arch_get_memused();
unsigned long arch_get_memtotal();
