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

#include "fastdate.h"

#define LINEALLOC 4096
#define REGEX_REAL "[-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?"

static FILE * fp = NULL;
static char line[LINEALLOC];
static long lineno = 0;

static regex_t regexp_node_exp;
static regex_t regexp_mrca_exp;
static regex_t regexp_node_ln;
static regex_t regexp_mrca_ln;
static regex_t regexp_comment;
static regex_t regexp_ignore;

static regmatch_t pmatch[9];

/* regular expressions for matching file entries */
static char * re_node_exp = "^\\s*(\\w+)\\s+exp\\s*\\(\\s*" 
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*\\)\\s*$";

static char * re_mrca_exp = "^\\s*(\\w+)\\s+(\\w+)\\s+exp\\s*\\(\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*\\)\\s*$";


static char * re_node_ln = "^\\s*(\\w+)\\s+ln\\s*\\(\\s*" 
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*\\)\\s*$";

static char * re_mrca_ln = "^\\s*(\\w+)\\s+(\\w+)\\s+ln\\s*\\(\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*,\\s*"
                            "(" REGEX_REAL ")"
                            "\\s*\\)\\s*$";

static char * re_comment = "^#.*$";
static char * re_ignore = "^\\s*$";


static void regexp_init()
{
  if (regcomp(&regexp_node_exp, re_node_exp, REG_EXTENDED))
    fatal("Bad pattern");
  if (regcomp(&regexp_mrca_exp, re_mrca_exp, REG_EXTENDED))
    fatal("Bad pattern");
  if (regcomp(&regexp_node_ln, re_node_ln, REG_EXTENDED))
    fatal("Bad pattern");
  if (regcomp(&regexp_mrca_ln, re_mrca_ln, REG_EXTENDED))
    fatal("Bad pattern");
  if (regcomp(&regexp_comment, re_comment, REG_EXTENDED))
    fatal("Bad pattern");
  if (regcomp(&regexp_ignore, re_ignore, REG_EXTENDED))
    fatal("Bad pattern");
}

static void regexp_free()
{
  regfree(&regexp_node_exp); 
  regfree(&regexp_mrca_exp); 
  regfree(&regexp_node_ln); 
  regfree(&regexp_mrca_ln); 
  regfree(&regexp_comment); 
  regfree(&regexp_ignore); 
}

static void priorfile_open()
{
  fp = fopen(opt_priorfile, "r");
  if (!fp)
    fatal("Cannot open file %s", opt_priorfile);
}

static tree_node_t * query_node(char * label)
{
  ENTRY query;
  ENTRY * found = NULL;

  query.key = label;
  found = hsearch(query,FIND);
  if (!found)
    return NULL;

  return (tree_node_t *)(found->data);
}

static void set_node_exp_prior(long * fossils_count,long * extinct_leaves_count)
{
  int len;
  tree_node_t * node;

  /* node label */
  len = pmatch[1].rm_eo - pmatch[1].rm_so;
  char * label = xstrndup(line+pmatch[1].rm_so, len);

  /* mean */
  len = pmatch[2].rm_eo - pmatch[2].rm_so;
  char * mean = xstrndup(line+pmatch[2].rm_so, len);

  /* offset */
  len = pmatch[4].rm_eo - pmatch[4].rm_so;
  char * offset = xstrndup(line+pmatch[4].rm_so, len);

  node = query_node(label);
  if (!node)
    fatal("Node %s does not exist (line %d of %s)", 
          label, lineno, opt_priorfile);

  if (node->prior_lineno)
    fatal("Error: line %d of file %s assigns a prior to node %s which already "
          "has a prior from line %d", 
          lineno, opt_priorfile, label, node->prior_lineno);


  exp_params_t * params = (exp_params_t *)xmalloc(sizeof(exp_params_t));

  params->mean   = atof(mean);
  params->offset = atof(offset);  

  node->prior        = NODEPRIOR_EXP;
  node->prior_params = params;
  node->prior_lineno = lineno;

  if (!opt_quiet)
    printf ("  exp(%s) on node (%s) with offset %s\n", mean, label, offset); 

  if (!node->left)
    *extinct_leaves_count = *extinct_leaves_count + 1;
  else
    *fossils_count = *fossils_count + 1;

  free(label);
  free(mean);
  free(offset);
}

/* set log-normal prior for a node */
static void set_node_ln_prior(long * fossils_count, long * extinct_leaves_count)
{
  int len;
  tree_node_t * node;

  /* node label */
  len = pmatch[1].rm_eo - pmatch[1].rm_so;
  char * label = xstrndup(line+pmatch[1].rm_so, len);

  /* mean */
  len = pmatch[2].rm_eo - pmatch[2].rm_so;
  char * mean = xstrndup(line+pmatch[2].rm_so, len);

  /* standard deviation */
  len = pmatch[4].rm_eo - pmatch[4].rm_so;
  char * stdev = xstrndup(line+pmatch[4].rm_so, len);

  /* offset */
  len = pmatch[6].rm_eo - pmatch[6].rm_so;
  char * offset = xstrndup(line+pmatch[6].rm_so, len);

  node = query_node(label);
  if (!node)
    fatal("Node %s does not exist (line %d of %s)", 
          label, lineno, opt_priorfile);

  if (node->prior_lineno)
    fatal("Error: line %d of file %s assigns a prior to node %s which already "
          "has a prior from line %d",
          lineno, opt_priorfile, label, node->prior_lineno);

  ln_params_t * params = (ln_params_t *)xmalloc(sizeof(ln_params_t));
  params->mean   = atof(mean);
  params->stdev  = atof(stdev);
  params->offset = atof(offset);

  node->prior        = NODEPRIOR_LN;
  node->prior_params = params;
  node->prior_lineno = lineno;

  if (!opt_quiet)
    printf ("  ln(%s,%s) on node (%s) with offset %s\n",
            mean, stdev, label, offset);

  if (!node->left)
    *extinct_leaves_count = *extinct_leaves_count + 1;
  else
    *fossils_count = *fossils_count + 1;

  free(label);
  free(mean);
  free(stdev);
  free(offset);
}

static void set_mrca_exp_prior(long * fossils_count)
{
  int len;
  tree_node_t * node1;
  tree_node_t * node2;
  tree_node_t * lca;

  /* node label */
  len = pmatch[1].rm_eo - pmatch[1].rm_so;
  char * tip1 = xstrndup(line+pmatch[1].rm_so, len);

  /* node label */
  len = pmatch[2].rm_eo - pmatch[2].rm_so;
  char * tip2 = xstrndup(line+pmatch[2].rm_so, len);

  /* mean */
  len = pmatch[3].rm_eo - pmatch[3].rm_so;
  char * mean = xstrndup(line+pmatch[3].rm_so, len);

  /* mean */
  len = pmatch[5].rm_eo - pmatch[5].rm_so;
  char * offset = xstrndup(line+pmatch[5].rm_so, len);

  /* match tip1 to its node */
  node1 = query_node(tip1);
  if (!node1)
    fatal("Tip %s does not exist (line %d of %s)", tip1, lineno, opt_priorfile);

  /* match tip1 to its node */
  node2 = query_node(tip2);
  if (!node2)
    fatal("Tip %s does not exist (line %d of %s)", tip2, lineno, opt_priorfile);

  /* find their lowest common ancestor */
  lca = lca_compute(node1, node2);

  if (lca->prior_lineno)
    fatal("Error: line %d of file %s assigns a prior to the mRCA of %s and %s, "
          "which already has a prior from line %d", 
          lineno, opt_priorfile, tip1, tip2, lca->prior_lineno);

  /* set distribution parameters at concrete node */
  exp_params_t * params = (exp_params_t *)xmalloc(sizeof(exp_params_t));
  params->mean   = atof(mean);
  params->offset = atof(offset);

  lca->prior        = NODEPRIOR_EXP;
  lca->prior_params = params;
  lca->prior_lineno = lineno;

  if (!opt_quiet)
    printf ("  exp(%s) on mrca (%s,%s) with offset %s\n",
            mean, tip1, tip2, offset); 

  *fossils_count = *fossils_count + 1;

  free(tip1);
  free(tip2);
  free(mean);
  free(offset);
}

static void set_mrca_ln_prior(long * fossils_count)
{
  int len;
  tree_node_t * node1;
  tree_node_t * node2;
  tree_node_t * lca;

  /* node label */
  len = pmatch[1].rm_eo - pmatch[1].rm_so;
  char * tip1 = xstrndup(line+pmatch[1].rm_so, len);

  /* node label */
  len = pmatch[2].rm_eo - pmatch[2].rm_so;
  char * tip2 = xstrndup(line+pmatch[2].rm_so, len);

  /* mean */
  len = pmatch[3].rm_eo - pmatch[3].rm_so;
  char * mean = xstrndup(line+pmatch[3].rm_so, len);

  /* standard deviation */
  len = pmatch[5].rm_eo - pmatch[5].rm_so;
  char * stdev = xstrndup(line+pmatch[5].rm_so, len);

  /* offset */
  len = pmatch[7].rm_eo - pmatch[7].rm_so;
  char * offset = xstrndup(line+pmatch[7].rm_so, len);

  /* match tip1 to its node */
  node1 = query_node(tip1);
  if (!node1)
    fatal("Tip %s does not exist (line %d of %s)", tip1, lineno, opt_priorfile);

  /* match tip1 to its node */
  node2 = query_node(tip2);
  if (!node2)
    fatal("Tip %s does not exist (line %d of %s)", tip2, lineno, opt_priorfile);

  /* find their lowest common ancestor */
  lca = lca_compute(node1, node2);

  if (lca->prior_lineno)
    fatal("Error: line %d of file %s assigns a prior to the mRCA of %s and %s, "
          "which already has a prior from line %d", 
          lineno, opt_priorfile, tip1, tip2, lca->prior_lineno);

  /* set distribution parameters at concrete node */
  ln_params_t * params = (ln_params_t *)xmalloc(sizeof(ln_params_t));
  params->mean   = atof(mean);
  params->stdev  = atof(stdev);
  params->offset = atof(offset);

  lca->prior        = NODEPRIOR_LN;
  lca->prior_params = params;
  lca->prior_lineno = lineno;

  if (!opt_quiet)
    printf ("  ln(%s,%s) on mrca(%s,%s) with offset %s\n",
            mean, stdev, tip1, tip2, offset);

  *fossils_count = *fossils_count + 1;

  free(tip1);
  free(tip2);
  free(mean);
  free(stdev);
  free(offset);
}

static void priorfile_parse(long * fossils_count, long * extinct_leaves_count)
{
  while(fgets(line, LINEALLOC, fp))
  {
    ++lineno;

    if (regexec(&regexp_node_exp, line, 6, pmatch, 0) == 0)
      set_node_exp_prior(fossils_count,extinct_leaves_count);
    else if (regexec(&regexp_mrca_exp, line, 7, pmatch, 0) == 0)
      set_mrca_exp_prior(fossils_count);
    else if (regexec(&regexp_node_ln, line, 7, pmatch, 0) == 0)
      set_node_ln_prior(fossils_count,extinct_leaves_count);
    else if (regexec(&regexp_mrca_ln, line, 8, pmatch, 0) == 0)
      set_mrca_ln_prior(fossils_count);
    else if (regexec(&regexp_comment, line, 0, NULL, 0) == 0)
      continue;
    else if (regexec(&regexp_ignore, line, 0, NULL, 0) == 0)
      continue;
    else
      fatal("Invalid syntax on line %d of file %s", lineno, opt_priorfile);
  }
}

static void priorfile_close(void)
{
  fclose(fp); 
}

static void hashtable_init(unsigned int node_count, tree_node_t ** nodes)
{
  ENTRY entry;
  unsigned int i;

  hcreate(node_count * 2);

  for (i = 0; i < node_count; ++i)
  {
    if (!nodes[i]->label)
      continue;
      
    entry.key  = nodes[i]->label;
    entry.data = (void *)nodes[i];
    hsearch(entry,ENTER);
  }
}

void set_node_priors(tree_node_t * root, 
                     long * fossils_count, 
                     long * extinct_leaves_count)
{
  unsigned int node_count = 2*root->leaves - 1;

  *fossils_count = 0;
  *extinct_leaves_count = 0;

  /* create a linear list of tree node pointers */
  tree_node_t ** nodes = (tree_node_t **)xmalloc(node_count * 
                                                 sizeof(tree_node_t *));
  tree_traverse(root, nodes);

  /* create hash table with node labels */
  hashtable_init(node_count, nodes);
    
  /* initialize LCA computation */
  lca_init(root);

  /* compile regular expressions */
  regexp_init();
  
  /* open prior file */
  priorfile_open();

  /*  parse prior file */
  priorfile_parse(fossils_count, extinct_leaves_count);

  /* close prior file */
  priorfile_close();

  /* destroy initialized LCA data */
  lca_destroy();

  /* destroy linear list of tree node pointers */
  free(nodes);

  /* deallocate compiled regular expression structures */
  regexp_free();

  /* destroy hash table */
  hdestroy();
}
