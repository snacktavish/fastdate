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

static void hashtable_create(tree_node_t * root)
{
  ENTRY entry;
  size_t i;

  tree_node_t ** node_list = (tree_node_t **)xmalloc((size_t)(root->leaves) *
                                                 sizeof(tree_node_t *));
  rtree_query_tipnodes(root, node_list);

  hcreate((size_t)(2*root->leaves));

  for (i=0; i < (size_t)(root->leaves); ++i)
  {
    if (!node_list[i]->label)
      continue;
      
    entry.key  = node_list[i]->label;
    entry.data = (void *)node_list[i];
    hsearch(entry,ENTER);
  }

  free(node_list);
}

static tree_node_t ** hashtable_query(char * tipstring,
                                     unsigned int * taxa_count)
{
  unsigned int k;
  unsigned int i;
  size_t taxon_len;
  unsigned int commas_count = 0;
  ENTRY * found = NULL;
  char * taxon;

  /* compute comma count */
  for (i = 0; i < strlen(tipstring); ++i)
    if (tipstring[i] == ',')
      commas_count++;

  /* allocate list of nodes that will be returned */
  tree_node_t ** out_node_list = (tree_node_t **)xmalloc((commas_count+1) *
                                                         sizeof(tree_node_t *));

  char * s = tipstring;
  k = 0;
  while (*s)
  {
    /* get next tip */
    taxon_len = strcspn(s, ",");
    if (!taxon_len)
      fatal("Erroneous prune list format (double comma)/taxon missing");

    taxon = strndup(s, taxon_len);

    /* search tip in hash table */
    ENTRY query;
    query.key = taxon;
    found = NULL;
    found = hsearch(query,FIND);
    
    if (!found)
      fatal("Cannot find taxon with node label (%s) in tree", taxon);

    /* store pointer in output list */
    out_node_list[k++] = (tree_node_t *)(found->data);

    /* free tip label, and move to the beginning of next tip if available */
    free(taxon);
    s += taxon_len;
    if (*s == ',') 
      s += 1;
  }

  /* return number of tips in the list */
  *taxa_count = commas_count + 1;

  /* return tip node list */
  return out_node_list;
}

static void print_prior(tree_node_t ** tip_list, unsigned int count, int priorno, prior_t * prior)
{
  unsigned int i;

  exp_params_t * exp_params;
  ln_params_t * ln_params;
  uni_params_t * uni_params;

  switch(prior->dist)
  {
    case NODEPRIOR_EXP:
      exp_params = (exp_params_t *)(prior->params);
      printf ("  %d: line %d -- exp(%f) offset %f MRCA (%s",
              priorno,
              prior->lineno,
              exp_params->mean,
              exp_params->offset,
              tip_list[0]->label);
      break;
    case NODEPRIOR_LN:
      ln_params = (ln_params_t *)(prior->params);
      printf ("%d: line %d -- ln(%f,%f) offset %f MRCA (%s",
              priorno,
              prior->lineno,
              ln_params->mean,
              ln_params->stdev,
              ln_params->offset,
              tip_list[0]->label);
      break;
    case NODEPRIOR_UNI:
      uni_params = (uni_params_t *)(prior->params);
      printf ("%d: line %d -- uni(%f,%f) MRCA (%s",
              priorno,
              prior->lineno,
              uni_params->min_age,
              uni_params->max_age,
              tip_list[0]->label);
      break;
    default:
      assert(0);
  }
  for (i = 1; i < count; ++i)
    printf(",%s", tip_list[i]->label);
  printf(")\n");
}

static void process_priors(list_t * prior_list,
                           tree_node_t * root,
                           long * fossils_count,
                           long * extinct_leaves_count)
{
  prior_t * prior;
  unsigned int taxa_count;
  tree_node_t * lca;
  int i = 0;

  while (prior_list)
  {
    i++;
    /* get current prior from list */
    prior = prior_list->prior;

    /* convert comma-separated taxa list to node pointers */
    tree_node_t ** tip_list = hashtable_query(prior->taxa, &taxa_count);

    /* store fossil / extinct leaves count */
    if (taxa_count == 1)
      *extinct_leaves_count = *extinct_leaves_count + 1;
    else
      *fossils_count = *fossils_count + 1;

    /* get the node pointer of the tip_list LCA */
    if (taxa_count == 1)
      lca = tip_list[0];
    else
      lca = lca_compute(root, tip_list, taxa_count);
    
    /* if a prior on that LCA already exists, then error */
    if (lca->prior_lineno)
      fatal("Error: line %d of file %s assigns a prior to the mRCA of %s and %s, "
            "which already has a prior from line %d", 
            prior->lineno, opt_priorfile, tip_list[0]->label, tip_list[1]->label, lca->prior_lineno);

    /* set prior on the LCA */
    lca->prior = prior->dist;
    lca->prior_params = prior->params;
    lca->prior_lineno = prior->lineno;
    
    /* output prior info on screen */
    if (!opt_quiet)
    {
      print_prior(tip_list, taxa_count, i, prior);
    }

    /* move to the next prior and free tip_list */
    prior_list = prior_list->next;
    free(tip_list);
  }
}

static void dealloc_prior_list(list_t * list)
{
  if (!list) return;

  free(list->prior->taxa);
  free(list->prior);
  if (list->next)
    dealloc_prior_list(list->next);
  free(list);
}

void set_node_priors(tree_node_t * root, 
                     long * fossils_count, 
                     long * extinct_leaves_count)
{
  *fossils_count = 0;
  *extinct_leaves_count = 0;

  /* create hash table with node labels */
  hashtable_create(root);
    
  /* parse prior file */
  list_t * prior_list = yy_parse_nodeprior(opt_priorfile);

  /* set priors to nodes */
  process_priors(prior_list,
                 root,
                 fossils_count,
                 extinct_leaves_count);

  /* deallocate parsed prior list */
  dealloc_prior_list(prior_list);

  /* destroy hash table */
  hdestroy();
}
