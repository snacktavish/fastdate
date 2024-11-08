/*
    Copyright (C) 2015 Tomas Flouri, Emily Jane McTavish

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

static int indend_space = 4;

static void output_dated_tree_recursive(tree_node_t * node, FILE * fp_out);
static void output_um_tree_recursive(tree_node_t * node, FILE * fp_out);

static double interval_age = 0;

tree_node_t * yy_create_tree()
{
  tree_node_t * t = xrealloc(0, sizeof(tree_node_t));
  memset(t, 0, sizeof(tree_node_t));
  return t;
}

void yy_dealloc_tree(tree_node_t * tree)
{
  if (!tree) return;
  free(tree->label);
  free(tree->matrix);
  free(tree->matrix_sum);
  free(tree->matrix_left);
  free(tree->matrix_right);
  free(tree->interval_weights);
  free(tree->prior_params);
  yy_dealloc_tree(tree->left);
  yy_dealloc_tree(tree->right);
  free(tree);
}

static void print_tree_recurse(tree_node_t * tree, 
                               int indend_level, int * active_node_order)
{
  int i,j;

  if (!tree) return;

  for (i = 0; i < indend_level; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }
  printf("\n");

  for (i = 0; i < indend_level-1; ++i)
  {
    if (active_node_order[i])
      printf("|");
    else
      printf(" ");

    for (j = 0; j < indend_space-1; ++j)
      printf(" ");
  }

  printf("+");
  for (j = 0; j < indend_space-1; ++j)
    printf ("-");
  if (tree->left || tree->right) printf("+");

  if (!tree->left && !tree->prior)
  {
    printf (" %s\n", tree->label);
  }
  else
  {
    printf (" %s (age: %f", tree->label, tree->interval_line * interval_age);
    if (opt_cred_interval)
      printf (", %.2f CR: [%f,%f]", opt_cred_interval,
                                    tree->cr_minage,
                                    tree->cr_maxage);
    printf(")\n");
  }

  if (active_node_order[indend_level-1] == 2) 
    active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left, indend_level+1, active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right, indend_level+1, active_node_order);

}

static unsigned int tree_indend_level(tree_node_t * tree, unsigned int indend)
{
  if (!tree) return indend;

  unsigned int a, b;

  a = tree_indend_level(tree->left, indend+1);
  b = tree_indend_level(tree->right, indend+1);

  return MAX(a,b);
}

void show_ascii_tree(tree_node_t * tree)
{
  unsigned int indend_max = tree_indend_level(tree,0);
  int * active_node_order = (int *)malloc((indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  printf (" %s (age: %f", tree->label, tree->interval_line * interval_age);
  if (opt_cred_interval)
    printf (", %.2f CR: [%f,%f]", opt_cred_interval,
                                  tree->cr_minage,
                                  tree->cr_maxage);
  printf(")\n");

  print_tree_recurse(tree->left, 1, active_node_order);
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

long set_node_heights(tree_node_t * root)
{
  if (!root) return (0);

  long a,b;

  a = set_node_heights(root->left);
  b = set_node_heights(root->right);

  root->height = MAX(a,b) + 1;
                     
  return (root->height);
}

static void output_dated_tree_recursive(tree_node_t * node, FILE * fp_out)
{
  if (!node->left || !node->right)
  {
    if (!node->prior)
      fprintf(fp_out, "%s:%f", node->label, node->length);
    else
    {
      fprintf(fp_out, "%s[", node->label);
      fprintf(fp_out, "&map_age=%f", node->interval_line * interval_age);
      if (opt_cred_interval)
        fprintf(fp_out, ",lower=%f,upper=%f", node->cr_minage, node->cr_maxage);
      fprintf(fp_out, "]:%f", node->length);
    }
  }
  else
  {
    fprintf(fp_out, "(");
    output_dated_tree_recursive(node->left, fp_out);
    fprintf(fp_out, ",");
    output_dated_tree_recursive(node->right, fp_out);

    fprintf(fp_out, ")%s[", node->label ? node->label : "");
    fprintf(fp_out, "&map_age=%f", node->interval_line * interval_age);
    if (opt_cred_interval)
      fprintf(fp_out, ",lower=%f,upper=%f", node->cr_minage, node->cr_maxage);
    fprintf(fp_out, "]:%f", node->length);
  }
}

static void output_um_tree_recursive(tree_node_t * node, FILE * fp_out)
{
  if (!node->left || !node->right)
  {
    if (!node->prior)
      fprintf(fp_out, "%s:%f", node->label,
              (node->parent->interval_line - node->interval_line) * interval_age);
    else
    {
      fprintf(fp_out, "%s[", node->label);
      fprintf(fp_out, "&map_age=%f", node->interval_line * interval_age);
      if (opt_cred_interval)
        fprintf(fp_out, ",lower=%f,upper=%f", node->cr_minage, node->cr_maxage);
      fprintf(fp_out, "]:%f",
              (node->parent->interval_line - node->interval_line)*interval_age);
    }
  }
  else
  {
    fprintf(fp_out, "(");
    output_um_tree_recursive(node->left, fp_out);
    fprintf(fp_out, ",");
    output_um_tree_recursive(node->right, fp_out);
    fprintf(fp_out, ")%s[", node->label ? node->label : "");
    fprintf(fp_out, "&map_age=%f", node->interval_line * interval_age);
    if (opt_cred_interval)
      fprintf(fp_out, ",lower=%f,upper=%f", node->cr_minage, node->cr_maxage);
    if (node->parent)
      fprintf(fp_out, "]:%f",
              (node->parent->interval_line - node->interval_line)*interval_age);
    else
      fprintf(fp_out, "]:0.0");
  }
}

static void traverse_iterative(tree_node_t * node,
                               int * index,
                               tree_node_t ** outbuffer)
{
  if (!node->left)
  {
    outbuffer[*index] = node;
    *index = *index + 1;
    return;
  }

  traverse_iterative(node->left,  index, outbuffer);
  traverse_iterative(node->right, index, outbuffer);

  outbuffer[*index] = node;
  *index = *index + 1;
}

int tree_traverse(tree_node_t * root, tree_node_t ** outbuffer)
{
  int index = 0;

  traverse_iterative(root, &index, outbuffer);
  return index;
}

static void rtree_query_tipnodes_recursive(tree_node_t * node,
                                           tree_node_t ** node_list,
                                           int * index)
{
  if (!node) return;

  if (!node->left)
  {
    node_list[*index] = node;
    *index = *index + 1;
    return;
  }

  rtree_query_tipnodes_recursive(node->left,  node_list, index);
  rtree_query_tipnodes_recursive(node->right, node_list, index);
}

int rtree_query_tipnodes(tree_node_t * root,
                         tree_node_t ** node_list)
{
  int index = 0;

  if (!root) return 0;
  if (!root->left)
  {
    node_list[index++] = root;
    return index;
  }

  rtree_query_tipnodes_recursive(root->left,  node_list, &index);
  rtree_query_tipnodes_recursive(root->right, node_list, &index);

  return index;
}

void write_newick_tree(tree_node_t * node)
{
  FILE * fp_out = fopen(opt_outfile, "w");

  if (!fp_out)
    fatal("Unable to open output file for writing");

  if (opt_max_age)
    interval_age = opt_max_age / (opt_grid_intervals - 1);
  else
    interval_age = 1.0 / (opt_grid_intervals - 1);
    
  if (opt_outform == OUTPUT_ULTRAMETRIC)
    output_um_tree_recursive(node, fp_out);
  else if (opt_outform == OUTPUT_DATED)
    output_dated_tree_recursive(node, fp_out);
  else
    fatal("Internal error while selecting output format");
    
  fprintf(fp_out, ";\n");

  fclose(fp_out);
}
