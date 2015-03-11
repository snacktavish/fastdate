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

static int indend_space = 4;

static void output_tree_recursive(tree_node_t * node, FILE * fp_out);

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

  printf (" %s:%f (age: %d)\n", tree->label, tree->length, tree->interval_line);

  if (active_node_order[indend_level-1] == 2) active_node_order[indend_level-1] = 0;

  active_node_order[indend_level] = 1;
  print_tree_recurse(tree->left, indend_level+1, active_node_order);
  active_node_order[indend_level] = 2;
  print_tree_recurse(tree->right, indend_level+1, active_node_order);

}

int tree_indend_level(tree_node_t * tree, int indend)
{
  if (!tree) return indend;

  return MAX(tree_indend_level(tree->left, indend+1),
             tree_indend_level(tree->right, indend+1));

}

void show_ascii_tree(tree_node_t * tree)
{
  int indend_max = tree_indend_level(tree,0);
  int * active_node_order = (int *)malloc((indend_max+1) * sizeof(int));
  active_node_order[0] = 1;
  active_node_order[1] = 1;

  printf (" %s:%f (age: %d)\n", tree->label, tree->length, tree->interval_line);
  print_tree_recurse(tree->left, 1, active_node_order);
  print_tree_recurse(tree->right, 1, active_node_order);
  free(active_node_order);
}

int set_node_heights(tree_node_t * root)
{
  if (!root) return (0);

  root->height = MAX(set_node_heights(root->left), 
                     set_node_heights(root->right)) + 1;

  return (root->height);

}

static void output_tree_recursive(tree_node_t * node, FILE * fp_out)
{
  if (!node->left || !node->right)
    fprintf(fp_out, "%s:[&age=%d]%f", node->label, node->interval_line, node->length);
  else
  {
    fprintf(fp_out, "(");
    output_tree_recursive(node->left, fp_out);
    fprintf(fp_out, ",");
    output_tree_recursive(node->right, fp_out);
    fprintf(fp_out, ")%s:[&age=%d]%f", node->label ? node->label : "", 
                    node->interval_line, node->length);
  }
}

void write_newick_tree(tree_node_t * node)
{
  FILE * fp_out = fopen(opt_outfile, "w");
  if (!fp_out)
    fatal("Unable to open output file for writing");

  output_tree_recursive(node, fp_out);
  fprintf(fp_out, ";");

  fclose(fp_out);
}
