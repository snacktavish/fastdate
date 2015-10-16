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

/* fill path with nodes of the path tip to root */
static void fill_path(tree_node_t ** path, int * path_len, tree_node_t * tip)
{
  int i = 0;

  while (tip)
  {
    path[i++] = tip;
    tip = tip->parent;
  }

  *path_len = i;
}

tree_node_t * lca_compute(tree_node_t * root,
                          tree_node_t ** tip_nodes,
                          unsigned int count)
{
  unsigned int i;
  tree_node_t *** path;

  assert(count >= 2);

  path = (tree_node_t ***)xmalloc((size_t)count *
                                  sizeof(tree_node_t **));
  int * path_len = (int *)xmalloc((size_t)count * sizeof(int));

  /* fill paths */
  for (i = 0; i < count; ++i)
  {
    path[i] = (tree_node_t **)xmalloc((size_t)(root->height+1) *
                                      sizeof(tree_node_t *));
    
    fill_path(path[i], &(path_len[i]), tip_nodes[i]);
  }

  /* find LCA */
  tree_node_t * lca = NULL;
  while (!lca)
  {
    for (i = 0; i < count; ++i)
      --path_len[i];

    for (i = 1; i < count; ++i)
    {
      if (path[i-1][path_len[i-1]] != path[i][path_len[i]])
      {
        lca = path[i][path_len[i]+1];
        break;
      }
    }
  }

  /* free allocated memory */
  for (i = 0; i < count; ++i)
    free(path[i]);
  free(path);
  free(path_len);

  return lca;
}
