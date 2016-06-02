/*
    Copyright (C) 2016 Tomas Flouri

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

static int precision = 7;
static long default_stroke_width = 3;
static long legend_spacing = 10;
static double scaler = 0;
static double max_tiplabel_size = 0;
static double max_tree_len = 0;
static double canvas_width;
static double interval_age = 0;
static double max_age = 0;
static FILE * svg_fp;

static void svg_rect_tt(double x,
                        double y,
                        double width,
                        double height,
                        const char * tooltip_id)
{
  fprintf(svg_fp, "<rect x=\"%f\" y=\"%f\" width=\"%f\" fill=\"#999999\" "
          "height=\"%f\" stroke=\"#999999\" visibility=\"hidden\" "
          "stroke-width=\"1\">\n"
          "<set attributeName=\"visibility\" from=\"hidden\" to=\"visible\" "
          "begin=\"%s.mouseover\" end=\"%s.mouseout\"/></rect>\n",
          x, 
          y,
          width,
          height,
          tooltip_id,
          tooltip_id);
}

static void svg_line(double x1,
                     double y1,
                     double x2,
                     double y2,
                     double stroke_width)
{
  fprintf(svg_fp,
          "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
          "stroke=\"#31a354\" stroke-width=\"%f\" />\n",
          x1, y1, x2, y2, stroke_width);
}

static void svg_circle(double cx, double cy, double r, const char * tooltip_id)
{
  if (!tooltip_id)
    fprintf(svg_fp,
            "<circle cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"#31a354\" "
            "stroke=\"#31a354\" />\n",
            cx, cy, r);
  else
    fprintf(svg_fp,
            "<circle id=\"%s\" cx=\"%f\" cy=\"%f\" r=\"%f\" fill=\"#31a354\" "
            "stroke=\"#31a354\" />\n",
            tooltip_id, cx, cy, r);
}

static void svg_tooltip_ages(double x,
                             double y,
                             double map_age,
                             double cr_minage,
                             double cr_maxage,
                             const char * tooltip_id)
{
  fprintf(svg_fp,
          "<text x=\"%f\" y=\"%f\" font-size=\"12\" font-family=\"Arial;\" "
          "fill=\"black\" visibility=\"hidden\">"
          "map_age: %.3f cr_min_age: %.3f cr_max_age: %.3f\n"
          "<set attributeName=\"visibility\" from=\"hidden\" to=\"visible\" "
          "begin=\"%s.mouseover\" end=\"%s.mouseout\"/></text>\n",
          x, y,
          map_age, cr_minage, cr_maxage,
          tooltip_id, tooltip_id);
}

static double xcoord(long interval_line)
{
  return opt_svg_marginleft +
         (opt_grid_intervals - 1 - interval_line)*interval_age*scaler;
}

static void svg_rtree_plot(tree_node_t * node)
{
  double y;
  static int tip_occ = 0;
  static int tooltip_count = 0;
  char * tooltip_id;

  /* traverse tree in post-order */
  if (node->left)
  {
    svg_rtree_plot(node->left);
    svg_rtree_plot(node->right);
  }

  if (node->parent)
  {
    /* any node that is not the root */

    double x,px;

    x  = xcoord(node->interval_line);
    px = xcoord(node->parent->interval_line);

    if (!node->left)
    {
      y = tip_occ * opt_svg_tipspace + opt_svg_margintop + legend_spacing;
      tip_occ++;
    }
    else
    {
      double ly,ry;

      /* get the y coordinates of its two children */
      ly = node->left->svg_y;
      ry = node->right->svg_y;

      /* set current node y coordinate to the middle of its two childrens */
      y = (ly + ry) / 2.0;

      /* print vertical lines
           
           |
           |
           *
           |
           |
      */
      svg_line(x,
               ly,
               x,
               ry,
               (double)default_stroke_width);
      
      ++tooltip_count;

      if (asprintf(&tooltip_id, "tooltip%d", tooltip_count) == -1)
        fatal("Unable to allocate enough memory.");

      /* print node */
      svg_circle(x,y,opt_svg_inner_radius,tooltip_id);

      /* print tooltip for node ages */
      svg_tooltip_ages(x+10, y+10,
                       node->interval_line*interval_age,
                       node->cr_minage, node->cr_maxage,
                       tooltip_id);

      /* print tooltip rectangle for credible interval */
      svg_rect_tt(opt_svg_marginleft + (max_age - node->cr_maxage)*scaler,
                  y-1,
                  (node->cr_maxage-node->cr_minage)*scaler,
                  2,
                  tooltip_id);

      free(tooltip_id);
    }

    /* print horizontal lines

           +---
           
           *
           
           +---
    */

    svg_line(px,
             y,
             x,
             y,
             (double)default_stroke_width);

    
    /* store the y coordinate of the currrent node */
    node->svg_y = y;

    /* if it is a tip then print its label */
    if (!node->left)
    {
      fprintf(svg_fp,
              "<text x=\"%f\" y=\"%f\" "
              "font-size=\"%ld\" font-family=\"Arial;\">%s</text>\n",
              x+5,
              y+opt_svg_fontsize/3.0,
              opt_svg_fontsize,
              node->label);
    }
    else
      fprintf(svg_fp, "\n");
  }
  else
  {
    /* now do the same as before for the root node */

    double ly,ry,x;

    /* get the y coordinates of the root's two children */
    ly = node->left->svg_y;
    ry = node->right->svg_y;

    /* place root's y coordinate to the middle of its two children */
    y = (ly + ry) / 2.0;

    x = xcoord(node->interval_line);

    /* print vertical lines for root case */
    svg_line(x,
             ly,
             x,
             ry,
             (double)default_stroke_width);

    ++tooltip_count;

    if (asprintf(&tooltip_id, "tooltip%d", tooltip_count) == -1)
      fatal("Unable to allocate enough memory.");

    svg_circle(x,y,opt_svg_inner_radius,tooltip_id);

    svg_tooltip_ages(x+10, y+10,
                     node->interval_line*interval_age,
                     node->cr_minage, node->cr_maxage,
                     tooltip_id);

    svg_rect_tt(opt_svg_marginleft+(max_age - node->cr_maxage)*scaler,
                y-1,
                (node->cr_maxage-node->cr_minage)*scaler,
                2,
                tooltip_id);

    free(tooltip_id);
  }
}

static void rtree_scaler_init(tree_node_t * root)
{
  double length;
  double tiplabel_size;
  unsigned int i;

  /* allocate space for retrieving the tips */
  tree_node_t ** node_list = (tree_node_t **)xmalloc((size_t)(root->leaves)*
                                                     sizeof(tree_node_t *));

  /* retrieve tips */
  rtree_query_tipnodes(root, node_list);

  /* compute max length of tree which in our case (ultrametric trees) is
     equivalent to max_age */
  length = (opt_grid_intervals-1) * interval_age;
  max_tree_len = length;

  /* compute the tip label that takes the most space in the SVG */
  max_tiplabel_size = (opt_svg_fontsize / 1.5) * 
                      (node_list[0]->label ? strlen(node_list[0]->label) : 0);
  for (i = 1; i < root->leaves; ++i)
  {
    tiplabel_size = (opt_svg_fontsize / 1.5) * 
                    (node_list[i]->label ? strlen(node_list[i]->label) : 0);
    if (tiplabel_size > max_tiplabel_size)
      max_tiplabel_size = tiplabel_size;
  }

  /* compute the scaler such that the SVG fits the canvas width */
  scaler = (canvas_width - max_tiplabel_size) / max_tree_len;

  /* deallocate node list */
  free(node_list);
}

static void svg_rtree_init(tree_node_t * root)
{
  int i;
  long svg_height;

  canvas_width = opt_svg_width - opt_svg_marginleft - opt_svg_marginright;

  /* initialize pixel scaler (scaler) and compute max tree 
     length (max_tree_len) */
  rtree_scaler_init(root);

  /* compute the height of the tree in pixels based on opt_svg_tipspace */
  svg_height = opt_svg_margintop + legend_spacing + opt_svg_marginbottom + 
               opt_svg_tipspace * root->leaves;


  /* print svg header tag with dimensions and grey border */
  fprintf(svg_fp, "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"%ld\" "
          "height=\"%ld\" style=\"border: 1px solid #cccccc;\">\n",
          opt_svg_width,
          svg_height);

  /* draw legend */
  if (opt_svg_showlegend)
  {
    /* print the segment with a stroke of 3 */
    svg_line(opt_svg_marginleft,
             opt_svg_margintop / 2,
             (canvas_width - max_tiplabel_size)*opt_svg_legend_ratio +
               opt_svg_marginleft,
             opt_svg_margintop / 2,
             3);

    /* print the length for that segment */
    fprintf(svg_fp, "<text x=\"%f\" y=\"%f\" font-size=\"%ld\" "
            "font-family=\"Arial;\">%.*f</text>\n",
            (canvas_width - max_tiplabel_size)*opt_svg_legend_ratio +
              opt_svg_marginleft + 5,
            opt_svg_margintop/2.0+10-opt_svg_fontsize/3.0,
            (long)opt_svg_fontsize, precision,
            max_tree_len * opt_svg_legend_ratio);
  }

  /* uncomment to print a dashed border to indicate margins */
  
  /*
  fprintf(svg_fp, "<rect x=\"%ld\" y=\"%ld\" width=\"%ld\" fill=\"none\" "
          "height=\"%ld\" stroke=\"#999999\" stroke-dasharray=\"5,5\" "
          "stroke-width=\"1\" />\n",
          opt_svg_marginleft, 
          opt_svg_margintop + legend_spacing, 
          svg_width - opt_svg_marginleft - opt_svg_marginright,
          svg_height - opt_svg_margintop - legend_spacing - opt_svg_marginbottom);
  */

  /* print grid lines */
  for (i = 0; i < opt_grid_intervals; ++i)
  fprintf(svg_fp,
          "<line x1=\"%f\" y1=\"%f\" x2=\"%f\" y2=\"%f\" "
          "stroke=\"#333333\" stroke-width=\"%ld\" title=\"line\" "
          "stroke-dasharray=\"1,1\" />\n",
          (double)(opt_svg_marginleft+i*interval_age*scaler),
          (double)(opt_svg_margintop + legend_spacing),
          (double)(opt_svg_marginleft + i*interval_age*scaler),
          (double)(opt_svg_margintop + (root->leaves-1)*opt_svg_tipspace +
                   legend_spacing),
          (long)1);
          
          
  /* plot the graph */
  svg_rtree_plot(root);

  /* close the svg tag */
  fprintf(svg_fp, "</svg>\n");
}

void cmd_svg(tree_node_t * root)
{
  if (!opt_outfile)
    fatal("An output file must be specified");

  if (!opt_quiet)
    printf("Creating SVG...\n");

  char * filename;
  if (asprintf(&filename, "%s.%s", opt_outfile, "svg") == -1)
    fatal("Unable to allocate enough memory.");

  if (opt_max_age)
    max_age = opt_max_age;
  else
    max_age = 1;

  interval_age = max_age / (opt_grid_intervals - 1);

  svg_fp = fopen(filename, "w");
  if (!svg_fp)
    fatal("Cannot write to file %s", opt_outfile);

  free(filename);
  svg_rtree_init(root);

  fclose(svg_fp);
}
