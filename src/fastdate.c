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

static char * progname;
static char progheader[80];
static char * cmdline;

/* number of mandatory options for the user to input */
static const char mandatory_options_count = 2;
static const char * mandatory_options_list = " --tree-file --out-file";

/* input/output params */
char * opt_treefile;
char * opt_outfile;

/* discretization */
long opt_grid_intervals;

/* birth-death model params */
double opt_birth_rate;
double opt_death_rate;

/* edge rate params */
double opt_edgerate_mean;
double opt_edgerate_var;

/* other commands */
long opt_help;
long opt_version;
long opt_divtimes;
long opt_show_tree;


static struct option long_options[] =
{
  {"help",               no_argument,       0, 0 },  /*  0 */
  {"version",            no_argument,       0, 0 },  /*  1 */
  {"tree-file",          required_argument, 0, 0 },  /*  2 */
  {"show-tree",          no_argument,       0, 0 },  /*  3 */
  {"out-file",           required_argument, 0, 0 },  /*  4 */
  {"grid-intervals",     required_argument, 0, 0 },  /*  5 */
  {"birth-rate",         required_argument, 0, 0 },  /*  6 */      
  {"death-rate",         required_argument, 0, 0 },  /*  7 */
  {"edge-rate-mean",     required_argument, 0, 0 },  /*  8 */
  {"edge-rate-variance", required_argument, 0, 0 },  /*  9 */
  {"divtimes",           no_argument,       0, 0 },  /* 10 */
  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  progname = argv[0];

  opt_help           = 0;
  opt_version        = 0;
  opt_divtimes       = 0;
  opt_show_tree      = 0;

  opt_treefile       = NULL;
  opt_outfile        = NULL;

  /* set some useful defaults here */
  opt_grid_intervals = 0;
  opt_birth_rate     = 0;
  opt_death_rate     = 0;
  opt_edgerate_mean  = 0;
  opt_edgerate_var   = 0;

  while ((c = getopt_long_only(argc, argv, "", long_options, &option_index)) == 0)
  {
    switch (option_index)
    {
      case 0:
        opt_help = 1;
        break;

      case 1:
        opt_version = 1;
        break;

      case 2:
        free(opt_treefile);
        opt_treefile = optarg;
        break;

      case 3:
        opt_show_tree = 1;
        break;

      case 4:
        free(opt_outfile);
        opt_outfile = optarg;
        break;
      
      case 5:
        opt_grid_intervals = atol(optarg);
        break;

      case 6:
        opt_birth_rate = atof(optarg);
        break;

      case 7:
        opt_death_rate = atof(optarg);
        break;
      
      case 8:
        opt_edgerate_mean = atof(optarg);
        break;

      case 9:
        opt_edgerate_var  = atof(optarg);
        break;

      case 10:
        opt_divtimes = 1;
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for --divtimes mandatory options */
  if (opt_treefile)
    mand_options++;
  if (opt_outfile)
    mand_options++;

  /* check for number of independent commands selected */
  if (opt_divtimes)
    commands++;
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");
  
  /* if --divtimes check for mandatory options */
  if (opt_divtimes)
    if (mand_options != mandatory_options_count)
        fatal("Mandatory options for --divtimes are:\n\n%s", 
              mandatory_options_list);


  /* if no command specified, turn on --help */
  if (!commands) opt_help = 1; 
}

void cmd_help()
{
  fprintf(stderr,
          "Usage: %s [OPTIONS]\n", progname);
  fprintf(stderr,
          "\n"
          "General options:\n"
          "  --help                         display help information.\n"
          "  --version                      display version information.\n"
          "  --show-tree                    display an ASCII version of the tree.\n"
          "  --divtimes                     perform divergence time estimations.\n"
          "  --grid-intervals INT           number of grid intervals to use.\n"
          "Input and output options:\n"
          "  --tree-file FILENAME           tree file in newick format.\n"
          "  --out-file FILENAME            output file name.\n"
          "Model parameters:\n"
          "  --birth-rate REAL              birth rate for Birth/Death model.\n"
          "  --death-rate REAL              death rate for Birth/Death model.\n"
          "  --edge-rate-mean REAL          Mean value of edge rate model.\n"
          "  --edge-rate-variance REAL      Variance value for edge rate model.\n"
         );
}

void cmd_divtimes()
{
  /* parse tree */
  fprintf(stdout, "Parsing tree file...\n");
  tree_node_t * tree = yy_parse_tree(opt_treefile);
  if (!tree)
    fatal("Tree is probably not binary.\n");


  dp(tree);

  if (opt_show_tree)
    show_ascii_tree(tree);

  fprintf(stdout, "Writing tree file...\n");
  write_newick_tree(tree);
  fprintf(stdout, "Done\n");

  yy_dealloc_tree(tree);

}


void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc(len + argc + 1);
  cmdline[0] = 0;

  for (i = 0; i < argc; ++i)
  {
    strcat(cmdline, argv[i]);
    strcat(cmdline, " ");
  }
}

void fillheader()
{
  snprintf(progheader, 80,
           "%s %s_%s",
           PROG_NAME, PROG_VERSION, PROG_ARCH);
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/speed-dating\n");
  fprintf(stdout,"\n");
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_divtimes)
  {
    cmd_divtimes();
  }

  free(cmdline);
  return (0);
}
