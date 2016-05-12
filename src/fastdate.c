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

static char * progname;
static char progheader[80];
static char * cmdline;

/* number of mandatory options for the user to input */

static const char mand_options_count = 2;
static const char * mand_options_list = " --tree_file\n --out_file\n"
;

/* options */

int opt_quiet;
int opt_outform;
int opt_help;
int opt_version;
int opt_method_relative;
int opt_method_nodeprior;
int opt_method_tipdates;
int opt_showtree;
int opt_fixgamma;
int opt_mu_scale;
char * opt_treefile;
char * opt_outfile;
char * opt_priorfile;
long opt_seed;
long opt_sample;
long opt_grid_intervals;
long opt_threads;
double opt_max_age;
double opt_lambda;
double opt_mu;
double opt_rho;
double opt_psi;
double opt_rate_mean;
double opt_rate_var;
double opt_cred_interval;

unsigned int opt_parameters_bitv;

static struct option long_options[] =
{
  {"help",               no_argument,       0, 0 },  /*  0 */
  {"version",            no_argument,       0, 0 },  /*  1 */
  {"tree_file",          required_argument, 0, 0 },  /*  2 */
  {"show_tree",          no_argument,       0, 0 },  /*  3 */
  {"out_file",           required_argument, 0, 0 },  /*  4 */
  {"out_form",           required_argument, 0, 0 },  /*  5 */
  {"grid",               required_argument, 0, 0 },  /*  6 */
  {"bd_lambda",          required_argument, 0, 0 },  /*  7 */
  {"bd_mu",              required_argument, 0, 0 },  /*  8 */
  {"bd_rho",             required_argument, 0, 0 },  /*  9 */
  {"bd_psi",             required_argument, 0, 0 },  /* 10 */
  {"rate_mean",          required_argument, 0, 0 },  /* 11 */
  {"rate_variance",      required_argument, 0, 0 },  /* 12 */
  {"max_age",            required_argument, 0, 0 },  /* 13 */
  {"method_relative",    no_argument,       0, 0 },  /* 14 */
  {"method_nodeprior",   no_argument,       0, 0 },  /* 15 */
  {"method_tipdates",    no_argument,       0, 0 },  /* 16 */
  {"prior_file",         required_argument, 0, 0 },  /* 17 */
  {"quiet",              no_argument,       0, 0 },  /* 18 */
  {"threads",            required_argument, 0, 0 },  /* 19 */
  {"seed",               required_argument, 0, 0 },  /* 20 */
  {"sample",             required_argument, 0, 0 },  /* 21 */
  {"cred_interval",      required_argument, 0, 0 },  /* 22 */
  {"opt_fix_gamma",      no_argument,       0, 0 },  /* 23 */

  { 0, 0, 0, 0 }
};

void args_init(int argc, char ** argv)
{
  int option_index = 0;
  int c;
  int mand_options = 0;

  /* set defaults */

  progname = argv[0];

  opt_help = 0;
  opt_version = 0;
  opt_method_relative = 0;
  opt_showtree = 0;
  opt_treefile = NULL;
  opt_outfile = NULL;
  opt_priorfile = NULL;
  opt_grid_intervals = 1000;
  opt_quiet = 0;
  opt_max_age = 1;
  opt_threads = 1;
  opt_outform = OUTPUT_ULTRAMETRIC;
  opt_seed = 0;
  opt_sample = 0;
  opt_cred_interval = 0;
  opt_parameters_bitv = 0;
  opt_fixgamma = 0;
  opt_mu_scale = 0;


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
        opt_showtree = 1;
        break;

      case 4:
        free(opt_outfile);
        opt_outfile = optarg;
        break;

      case 5:
        if (strcasecmp(optarg, "ultrametric") == 0)
          opt_outform = OUTPUT_ULTRAMETRIC;
        else if (strcasecmp(optarg, "dated") == 0)
          opt_outform = OUTPUT_DATED;
        else
          fatal("Unrecognized argument for --output-form");
        break;


      case 6:
        opt_grid_intervals = atol(optarg);
        if (opt_grid_intervals < 0)
          fatal("  --grid must be a positive value");
        break;

      case 7:
        opt_lambda = atof(optarg);
        if (opt_lambda < 0)
          fatal("  --bd_lambda must be a positive value");
        break;

      case 8:
        opt_mu = atof(optarg);
        if (opt_mu < 0)
          fatal("  --bd_mu must be a positive value");
        break;

      case 9:
        opt_rho = atof(optarg);
        if (opt_rho < 0 || opt_rho > 1)
          fatal("  --bd_rho is a probability and must be between <0,1>");
        break;

      case 10:
        opt_psi = atof(optarg);
        if (opt_psi <= 0)
          fatal("  --bd_psi must be a positive value");
        break;

      case 11:
        opt_rate_mean = atof(optarg);
        break;

      case 12:
        opt_rate_var  = atof(optarg);
        break;
      case 13:
        opt_max_age = atof(optarg);
        break;

      case 14:
        opt_method_relative = 1;
        break;

      case 15:
        opt_method_nodeprior = 1;
        break;

      case 16:
        opt_method_tipdates = 1;
        break;

      case 17:
        opt_priorfile = optarg;
        break;

      case 18:
        opt_quiet = 1;
        break;

      case 19:
        opt_threads = atol(optarg);
        break;

      case 20:
        opt_seed = atol(optarg);
        break;

      case 21:
        opt_sample = atol(optarg);
        break;

      case 22:
        opt_cred_interval = atof(optarg);
        break;
      case 23:
        opt_fixgamma = 1;
        break;

      default:
        fatal("Internal error in option parsing");
    }
  }

  if (c != -1)
    exit(EXIT_FAILURE);

  int commands  = 0;

  /* check for mandatory options */
  if (opt_treefile)
    mand_options++;
  if (opt_outfile)
    mand_options++;

  /* check for non-fixed parameters */
  if (!(opt_lambda || opt_mu))
  {
    opt_mu_scale = 1;
  }
  if (!opt_lambda)
  {
    opt_lambda =  opt_mu + DEFAULT_LAMBDA;
    opt_parameters_bitv |= PARAM_LAMBDA;
  }
  if (!opt_mu)
  {
    opt_mu = DEFAULT_MU;
    /* check if lambda has been fixed by the user to a very low number */
    if (opt_mu > opt_lambda)
      opt_mu = opt_lambda/10;
    opt_parameters_bitv |= PARAM_MU;
  }
  if (!opt_rho)
  {
    opt_rho = DEFAULT_RHO;
    opt_parameters_bitv |= PARAM_RHO;
  }
  if (opt_fixgamma)
  {
    /* if fixgamma option is set,
     * we optimize both mean and variance together */
    opt_parameters_bitv |= PARAM_RATE_MEAN;
    if (opt_rate_var || opt_rate_mean)
        fatal("  --opt_fix_gamma is incompatible with --rate_mean and --rate_var");
    opt_rate_var = DEFAULT_RATE;
    opt_rate_mean = DEFAULT_RATE;
  }
  else
  {
    if (!opt_rate_mean)
    {
      opt_rate_mean = DEFAULT_RATE;
      opt_parameters_bitv |= PARAM_RATE_MEAN;
    }
    if (!opt_rate_var)
    {
      opt_rate_var = DEFAULT_RATE;
      opt_parameters_bitv |= PARAM_RATE_VAR;
    }
  }

  /* check for number of independent commands selected */
  if (opt_method_relative)
    commands++;
  if (opt_method_nodeprior)
    commands++;
  if (opt_method_tipdates)
    commands++;
  if (opt_version)
    commands++;
  if (opt_help)
    commands++;

  /* if more than one independent command, fail */
  if (commands > 1)
    fatal("More than one command specified");

  /* all method specific checks */
  if (opt_method_relative || opt_method_nodeprior || opt_method_tipdates)
  {
    /* check for mandatory options */
    if (mand_options < mand_options_count)
      fatal("Mandatory options are:\n\n%s",
            mand_options_list);

    if (opt_lambda <= 0)
      fatal("  --bd_lambda must be a positive value");
    if (opt_mu <= 0)
      fatal("  --bd_mu must be a positive value");
    if (opt_lambda <= opt_mu)
      fatal("  --bd_lambda must be greater than --bd_mu");
    if (opt_rho <= 0 || opt_rho > 1)
      fatal("  --bd_rho must be a value between 0 and 1");
  }

  /* --method_relative specific checks */
  if (opt_method_relative)
  {
    if (opt_priorfile)
      fatal("Option --prior_file can only be used with --method_nodeprior or --method_tipdates");
    if (opt_max_age)
      fatal("Option --max_age can only be used with --method_nodeprior "
            "and --method_tipdates");
  }

  /* method_tipdates specific checks */
  if (opt_method_tipdates)
  {
    if (opt_psi == 0)
    {
      opt_psi = DEFAULT_PSI;
      opt_parameters_bitv |= PARAM_PSI;
    }
    if (opt_psi < 0)
      fatal(" Option --bd_psi must be a positive value");
  }

  if (opt_method_nodeprior)
  {
    if (!opt_priorfile)
      fatal("Method --method_nodeprior requires a file with node prior "
            "information to be specified with the --prior_file switch");
    if (!opt_max_age)
      fatal("Method --method_nodeprior requires that --max_age is defined");
  }

  if (opt_method_tipdates)
  {
    if (!opt_priorfile)
      fatal("Method --method_tipdates requires a file with the ages of some "
            "tips to be specified with the --prior_file switch");
    if (!opt_max_age)
      fatal("Method --method_tipdates requires that --max_age is defined");
  }

  /* if no command specified, turn on --help */
  if (!commands) opt_help = 1;

  if ((opt_threads < 0) || (opt_threads > 1024))
    fatal("The argument to --threads must be in the range 0 (default) to 1024");

  if (opt_threads == 0)
    opt_threads = sysconf(_SC_NPROCESSORS_ONLN);

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
          "  --show_tree                    display an ASCII version of the computed tree.\n"
          "  --method_relative              perform divergence time estimations given no fossil information is available.\n"
          "  --method_nodeprior             perform divergence time estimations given prior information is available for certain nodes.\n"
          "  --method_tipdates              perform divergence time estimation given priors on ages of certain tips.\n"
          "  --quiet                        only output warnings and fatal errors to stderr.\n"
          "  --threads INT                  number of threads to use, zero for all cores (default: 0).\n"
          "  --grid INT                     number of grid intervals to use (default: 1000).\n"
          "  --seed INT                     seed for pseudo-random number generator.\n"
          "  --sample INT                   draw INT trees from an approximation of the posterior probability distribution.\n"
          "Input and output options:\n"
          "  --tree_file FILENAME           tree file in newick format.\n"
          "  --prior_file FILENAME          file containing node priors (for --method_nodeprior).\n"
          "  --out_file FILENAME            output file name.\n"
          "  --out_form STRING              format of the output tree. Can be either 'dated' or"
                                            "'ultrametric' (default: 'ultrametric').\n"
          "Model parameters:\n"
          "  --bd_lambda REAL               birth rate for Birth/Death model (default: optimized).\n"
          "  --bd_mu REAL                   death rate for Birth/Death model (default: optimized).\n"
          "  --bd_rho REAL                  probability of sampling individuals. (default: optimized)\n"
          "  --bd_psi REAL                  fossil sample rate\n"
          "  --rate_mean REAL               mean value of edge rate model (default: optimized).\n"
          "  --rate_variance REAL           variance value for edge rate model (default: optimized).\n"
          "  --opt_fix_gamma                assume variance and mean for edge rate model are reciprocal to each other.\n"
          "  --max_age                      max age of the grid when using methods --method_nodeprior or --method_tipdates.\n"
         );
}

void cmd_method_relative()
{
  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");
  tree_node_t * tree = yy_parse_tree(opt_treefile);
  if (!tree)
    fatal("Tree is probably not binary.\n");


  dp(tree);

  if (!opt_quiet)
    fprintf(stdout, "Writing tree file...\n");

  if (opt_sample)
    sample(tree);

  if (opt_cred_interval)
    credible(tree);

  write_newick_tree(tree);

  if (opt_showtree)
    show_ascii_tree(tree);

  if (!opt_quiet)
    fprintf(stdout, "Done\n");

  yy_dealloc_tree(tree);
}

void cmd_method_nodeprior()
{
  long fossils_count, extinct_leaves_count;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");
  tree_node_t * tree = yy_parse_tree(opt_treefile);
  if (!tree)
    fatal("Tree is probably not binary.\n");

  /* set node priors */
  if (!opt_quiet)
    fprintf(stdout, "Setting node priors...\n");
  set_node_priors(tree, &fossils_count, &extinct_leaves_count);
  assert(extinct_leaves_count == 0);
  dp(tree);

  if (!opt_quiet)
    fprintf(stdout, "Writing tree file...\n");

  if (opt_sample)
    sample(tree);

  if (opt_cred_interval)
    credible(tree);

  write_newick_tree(tree);

  if (opt_showtree)
    show_ascii_tree(tree);

  if (!opt_quiet)
    fprintf(stdout, "Done\n");

  yy_dealloc_tree(tree);
}

void cmd_method_tipdates()
{
  long fossils_count, extinct_leaves_count;

  /* parse tree */
  if (!opt_quiet)
    fprintf(stdout, "Parsing tree file...\n");
  tree_node_t * tree = yy_parse_tree(opt_treefile);
  if (!tree)
    fatal("Tree is probably not binary.\n");

  /* set tip priors */
  if (!opt_quiet)
    fprintf(stdout, "Setting tip priors...\n");
  set_node_priors(tree, &fossils_count, &extinct_leaves_count);
  assert(extinct_leaves_count > 0);
  dp(tree);

  if (!opt_quiet)
    fprintf(stdout, "Writing tree file...\n");

  if (opt_sample)
    sample(tree);

  if (opt_cred_interval)
    credible(tree);

  write_newick_tree(tree);

  if (opt_showtree)
    show_ascii_tree(tree);

  if (!opt_quiet)
    fprintf(stdout, "Done\n");

  yy_dealloc_tree(tree);
}

void getentirecommandline(int argc, char * argv[])
{
  int len = 0;
  int i;

  for (i = 0; i < argc; ++i)
    len += strlen(argv[i]);

  cmdline = (char *)xmalloc((size_t)(len + argc + 1));
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
           "%s %s_%s, %1.fGB RAM, %ld cores",
           PROG_NAME, PROG_VERSION, PROG_ARCH,
           arch_get_memtotal() / 1024.0 / 1024.0 / 1024.0,
           sysconf(_SC_NPROCESSORS_ONLN));
}

void show_header()
{
  fprintf(stdout, "%s\n", progheader);
  fprintf(stdout, "https://github.com/xflouris/speed-dating\n");
  fprintf(stdout,"\n");
}

void set_seed()
{
  if (opt_seed)
    srand48(opt_seed);
  else
    srand48(time(NULL));
}

int main (int argc, char * argv[])
{
  fillheader();
  getentirecommandline(argc, argv);

  args_init(argc, argv);

  show_header();

  set_seed();

  if (opt_help)
  {
    cmd_help();
  }
  else if (opt_method_relative)
  {
    cmd_method_relative();
  }
  else if (opt_method_nodeprior)
  {
    cmd_method_nodeprior();
  }
  else if (opt_method_tipdates)
  {
    cmd_method_tipdates();
  }

  free(cmdline);
  return (0);
}
