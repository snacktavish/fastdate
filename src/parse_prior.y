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
%{
#include "fastdate.h"

extern int prior_lex();
extern FILE * prior_in;
extern void prior_lex_destroy();
extern int prior_lineno;

void prior_error(list_t ** list, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

static char * concat_del(char * s1, char * s2)
{
  char * new = (char *)malloc(strlen(s1) + strlen(s2) + 2);

  strcpy(new,s1);
  strcat(new,",");
  strcat(new,s2);
  free(s1);
  free(s2);

  return new;
}

static list_t * yy_create_list()
{
  list_t * list = xrealloc(0, sizeof(list_t));
  memset(list,0,sizeof(list_t));
  return list;
}

static void yy_dealloc_list(list_t * list)
{
  if (!list) return;

  free(list->prior->params);
  free(list->prior->taxa);
  free(list->prior);
  if (list->next);
    yy_dealloc_list(list->next);
  free(list);
}
%}

%union
{
  char * s;
  char * d;
  struct list_s * list;
  struct prior_s * prior;
}

%error-verbose
%parse-param {struct list_s ** list}
%destructor { yy_dealloc_list($$); } prior_list

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING EXP LN UNI
%token<d> NUMBER
%type<s> taxa
%type<s> dist
%type<d> number 
%type<list> prior_list
%type<prior> prior
%start input
%%

input: prior_list 
{
  *list = $1;
}

prior_list: prior prior_list
{
  $$ = yy_create_list();
  
  $$->prior = $1;
  $$->next = $2;
}
        | prior
{
  $$ = yy_create_list();

  $$->prior = $1;
  $$->next = NULL;
};

prior: taxa dist OPAR number COMMA number CPAR
{
  $$ = (prior_t *)calloc(1,sizeof(prior_t));
  $$->params = (exp_params_t *)malloc(sizeof(exp_params_t));
  $$->lineno = prior_lineno;

  if (!strcasecmp($2,"exp"))
  {
    $$->dist = NODEPRIOR_EXP;
    ((exp_params_t *)($$->params))->mean = atof($4);
    ((exp_params_t *)($$->params))->offset = atof($6);
  }
  else if (!strcasecmp($2,"uni"))
  {
    $$->dist = NODEPRIOR_UNI;
    ((uni_params_t *)($$->params))->min_age = atof($4);
    ((uni_params_t *)($$->params))->max_age = atof($6);
  }
  else
  {
    assert(0);
  }

  $$->taxa = $1;

  free($2);
  free($4);
  free($6);
}
    | taxa dist OPAR number COMMA number COMMA number CPAR
{
  $$ = (prior_t *)calloc(1,sizeof(prior_t));
  $$->lineno = prior_lineno;
  $$->taxa = $1;

  if (!strcasecmp($2,"ln"))
  {
    /*  setup parasmeters */
    $$->dist = NODEPRIOR_LN;
    $$->params = (ln_params_t *)malloc(sizeof(ln_params_t));
    ((ln_params_t *)($$->params))->mean = atof($4);
    ((ln_params_t *)($$->params))->stdev = atof($6);
    ((ln_params_t *)($$->params))->offset = atof($8);
  }
  else if (!strcasecmp($2,"norm"))
  {
    $$->dist = NODEPRIOR_NORM;
    $$->params = (norm_params_t *)malloc(sizeof(norm_params_t));
    ((norm_params_t *)($$->params))->mean = atof($4);
    ((norm_params_t *)($$->params))->variance = atof($6);
    ((norm_params_t *)($$->params))->offset = atof($8);
  }
  else
  {
    assert(0);
  }

  free($2);
  free($4);
  free($6);
  free($8);
};

taxa: STRING { $$ = $1; }
    | NUMBER { $$ = $1; }
    | STRING COMMA taxa { $$ = concat_del($1,$3); }
    | NUMBER COMMA taxa { $$ = concat_del($1,$3); };

number: NUMBER { $$=$1;};

dist: STRING { $$ = $1; };

%%

list_t * yy_parse_nodeprior(const char * filename)
{
  struct list_s * list; 

  prior_lineno = 1;
  prior_in = fopen(filename, "r");
  if (!prior_in)
  {
    fatal("Cannot open file %s", filename);
  }
  else if (prior_parse(&list))
  {
    fatal("Error while parsing file %s", filename);
  }
  
  if (prior_in) fclose(prior_in);

  prior_lex_destroy();

  return list;
}
