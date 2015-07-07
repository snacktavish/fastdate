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

#define HASHTABLE_SIZE 5000

extern int yylex();
extern FILE * yyin;
extern void yylex_destroy();

void yyerror(tree_node_t * tree, const char * s) 
{
  fprintf(stderr, "%s.\n", s);
}

static void hash_table_create()
{
  hcreate(HASHTABLE_SIZE);
}

/* insert the taxon/node label in the hashtable */
static void hash_table_insert(char * label)
{
  ENTRY entry;

  entry.key = label;
  entry.data = NULL;
  hsearch(entry,ENTER);

}

/* check if a taxon/node label with thte same name already exists */
static void hash_table_find(char * label)
{
  ENTRY entry;

  entry.key = label;
  if (hsearch(entry,FIND))
    fatal("Duplicate taxon or node label (%s) in file %s", label, opt_treefile);
}

static void hash_table_destroy()
{
  hdestroy();
}

%}


%union
{
  char * s;
  char * d;
  struct tree_noderec * tree;
}

%error-verbose
%parse-param {struct tree_noderec * tree}
%destructor { yy_dealloc_tree($$); } subtree

%token OPAR
%token CPAR
%token COMMA
%token COLON SEMICOLON 
%token<s> STRING
%token<d> NUMBER
%type<s> label optional_label
%type<d> number optional_length
%type<tree> subtree
%start input
%%

input: OPAR subtree COMMA subtree CPAR optional_label optional_length SEMICOLON
{
  //tree_node_t * tree = yy_create_tree();
  tree->left   = $2;
  tree->right  = $4;
  tree->label  = $6;
  tree->length = $7 ? atof($7) : 0;
  tree->leaves = $2->leaves + $4->leaves;
  tree->height = ($2->height > $4->height) ? $2->height + 1 : $4->height + 1;
  tree->parent = NULL;

  tree->matrix_left  = NULL;
  tree->matrix_right = NULL;
  tree->matrix       = NULL;

  $2->parent = tree;
  $4->parent = tree;

  free($7);

  /* check for duplicate taxa/node labels */
  if (tree->label)
  {
    hash_table_find(tree->label);
    hash_table_insert(tree->label);
  }
};

subtree: OPAR subtree COMMA subtree CPAR optional_label optional_length
{
  $$ = yy_create_tree();
  $$->left   = $2;
  $$->right  = $4;
  $$->label  = $6;
  $$->length = $7 ? atof($7) : 0;
  $$->leaves = $2->leaves + $4->leaves;
  $$->height = ($2->height > $4->height) ? $2->height + 1 : $4->height + 1;

  $$->matrix_left  = NULL;
  $$->matrix_right = NULL;
  $$->matrix       = NULL;
  
  $2->parent = $$;
  $4->parent = $$;

  free($7);

  /* check for duplicate taxa/node labels */
  if ($$->label)
  {
    hash_table_find($$->label);
    hash_table_insert($$->label);
  }

}
       | label optional_length
{
  $$ = yy_create_tree();
  $$->label  = $1;
  $$->length = $2 ? atof($2) : 0;
  $$->left   = NULL;
  $$->right  = NULL;
  $$->leaves = 1;
  $$->height = 0;

  $$->matrix_left  = NULL;
  $$->matrix_right = NULL;
  $$->matrix       = NULL;

  free($2);
  
  /* check for duplicate taxa/node labels */
  hash_table_find($$->label);
  hash_table_insert($$->label);
};

 
optional_label:  { $$ = NULL;} | label  {$$ = $1;};
optional_length: { $$ = NULL;} | COLON number {$$ = $2;};
label: STRING    { $$=$1;};
number: NUMBER   { $$=$1;};

%%

tree_node_t * yy_parse_tree(const char * filename)
{
  struct tree_noderec * tree;

  tree = yy_create_tree();
  hash_table_create();

  yyin = fopen(filename, "r");
  if (!yyin)
  {
    fatal("Cannot open file %s", filename);
  }
  else if (yyparse(tree))
  {
    fatal("Cannot parse tree file %s (maybe non-binary?)", filename);
  }
  
  if (yyin) fclose(yyin);

  yylex_destroy();
  hash_table_destroy();

  return tree;
}
