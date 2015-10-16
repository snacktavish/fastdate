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
#include "parse_newick.h"
#include "fastdate.h"

static size_t string_length = 0;

static char * append(size_t * current_len, const char * s2, size_t append_len)
{
  newick_lval.s = (char *)xrealloc(newick_lval.s, sizeof(char) * (*current_len+append_len+1));
  strncpy(&newick_lval.s[*current_len], s2, append_len);
  newick_lval.s[*current_len + append_len]=0;
  *current_len += append_len;
  return newick_lval.s;
}

%}
%option noyywrap
%option nounput
%option noinput
%s apos
%s quot

%%
<quot>{
\\\"                                         append(&string_length, "\"", 1);
\'                                           append(&string_length, "\'", 1);
\"                                           {BEGIN(INITIAL);return STRING;}
}

<apos>{
\\\'                                         append(&string_length, "\'", 1);
\"                                           append(&string_length, "\"", 1);
\'                                           {BEGIN(INITIAL);return STRING;}
}

<apos,quot>{
\\n                                          append(&string_length, "\n", 1);
\\t                                          append(&string_length, "\t", 1);
\\\\                                         append(&string_length, "\\", 1);
([^\"\'\\]|\n)+                              append(&string_length, newick_text, newick_leng);
}

\:                                           return COLON;
\;                                           return SEMICOLON;
\)                                           return CPAR;
\(                                           return OPAR;
\,                                           return COMMA;
\"                                           {string_length =0; BEGIN(quot);}
\'                                           {string_length =0; BEGIN(apos);}
[a-zA-Z_][a-zA-Z_0-9]*                       {newick_lval.s = xstrndup(newick_text, newick_leng); return STRING;}
[\+|\-]?[0-9]+                               {newick_lval.d = xstrndup(newick_text, newick_leng); return NUMBER;}
[\+|\-]?[0-9]+\.[0-9]+([e|E][\+|\-]?[0-9]+)? {newick_lval.d = xstrndup(newick_text, newick_leng); return NUMBER;}
[ \t\n\r]                                    ;
.                                            {fprintf(stderr, "Syntax error (%c)\n", newick_text[0]); exit(EXIT_FAILURE);}
%%
