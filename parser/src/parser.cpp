// A Bison parser, made by GNU Bison 3.5.1.

// Skeleton implementation for Bison LALR(1) parsers in C++

// Copyright (C) 2002-2015, 2018-2020 Free Software Foundation, Inc.

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// As a special exception, you may create a larger work that contains
// part or all of the Bison parser skeleton and distribute that work
// under terms of your choice, so long as that work isn't itself a
// parser generator using the skeleton or a modified version thereof
// as a parser skeleton.  Alternatively, if you modify or redistribute
// the parser skeleton itself, you may (at your option) remove this
// special exception, which will cause the skeleton and the resulting
// Bison output files to be licensed under the GNU General Public
// License without this special exception.

// This special exception was added by the Free Software Foundation in
// version 2.2 of Bison.

// Undocumented macros, especially those whose name start with YY_,
// are private implementation details.  Do not rely on them.





#include "parser.h"


// Unqualified %code blocks.
#line 24 "parser/parser.yy"

#include "driver.h"

#line 49 "parser/src/parser.cpp"


#ifndef YY_
# if defined YYENABLE_NLS && YYENABLE_NLS
#  if ENABLE_NLS
#   include <libintl.h> // FIXME: INFRINGES ON USER NAME SPACE.
#   define YY_(msgid) dgettext ("bison-runtime", msgid)
#  endif
# endif
# ifndef YY_
#  define YY_(msgid) msgid
# endif
#endif

// Whether we are compiled with exception support.
#ifndef YY_EXCEPTIONS
# if defined __GNUC__ && !defined __EXCEPTIONS
#  define YY_EXCEPTIONS 0
# else
#  define YY_EXCEPTIONS 1
# endif
#endif

#define YYRHSLOC(Rhs, K) ((Rhs)[K].location)
/* YYLLOC_DEFAULT -- Set CURRENT to span from RHS[1] to RHS[N].
   If N is 0, then set CURRENT to the empty location which ends
   the previous symbol: RHS[0] (always defined).  */

# ifndef YYLLOC_DEFAULT
#  define YYLLOC_DEFAULT(Current, Rhs, N)                               \
    do                                                                  \
      if (N)                                                            \
        {                                                               \
          (Current).begin  = YYRHSLOC (Rhs, 1).begin;                   \
          (Current).end    = YYRHSLOC (Rhs, N).end;                     \
        }                                                               \
      else                                                              \
        {                                                               \
          (Current).begin = (Current).end = YYRHSLOC (Rhs, 0).end;      \
        }                                                               \
    while (false)
# endif


// Enable debugging if requested.
#if YYDEBUG

// A pseudo ostream that takes yydebug_ into account.
# define YYCDEBUG if (yydebug_) (*yycdebug_)

# define YY_SYMBOL_PRINT(Title, Symbol)         \
  do {                                          \
    if (yydebug_)                               \
    {                                           \
      *yycdebug_ << Title << ' ';               \
      yy_print_ (*yycdebug_, Symbol);           \
      *yycdebug_ << '\n';                       \
    }                                           \
  } while (false)

# define YY_REDUCE_PRINT(Rule)          \
  do {                                  \
    if (yydebug_)                       \
      yy_reduce_print_ (Rule);          \
  } while (false)

# define YY_STACK_PRINT()               \
  do {                                  \
    if (yydebug_)                       \
      yystack_print_ ();                \
  } while (false)

#else // !YYDEBUG

# define YYCDEBUG if (false) std::cerr
# define YY_SYMBOL_PRINT(Title, Symbol)  YYUSE (Symbol)
# define YY_REDUCE_PRINT(Rule)           static_cast<void> (0)
# define YY_STACK_PRINT()                static_cast<void> (0)

#endif // !YYDEBUG

#define yyerrok         (yyerrstatus_ = 0)
#define yyclearin       (yyla.clear ())

#define YYACCEPT        goto yyacceptlab
#define YYABORT         goto yyabortlab
#define YYERROR         goto yyerrorlab
#define YYRECOVERING()  (!!yyerrstatus_)

namespace yy {
#line 140 "parser/src/parser.cpp"


  /* Return YYSTR after stripping away unnecessary quotes and
     backslashes, so that it's suitable for yyerror.  The heuristic is
     that double-quoting is unnecessary unless the string contains an
     apostrophe, a comma, or backslash (other than backslash-backslash).
     YYSTR is taken from yytname.  */
  std::string
  parser::yytnamerr_ (const char *yystr)
  {
    if (*yystr == '"')
      {
        std::string yyr;
        char const *yyp = yystr;

        for (;;)
          switch (*++yyp)
            {
            case '\'':
            case ',':
              goto do_not_strip_quotes;

            case '\\':
              if (*++yyp != '\\')
                goto do_not_strip_quotes;
              else
                goto append;

            append:
            default:
              yyr += *yyp;
              break;

            case '"':
              return yyr;
            }
      do_not_strip_quotes: ;
      }

    return yystr;
  }


  /// Build a parser object.
  parser::parser (driver& drv_yyarg)
#if YYDEBUG
    : yydebug_ (false),
      yycdebug_ (&std::cerr),
#else
    :
#endif
      drv (drv_yyarg)
  {}

  parser::~parser ()
  {}

  parser::syntax_error::~syntax_error () YY_NOEXCEPT YY_NOTHROW
  {}

  /*---------------.
  | Symbol types.  |
  `---------------*/



  // by_state.
  parser::by_state::by_state () YY_NOEXCEPT
    : state (empty_state)
  {}

  parser::by_state::by_state (const by_state& that) YY_NOEXCEPT
    : state (that.state)
  {}

  void
  parser::by_state::clear () YY_NOEXCEPT
  {
    state = empty_state;
  }

  void
  parser::by_state::move (by_state& that)
  {
    state = that.state;
    that.clear ();
  }

  parser::by_state::by_state (state_type s) YY_NOEXCEPT
    : state (s)
  {}

  parser::symbol_number_type
  parser::by_state::type_get () const YY_NOEXCEPT
  {
    if (state == empty_state)
      return empty_symbol;
    else
      return yystos_[+state];
  }

  parser::stack_symbol_type::stack_symbol_type ()
  {}

  parser::stack_symbol_type::stack_symbol_type (YY_RVREF (stack_symbol_type) that)
    : super_type (YY_MOVE (that.state), YY_MOVE (that.location))
  {
    switch (that.type_get ())
    {
      case 61: // expr
        value.YY_MOVE_OR_COPY< AbsSyn::Expr * > (YY_MOVE (that.value));
        break;

      case 67: // path_formula
      case 68: // state_formula
        value.YY_MOVE_OR_COPY< AbsSyn::Formula * > (YY_MOVE (that.value));
        break;

      case 56: // modeType
        value.YY_MOVE_OR_COPY< AbsSyn::modeType > (YY_MOVE (that.value));
        break;

      case 55: // problemType
        value.YY_MOVE_OR_COPY< AbsSyn::problemType > (YY_MOVE (that.value));
        break;

      case 69: // transType
        value.YY_MOVE_OR_COPY< AbsSyn::transType > (YY_MOVE (that.value));
        break;

      case 52: // DOUBLE
      case 59: // number
        value.YY_MOVE_OR_COPY< double > (YY_MOVE (that.value));
        break;

      case 51: // INTEGER
        value.YY_MOVE_OR_COPY< int > (YY_MOVE (that.value));
        break;

      case 57: // doubleInterval
        value.YY_MOVE_OR_COPY< std::pair<double, double> > (YY_MOVE (that.value));
        break;

      case 58: // intInterval
        value.YY_MOVE_OR_COPY< std::pair<int, int> > (YY_MOVE (that.value));
        break;

      case 50: // IDENT
        value.YY_MOVE_OR_COPY< std::string > (YY_MOVE (that.value));
        break;

      case 62: // numList
        value.YY_MOVE_OR_COPY< std::vector<double> > (YY_MOVE (that.value));
        break;

      case 63: // matrixRow
      case 64: // intList
        value.YY_MOVE_OR_COPY< std::vector<int> > (YY_MOVE (that.value));
        break;

      case 60: // identList
        value.YY_MOVE_OR_COPY< std::vector<std::string> > (YY_MOVE (that.value));
        break;

      case 65: // _matrix
      case 66: // rowList
        value.YY_MOVE_OR_COPY< std::vector<std::vector<int>> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

#if 201103L <= YY_CPLUSPLUS
    // that is emptied.
    that.state = empty_state;
#endif
  }

  parser::stack_symbol_type::stack_symbol_type (state_type s, YY_MOVE_REF (symbol_type) that)
    : super_type (s, YY_MOVE (that.location))
  {
    switch (that.type_get ())
    {
      case 61: // expr
        value.move< AbsSyn::Expr * > (YY_MOVE (that.value));
        break;

      case 67: // path_formula
      case 68: // state_formula
        value.move< AbsSyn::Formula * > (YY_MOVE (that.value));
        break;

      case 56: // modeType
        value.move< AbsSyn::modeType > (YY_MOVE (that.value));
        break;

      case 55: // problemType
        value.move< AbsSyn::problemType > (YY_MOVE (that.value));
        break;

      case 69: // transType
        value.move< AbsSyn::transType > (YY_MOVE (that.value));
        break;

      case 52: // DOUBLE
      case 59: // number
        value.move< double > (YY_MOVE (that.value));
        break;

      case 51: // INTEGER
        value.move< int > (YY_MOVE (that.value));
        break;

      case 57: // doubleInterval
        value.move< std::pair<double, double> > (YY_MOVE (that.value));
        break;

      case 58: // intInterval
        value.move< std::pair<int, int> > (YY_MOVE (that.value));
        break;

      case 50: // IDENT
        value.move< std::string > (YY_MOVE (that.value));
        break;

      case 62: // numList
        value.move< std::vector<double> > (YY_MOVE (that.value));
        break;

      case 63: // matrixRow
      case 64: // intList
        value.move< std::vector<int> > (YY_MOVE (that.value));
        break;

      case 60: // identList
        value.move< std::vector<std::string> > (YY_MOVE (that.value));
        break;

      case 65: // _matrix
      case 66: // rowList
        value.move< std::vector<std::vector<int>> > (YY_MOVE (that.value));
        break;

      default:
        break;
    }

    // that is emptied.
    that.type = empty_symbol;
  }

#if YY_CPLUSPLUS < 201103L
  parser::stack_symbol_type&
  parser::stack_symbol_type::operator= (const stack_symbol_type& that)
  {
    state = that.state;
    switch (that.type_get ())
    {
      case 61: // expr
        value.copy< AbsSyn::Expr * > (that.value);
        break;

      case 67: // path_formula
      case 68: // state_formula
        value.copy< AbsSyn::Formula * > (that.value);
        break;

      case 56: // modeType
        value.copy< AbsSyn::modeType > (that.value);
        break;

      case 55: // problemType
        value.copy< AbsSyn::problemType > (that.value);
        break;

      case 69: // transType
        value.copy< AbsSyn::transType > (that.value);
        break;

      case 52: // DOUBLE
      case 59: // number
        value.copy< double > (that.value);
        break;

      case 51: // INTEGER
        value.copy< int > (that.value);
        break;

      case 57: // doubleInterval
        value.copy< std::pair<double, double> > (that.value);
        break;

      case 58: // intInterval
        value.copy< std::pair<int, int> > (that.value);
        break;

      case 50: // IDENT
        value.copy< std::string > (that.value);
        break;

      case 62: // numList
        value.copy< std::vector<double> > (that.value);
        break;

      case 63: // matrixRow
      case 64: // intList
        value.copy< std::vector<int> > (that.value);
        break;

      case 60: // identList
        value.copy< std::vector<std::string> > (that.value);
        break;

      case 65: // _matrix
      case 66: // rowList
        value.copy< std::vector<std::vector<int>> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    return *this;
  }

  parser::stack_symbol_type&
  parser::stack_symbol_type::operator= (stack_symbol_type& that)
  {
    state = that.state;
    switch (that.type_get ())
    {
      case 61: // expr
        value.move< AbsSyn::Expr * > (that.value);
        break;

      case 67: // path_formula
      case 68: // state_formula
        value.move< AbsSyn::Formula * > (that.value);
        break;

      case 56: // modeType
        value.move< AbsSyn::modeType > (that.value);
        break;

      case 55: // problemType
        value.move< AbsSyn::problemType > (that.value);
        break;

      case 69: // transType
        value.move< AbsSyn::transType > (that.value);
        break;

      case 52: // DOUBLE
      case 59: // number
        value.move< double > (that.value);
        break;

      case 51: // INTEGER
        value.move< int > (that.value);
        break;

      case 57: // doubleInterval
        value.move< std::pair<double, double> > (that.value);
        break;

      case 58: // intInterval
        value.move< std::pair<int, int> > (that.value);
        break;

      case 50: // IDENT
        value.move< std::string > (that.value);
        break;

      case 62: // numList
        value.move< std::vector<double> > (that.value);
        break;

      case 63: // matrixRow
      case 64: // intList
        value.move< std::vector<int> > (that.value);
        break;

      case 60: // identList
        value.move< std::vector<std::string> > (that.value);
        break;

      case 65: // _matrix
      case 66: // rowList
        value.move< std::vector<std::vector<int>> > (that.value);
        break;

      default:
        break;
    }

    location = that.location;
    // that is emptied.
    that.state = empty_state;
    return *this;
  }
#endif

  template <typename Base>
  void
  parser::yy_destroy_ (const char* yymsg, basic_symbol<Base>& yysym) const
  {
    if (yymsg)
      YY_SYMBOL_PRINT (yymsg, yysym);
  }

#if YYDEBUG
  template <typename Base>
  void
  parser::yy_print_ (std::ostream& yyo,
                                     const basic_symbol<Base>& yysym) const
  {
    std::ostream& yyoutput = yyo;
    YYUSE (yyoutput);
    symbol_number_type yytype = yysym.type_get ();
#if defined __GNUC__ && ! defined __clang__ && ! defined __ICC && __GNUC__ * 100 + __GNUC_MINOR__ <= 408
    // Avoid a (spurious) G++ 4.8 warning about "array subscript is
    // below array bounds".
    if (yysym.empty ())
      std::abort ();
#endif
    yyo << (yytype < yyntokens_ ? "token" : "nterm")
        << ' ' << yytname_[yytype] << " ("
        << yysym.location << ": ";
    switch (yytype)
    {
      case 50: // IDENT
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::string > (); }
#line 576 "parser/src/parser.cpp"
        break;

      case 51: // INTEGER
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < int > (); }
#line 582 "parser/src/parser.cpp"
        break;

      case 52: // DOUBLE
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < double > (); }
#line 588 "parser/src/parser.cpp"
        break;

      case 55: // problemType
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::problemType > (); }
#line 594 "parser/src/parser.cpp"
        break;

      case 56: // modeType
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::modeType > (); }
#line 600 "parser/src/parser.cpp"
        break;

      case 57: // doubleInterval
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::pair<double, double> > (); }
#line 606 "parser/src/parser.cpp"
        break;

      case 58: // intInterval
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::pair<int, int> > (); }
#line 612 "parser/src/parser.cpp"
        break;

      case 59: // number
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < double > (); }
#line 618 "parser/src/parser.cpp"
        break;

      case 60: // identList
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<std::string> > (); }
#line 624 "parser/src/parser.cpp"
        break;

      case 61: // expr
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::Expr * > (); }
#line 630 "parser/src/parser.cpp"
        break;

      case 62: // numList
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<double> > (); }
#line 636 "parser/src/parser.cpp"
        break;

      case 63: // matrixRow
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<int> > (); }
#line 642 "parser/src/parser.cpp"
        break;

      case 64: // intList
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<int> > (); }
#line 648 "parser/src/parser.cpp"
        break;

      case 65: // _matrix
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<std::vector<int>> > (); }
#line 654 "parser/src/parser.cpp"
        break;

      case 66: // rowList
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < std::vector<std::vector<int>> > (); }
#line 660 "parser/src/parser.cpp"
        break;

      case 67: // path_formula
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::Formula * > (); }
#line 666 "parser/src/parser.cpp"
        break;

      case 68: // state_formula
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::Formula * > (); }
#line 672 "parser/src/parser.cpp"
        break;

      case 69: // transType
#line 108 "parser/parser.yy"
                 { yyo << yysym.value.template as < AbsSyn::transType > (); }
#line 678 "parser/src/parser.cpp"
        break;

      default:
        break;
    }
    yyo << ')';
  }
#endif

  void
  parser::yypush_ (const char* m, YY_MOVE_REF (stack_symbol_type) sym)
  {
    if (m)
      YY_SYMBOL_PRINT (m, sym);
    yystack_.push (YY_MOVE (sym));
  }

  void
  parser::yypush_ (const char* m, state_type s, YY_MOVE_REF (symbol_type) sym)
  {
#if 201103L <= YY_CPLUSPLUS
    yypush_ (m, stack_symbol_type (s, std::move (sym)));
#else
    stack_symbol_type ss (s, sym);
    yypush_ (m, ss);
#endif
  }

  void
  parser::yypop_ (int n)
  {
    yystack_.pop (n);
  }

#if YYDEBUG
  std::ostream&
  parser::debug_stream () const
  {
    return *yycdebug_;
  }

  void
  parser::set_debug_stream (std::ostream& o)
  {
    yycdebug_ = &o;
  }


  parser::debug_level_type
  parser::debug_level () const
  {
    return yydebug_;
  }

  void
  parser::set_debug_level (debug_level_type l)
  {
    yydebug_ = l;
  }
#endif // YYDEBUG

  parser::state_type
  parser::yy_lr_goto_state_ (state_type yystate, int yysym)
  {
    int yyr = yypgoto_[yysym - yyntokens_] + yystate;
    if (0 <= yyr && yyr <= yylast_ && yycheck_[yyr] == yystate)
      return yytable_[yyr];
    else
      return yydefgoto_[yysym - yyntokens_];
  }

  bool
  parser::yy_pact_value_is_default_ (int yyvalue)
  {
    return yyvalue == yypact_ninf_;
  }

  bool
  parser::yy_table_value_is_error_ (int yyvalue)
  {
    return yyvalue == yytable_ninf_;
  }

  int
  parser::operator() ()
  {
    return parse ();
  }

  int
  parser::parse ()
  {
    int yyn;
    /// Length of the RHS of the rule being reduced.
    int yylen = 0;

    // Error handling.
    int yynerrs_ = 0;
    int yyerrstatus_ = 0;

    /// The lookahead symbol.
    symbol_type yyla;

    /// The locations where the error started and ended.
    stack_symbol_type yyerror_range[3];

    /// The return value of parse ().
    int yyresult;

#if YY_EXCEPTIONS
    try
#endif // YY_EXCEPTIONS
      {
    YYCDEBUG << "Starting parse\n";


    /* Initialize the stack.  The initial state will be set in
       yynewstate, since the latter expects the semantical and the
       location values to have been already stored, initialize these
       stacks with a primary value.  */
    yystack_.clear ();
    yypush_ (YY_NULLPTR, 0, YY_MOVE (yyla));

  /*-----------------------------------------------.
  | yynewstate -- push a new symbol on the stack.  |
  `-----------------------------------------------*/
  yynewstate:
    YYCDEBUG << "Entering state " << int (yystack_[0].state) << '\n';

    // Accept?
    if (yystack_[0].state == yyfinal_)
      YYACCEPT;

    goto yybackup;


  /*-----------.
  | yybackup.  |
  `-----------*/
  yybackup:
    // Try to take a decision without lookahead.
    yyn = yypact_[+yystack_[0].state];
    if (yy_pact_value_is_default_ (yyn))
      goto yydefault;

    // Read a lookahead token.
    if (yyla.empty ())
      {
        YYCDEBUG << "Reading a token: ";
#if YY_EXCEPTIONS
        try
#endif // YY_EXCEPTIONS
          {
            symbol_type yylookahead (yylex (drv));
            yyla.move (yylookahead);
          }
#if YY_EXCEPTIONS
        catch (const syntax_error& yyexc)
          {
            YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
            error (yyexc);
            goto yyerrlab1;
          }
#endif // YY_EXCEPTIONS
      }
    YY_SYMBOL_PRINT ("Next token is", yyla);

    /* If the proper action on seeing token YYLA.TYPE is to reduce or
       to detect an error, take that action.  */
    yyn += yyla.type_get ();
    if (yyn < 0 || yylast_ < yyn || yycheck_[yyn] != yyla.type_get ())
      {
        goto yydefault;
      }

    // Reduce or error.
    yyn = yytable_[yyn];
    if (yyn <= 0)
      {
        if (yy_table_value_is_error_ (yyn))
          goto yyerrlab;
        yyn = -yyn;
        goto yyreduce;
      }

    // Count tokens shifted since error; after three, turn off error status.
    if (yyerrstatus_)
      --yyerrstatus_;

    // Shift the lookahead token.
    yypush_ ("Shifting", state_type (yyn), YY_MOVE (yyla));
    goto yynewstate;


  /*-----------------------------------------------------------.
  | yydefault -- do the default action for the current state.  |
  `-----------------------------------------------------------*/
  yydefault:
    yyn = yydefact_[+yystack_[0].state];
    if (yyn == 0)
      goto yyerrlab;
    goto yyreduce;


  /*-----------------------------.
  | yyreduce -- do a reduction.  |
  `-----------------------------*/
  yyreduce:
    yylen = yyr2_[yyn];
    {
      stack_symbol_type yylhs;
      yylhs.state = yy_lr_goto_state_ (yystack_[yylen].state, yyr1_[yyn]);
      /* Variants are always initialized to an empty instance of the
         correct type. The default '$$ = $1' action is NOT applied
         when using variants.  */
      switch (yyr1_[yyn])
    {
      case 61: // expr
        yylhs.value.emplace< AbsSyn::Expr * > ();
        break;

      case 67: // path_formula
      case 68: // state_formula
        yylhs.value.emplace< AbsSyn::Formula * > ();
        break;

      case 56: // modeType
        yylhs.value.emplace< AbsSyn::modeType > ();
        break;

      case 55: // problemType
        yylhs.value.emplace< AbsSyn::problemType > ();
        break;

      case 69: // transType
        yylhs.value.emplace< AbsSyn::transType > ();
        break;

      case 52: // DOUBLE
      case 59: // number
        yylhs.value.emplace< double > ();
        break;

      case 51: // INTEGER
        yylhs.value.emplace< int > ();
        break;

      case 57: // doubleInterval
        yylhs.value.emplace< std::pair<double, double> > ();
        break;

      case 58: // intInterval
        yylhs.value.emplace< std::pair<int, int> > ();
        break;

      case 50: // IDENT
        yylhs.value.emplace< std::string > ();
        break;

      case 62: // numList
        yylhs.value.emplace< std::vector<double> > ();
        break;

      case 63: // matrixRow
      case 64: // intList
        yylhs.value.emplace< std::vector<int> > ();
        break;

      case 60: // identList
        yylhs.value.emplace< std::vector<std::string> > ();
        break;

      case 65: // _matrix
      case 66: // rowList
        yylhs.value.emplace< std::vector<std::vector<int>> > ();
        break;

      default:
        break;
    }


      // Default location.
      {
        stack_type::slice range (yystack_, yylen);
        YYLLOC_DEFAULT (yylhs.location, range, yylen);
        yyerror_range[1].location = yylhs.location;
      }

      // Perform the reduction.
      YY_REDUCE_PRINT (yyn);
#if YY_EXCEPTIONS
      try
#endif // YY_EXCEPTIONS
        {
          switch (yyn)
            {
  case 2:
#line 115 "parser/parser.yy"
                                {
					if (!drv.m->isProblemDefined())
					{
						yy::parser::error(yystack_[0].location, "Problem type must be defined");
						YYERROR;
					}
					
					if (!drv.m->isVarModeDefined())
						drv.m->setVarMode(AbsSyn::modeType::BOX);
					
					if (!drv.m->isParamModeDefined())
						drv.m->setParamMode(AbsSyn::modeType::BOX);
					
					if (drv.m->getIterations() < 0)
					{
						yy::parser::error(yystack_[0].location, "Iteration number must be defined");
						YYERROR;
					}
				}
#line 997 "parser/src/parser.cpp"
    break;

  case 3:
#line 135 "parser/parser.yy"
                {
			if (drv.m->getVarMode() == AbsSyn::modeType::BOX)
						drv.m->defaultDirections();
					
			if (drv.m->getParamMode() == AbsSyn::modeType::BOX)
				drv.m->defaultParamDirections();
			
			if (drv.m->getVarMode() != AbsSyn::modeType::POLY)
							drv.m->defaultTemplate();
			
			if (!drv.m->isTransModeDefined())
				drv.m->setTransMode(AbsSyn::transType::AFO);
			
			if(!drv.m->isAlphaDefined())
				drv.m->setAlpha(0.5);
		}
#line 1018 "parser/src/parser.cpp"
    break;

  case 4:
#line 151 "parser/parser.yy"
                      { yy::parser::error(yystack_[0].location, "Empty file"); YYERROR; }
#line 1024 "parser/src/parser.cpp"
    break;

  case 5:
#line 153 "parser/parser.yy"
                         {}
#line 1030 "parser/src/parser.cpp"
    break;

  case 6:
#line 154 "parser/parser.yy"
                                                                    {}
#line 1036 "parser/src/parser.cpp"
    break;

  case 7:
#line 157 "parser/parser.yy"
                                                {
							if (drv.m->isProblemDefined())
							{
								yy::parser::error(yystack_[0].location, "Problem has already been defined");
								YYERROR;
							}
							drv.m->setProblem(yystack_[1].value.as < AbsSyn::problemType > ());
						}
#line 1049 "parser/src/parser.cpp"
    break;

  case 8:
#line 166 "parser/parser.yy"
                                                {
							if (drv.m->isVarModeDefined())
							{
								yy::parser::error(yystack_[0].location, "Variable modality has already been defined");
								YYERROR;
							}
							drv.m->setVarMode(yystack_[1].value.as < AbsSyn::modeType > ());
						}
#line 1062 "parser/src/parser.cpp"
    break;

  case 9:
#line 175 "parser/parser.yy"
                                                {
							if (drv.m->isParamModeDefined())
							{
								yy::parser::error(yystack_[0].location, "Parameter modality has already been defined");
								YYERROR;
							}
							drv.m->setParamMode(yystack_[1].value.as < AbsSyn::modeType > ());
						}
#line 1075 "parser/src/parser.cpp"
    break;

  case 10:
#line 184 "parser/parser.yy"
                                                {
							if (drv.m->getIterations() >= 0)
							{
								yy::parser::error(yystack_[0].location, "Iteration number already defined");
								YYERROR;
							}
							
							drv.m->setIterations(yystack_[1].value.as < int > ());
						}
#line 1089 "parser/src/parser.cpp"
    break;

  case 11:
#line 194 "parser/parser.yy"
                         {}
#line 1095 "parser/src/parser.cpp"
    break;

  case 12:
#line 195 "parser/parser.yy"
                                                                    {}
#line 1101 "parser/src/parser.cpp"
    break;

  case 13:
#line 198 "parser/parser.yy"
                                                {
							if (drv.m->getVarMode() != AbsSyn::modeType::BOX)
							{
								yy::parser::error(yylhs.location, "Cannot define variable bounds if modality is not 'boxes'");
								YYERROR;
							}
							
							for (unsigned i = 0; i < yystack_[3].value.as < std::vector<std::string> > ().size(); i++)
							{
								if (drv.m->isSymbolDefined(yystack_[3].value.as < std::vector<std::string> > ()[i]))
								{
									yy::parser::error(yystack_[3].location, "Symbol '" + yystack_[3].value.as < std::vector<std::string> > ()[i] + "' already defined");
									YYERROR;
								}
								drv.m->addVariable(new AbsSyn::Variable(yystack_[3].value.as < std::vector<std::string> > ()[i]));
								drv.m->addBounds(yystack_[1].value.as < std::pair<double, double> > ().first, yystack_[1].value.as < std::pair<double, double> > ().second);
							}
						}
#line 1124 "parser/src/parser.cpp"
    break;

  case 14:
#line 217 "parser/parser.yy"
                                                {
							if (drv.m->getVarMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(yylhs.location, "If variable modality is 'boxes', bounds must be provided");
								YYERROR;
							}
							
							for (unsigned i = 0; i < yystack_[1].value.as < std::vector<std::string> > ().size(); i++)
							{
								if (drv.m->isSymbolDefined(yystack_[1].value.as < std::vector<std::string> > ()[i]))
								{
									yy::parser::error(yystack_[1].location, "Symbol '" + yystack_[1].value.as < std::vector<std::string> > ()[i] + "' already defined");
									YYERROR;
								}
								drv.m->addVariable(new AbsSyn::Variable(yystack_[1].value.as < std::vector<std::string> > ()[i]));
							}
						}
#line 1146 "parser/src/parser.cpp"
    break;

  case 15:
#line 235 "parser/parser.yy"
                                                {
							if (drv.m->getParamMode() != AbsSyn::modeType::BOX)
							{
								yy::parser::error(yylhs.location, "Cannot define parameter bounds if modality is not 'boxes'");
								YYERROR;
							}
							
							for (unsigned i = 0; i < yystack_[3].value.as < std::vector<std::string> > ().size(); i++)
							{
								if (drv.m->isSymbolDefined(yystack_[3].value.as < std::vector<std::string> > ()[i]))
								{
									yy::parser::error(yystack_[3].location, "Symbol '" + yystack_[3].value.as < std::vector<std::string> > ()[i] + "' already defined");
									YYERROR;
								}
								drv.m->addParameter(new AbsSyn::Parameter(yystack_[3].value.as < std::vector<std::string> > ()[i]));
								drv.m->addParamBounds(yystack_[1].value.as < std::pair<double, double> > ().first, yystack_[1].value.as < std::pair<double, double> > ().second);
							}
						}
#line 1169 "parser/src/parser.cpp"
    break;

  case 16:
#line 254 "parser/parser.yy"
                                                {
							for (unsigned i = 0; i < yystack_[1].value.as < std::vector<std::string> > ().size(); i++)
							{
								if (drv.m->isSymbolDefined(yystack_[1].value.as < std::vector<std::string> > ()[i]))
								{
									yy::parser::error(yystack_[1].location, "Symbol '" + yystack_[1].value.as < std::vector<std::string> > ()[i] + "' already defined");
									YYERROR;
								}
								drv.m->addParameter(new AbsSyn::Parameter(yystack_[1].value.as < std::vector<std::string> > ()[i]));
							}
						}
#line 1185 "parser/src/parser.cpp"
    break;

  case 17:
#line 266 "parser/parser.yy"
                                                {
							if (!yystack_[1].value.as < AbsSyn::Expr * > ()->isNumeric())
							{
								yy::parser::error(yystack_[2].location, "Expression defining constant must be numeric");
								YYERROR;
							}
							if (drv.m->isSymbolDefined(yystack_[3].value.as < std::string > ()))
							{
								yy::parser::error(yystack_[3].location, "Symbol '" + yystack_[3].value.as < std::string > () + "' already defined");
								YYERROR;
							}
							drv.m->addConstant(new AbsSyn::Constant(yystack_[3].value.as < std::string > (), yystack_[1].value.as < AbsSyn::Expr * > ()->evaluate()));
						}
#line 1203 "parser/src/parser.cpp"
    break;

  case 18:
#line 280 "parser/parser.yy"
                                                {
							if (!drv.m->isSymbolDefined(yystack_[4].value.as < std::string > ()))
							{
								yy::parser::error(yystack_[4].location, "Unknown symbol name '" + yystack_[4].value.as < std::string > () + "'");
								YYERROR;
							}
							
							if (!drv.m->isVarDefined(yystack_[4].value.as < std::string > ()))
							{
								yy::parser::error(yystack_[4].location, "'" + yystack_[4].value.as < std::string > () + "' is not a variable");
								YYERROR;
							}
							
							AbsSyn::Variable *v = drv.m->getVar(yystack_[4].value.as < std::string > ());
							if (v->isDynamicDefined())
							{
								yy::parser::error(yylhs.location, "Redefinition of dynamic for variable '" + yystack_[4].value.as < std::string > () + "'");
								YYERROR;
							}
							v->setDynamic(yystack_[1].value.as < AbsSyn::Expr * > ());
						}
#line 1229 "parser/src/parser.cpp"
    break;

  case 19:
#line 302 "parser/parser.yy"
                                                {
							if (drv.m->getProblem() != AbsSyn::problemType::SYNTH)
							{
								yy::parser::error(yylhs.location, "Specification not required when problem is 'reachability'");
								YYERROR;
							}
							
							if (!yystack_[1].value.as < AbsSyn::Formula * > ()->simplify())
							{
								yy::parser::error(yystack_[1].location, "Negations of UNTIL are not allowed");
								YYERROR;
							}
							
							drv.m->addSpec(yystack_[1].value.as < AbsSyn::Formula * > ());
						}
#line 1249 "parser/src/parser.cpp"
    break;

  case 20:
#line 318 "parser/parser.yy"
                         {}
#line 1255 "parser/src/parser.cpp"
    break;

  case 21:
#line 319 "parser/parser.yy"
                                                                                {}
#line 1261 "parser/src/parser.cpp"
    break;

  case 22:
#line 322 "parser/parser.yy"
                                                {
							if (drv.m->getVarMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(yylhs.location, "Cannot define directions if variable modality is 'boxes'");
								YYERROR;
							}
							
							if (drv.m->templateRows() != 0)
							{
								yy::parser::error(yylhs.location, "Cannot define directions after defining template matrix");
								YYERROR;
							}
						}
#line 1279 "parser/src/parser.cpp"
    break;

  case 23:
#line 336 "parser/parser.yy"
                                                {
							if (drv.m->getVarMode() != AbsSyn::modeType::POLY)
							{
								yy::parser::error(yylhs.location, "Template matrix can be provided only if variable modality is 'polytopes'");
								YYERROR;
							}
							
							if (drv.m->directionsNum() == 0)
							{
								yy::parser::error(yylhs.location, "Template matrix cannot be provided before directions");
								YYERROR;
							}
						}
#line 1297 "parser/src/parser.cpp"
    break;

  case 24:
#line 350 "parser/parser.yy"
                                                {
							if (drv.m->getParamMode() == AbsSyn::modeType::BOX)
							{
								yy::parser::error(yylhs.location, "Cannot define parameter directions if parameter modality is 'boxes'");
								YYERROR;
							}
						}
#line 1309 "parser/src/parser.cpp"
    break;

  case 25:
#line 359 "parser/parser.yy"
                                                {
							if (drv.m->getVarNum() != yystack_[4].value.as < std::vector<double> > ().size())
							{
								yy::parser::error(yystack_[4].location, "A direction must have as many components as the number of variables");
								YYERROR;
							}
							
							drv.m->addDirection(yystack_[4].value.as < std::vector<double> > (), yystack_[1].value.as < std::pair<double, double> > ().first, yystack_[1].value.as < std::pair<double, double> > ().second);
						}
#line 1323 "parser/src/parser.cpp"
    break;

  case 26:
#line 370 "parser/parser.yy"
                                                {
							if (yystack_[0].value.as < std::vector<std::vector<int>> > ()[0].size() != drv.m->getVarNum())
							{
								yy::parser::error(yystack_[0].location, "template matrix must have as many columns as the number of variables");
								YYERROR;
							}
							
							if (drv.m->templateRows() != 0)
							{
								yy::parser::error(yylhs.location, "Redefinition of template matrix");
								YYERROR;
							}
							
							for (unsigned i = 0; i < yystack_[0].value.as < std::vector<std::vector<int>> > ().size(); i++)
							{
								for (unsigned j = 0; j < yystack_[0].value.as < std::vector<std::vector<int>> > ()[i].size(); j++)
								{
									if (yystack_[0].value.as < std::vector<std::vector<int>> > ()[i][j] < 0 || yystack_[0].value.as < std::vector<std::vector<int>> > ()[i][j] >= drv.m->directionsNum())
									{
										yy::parser::error(yystack_[0].location, "Unknown direction");
										YYERROR;
									}
								}
							}
							
							drv.m->setTemplate(yystack_[0].value.as < std::vector<std::vector<int>> > ());
						}
#line 1355 "parser/src/parser.cpp"
    break;

  case 27:
#line 399 "parser/parser.yy"
                                                {
							if (drv.m->getParamNum() != yystack_[4].value.as < std::vector<double> > ().size())
							{
								yy::parser::error(yystack_[4].location, "A paramter direction must have as many components as the number of parameters");
								YYERROR;
							}
							
							drv.m->addParamDirection(yystack_[4].value.as < std::vector<double> > (), yystack_[1].value.as < std::pair<double, double> > ().first, yystack_[1].value.as < std::pair<double, double> > ().second);
						}
#line 1369 "parser/src/parser.cpp"
    break;

  case 28:
#line 409 "parser/parser.yy"
                        { yylhs.value.as < AbsSyn::problemType > () = AbsSyn::problemType::REACH; }
#line 1375 "parser/src/parser.cpp"
    break;

  case 29:
#line 410 "parser/parser.yy"
                                                        { yylhs.value.as < AbsSyn::problemType > () = AbsSyn::problemType::SYNTH; }
#line 1381 "parser/src/parser.cpp"
    break;

  case 30:
#line 412 "parser/parser.yy"
                      { yylhs.value.as < AbsSyn::modeType > () = AbsSyn::modeType::BOX; }
#line 1387 "parser/src/parser.cpp"
    break;

  case 31:
#line 413 "parser/parser.yy"
                                                { yylhs.value.as < AbsSyn::modeType > () = AbsSyn::modeType::PARAL; }
#line 1393 "parser/src/parser.cpp"
    break;

  case 32:
#line 414 "parser/parser.yy"
                                               { yylhs.value.as < AbsSyn::modeType > () = AbsSyn::modeType::POLY; }
#line 1399 "parser/src/parser.cpp"
    break;

  case 33:
#line 416 "parser/parser.yy"
                        { yylhs.value.as < std::vector<std::string> > () = vector<string>{yystack_[0].value.as < std::string > ()}; }
#line 1405 "parser/src/parser.cpp"
    break;

  case 34:
#line 417 "parser/parser.yy"
                                                              { yystack_[2].value.as < std::vector<std::string> > ().push_back(yystack_[0].value.as < std::string > ()); yylhs.value.as < std::vector<std::string> > () = yystack_[2].value.as < std::vector<std::string> > (); }
#line 1411 "parser/src/parser.cpp"
    break;

  case 35:
#line 420 "parser/parser.yy"
                                                                {
									if (!yystack_[3].value.as < AbsSyn::Expr * > ()->isNumeric())
									{
										yy::parser::error(yystack_[3].location, "Intervals require numeric expressions");
										YYERROR;
									}
									
									if (!yystack_[1].value.as < AbsSyn::Expr * > ()->isNumeric())
									{
										yy::parser::error(yystack_[1].location, "Intervals require numeric expressions");
										YYERROR;
									}
									
									int x1 = (int) yystack_[3].value.as < AbsSyn::Expr * > ()->evaluate();
									int x2 = (int) yystack_[1].value.as < AbsSyn::Expr * > ()->evaluate();
								
									if (x2 < x1)
									{
										yy::parser::error(yylhs.location, "Right endpoint must be greater than  or equal to the left one");
										YYERROR;
									}
									
									yylhs.value.as < std::pair<int, int> > () = pair<int, int>(x1, x2);
								}
#line 1440 "parser/src/parser.cpp"
    break;

  case 36:
#line 446 "parser/parser.yy"
                                                                {
									if (!yystack_[3].value.as < AbsSyn::Expr * > ()->isNumeric())
									{
										yy::parser::error(yystack_[3].location, "Intervals require numeric expressions");
										YYERROR;
									}
									
									if (!yystack_[1].value.as < AbsSyn::Expr * > ()->isNumeric())
									{
										yy::parser::error(yystack_[1].location, "Intervals require numeric expressions");
										YYERROR;
									}
									
									double x1 = yystack_[3].value.as < AbsSyn::Expr * > ()->evaluate();
									double x2 = yystack_[1].value.as < AbsSyn::Expr * > ()->evaluate();
									
									if (x2 < x1)
									{
										yy::parser::error(yylhs.location, "Right endpoint must be greater than or equal to the left one");
										YYERROR;
									}
									
									yylhs.value.as < std::pair<double, double> > () = pair<double, double>(x1, x2);
								}
#line 1469 "parser/src/parser.cpp"
    break;

  case 37:
#line 471 "parser/parser.yy"
                         { yylhs.value.as < double > () = yystack_[0].value.as < double > (); }
#line 1475 "parser/src/parser.cpp"
    break;

  case 38:
#line 472 "parser/parser.yy"
                                                  { yylhs.value.as < double > () = (double) yystack_[0].value.as < int > (); }
#line 1481 "parser/src/parser.cpp"
    break;

  case 39:
#line 474 "parser/parser.yy"
                                { yylhs.value.as < AbsSyn::Expr * > () = new AbsSyn::Expr(yystack_[0].value.as < double > ()); }
#line 1487 "parser/src/parser.cpp"
    break;

  case 40:
#line 476 "parser/parser.yy"
                                {
					if (!drv.m->isSymbolDefined(yystack_[0].value.as < std::string > ()))
					{
						yy::parser::error(yystack_[0].location, "unknown symbol name");
						YYERROR;
					}
				
					if (drv.m->isConstDefined(yystack_[0].value.as < std::string > ()))
						yylhs.value.as < AbsSyn::Expr * > () = new AbsSyn::Expr(drv.m->getConst(yystack_[0].value.as < std::string > ())->getValue());
					else
						yylhs.value.as < AbsSyn::Expr * > () = new AbsSyn::Expr(yystack_[0].value.as < std::string > ());
				}
#line 1504 "parser/src/parser.cpp"
    break;

  case 41:
#line 488 "parser/parser.yy"
                                                { yylhs.value.as < AbsSyn::Expr * > () = yystack_[2].value.as < AbsSyn::Expr * > ()->mul(yystack_[0].value.as < AbsSyn::Expr * > ()); }
#line 1510 "parser/src/parser.cpp"
    break;

  case 42:
#line 489 "parser/parser.yy"
                                                { yylhs.value.as < AbsSyn::Expr * > () = yystack_[2].value.as < AbsSyn::Expr * > ()->sum(yystack_[0].value.as < AbsSyn::Expr * > ()); }
#line 1516 "parser/src/parser.cpp"
    break;

  case 43:
#line 490 "parser/parser.yy"
                                                { yylhs.value.as < AbsSyn::Expr * > () = yystack_[2].value.as < AbsSyn::Expr * > ()->sub(yystack_[0].value.as < AbsSyn::Expr * > ()); }
#line 1522 "parser/src/parser.cpp"
    break;

  case 44:
#line 491 "parser/parser.yy"
                                                        { yylhs.value.as < AbsSyn::Expr * > () = yystack_[0].value.as < AbsSyn::Expr * > ()->neg(); }
#line 1528 "parser/src/parser.cpp"
    break;

  case 45:
#line 492 "parser/parser.yy"
                                               { yylhs.value.as < AbsSyn::Expr * > () = yystack_[1].value.as < AbsSyn::Expr * > ();}
#line 1534 "parser/src/parser.cpp"
    break;

  case 46:
#line 494 "parser/parser.yy"
                                  { yylhs.value.as < std::vector<int> > () = yystack_[1].value.as < std::vector<int> > (); }
#line 1540 "parser/src/parser.cpp"
    break;

  case 47:
#line 495 "parser/parser.yy"
                                                  { yy::parser::error(yylhs.location, "Matrix row cannot be empty"); YYERROR; }
#line 1546 "parser/src/parser.cpp"
    break;

  case 48:
#line 498 "parser/parser.yy"
                                        {
						if (!yystack_[0].value.as < AbsSyn::Expr * > ()->isNumeric())
						{
							yy::parser::error(yystack_[0].location, "Expression must be numeric only");
							YYERROR;
						}
						yylhs.value.as < std::vector<double> > () = std::vector<double>{yystack_[0].value.as < AbsSyn::Expr * > ()->evaluate()};
					}
#line 1559 "parser/src/parser.cpp"
    break;

  case 49:
#line 507 "parser/parser.yy"
                                        {
						if (!yystack_[0].value.as < AbsSyn::Expr * > ()->isNumeric())
						{
							yy::parser::error(yystack_[0].location, "Expression must be numeric only");
							YYERROR;
						}
						yystack_[2].value.as < std::vector<double> > ().push_back(yystack_[0].value.as < AbsSyn::Expr * > ()->evaluate());
						yylhs.value.as < std::vector<double> > () = yystack_[2].value.as < std::vector<double> > ();
					}
#line 1573 "parser/src/parser.cpp"
    break;

  case 50:
#line 518 "parser/parser.yy"
                                        {
						if (yystack_[0].value.as < int > () < 0 || yystack_[0].value.as < int > () >= drv.m->directionsNum())
						{
							yy::parser::error(yystack_[0].location, "Unknown direction");
							YYERROR;
						}
						yylhs.value.as < std::vector<int> > () = std::vector<int>{yystack_[0].value.as < int > ()};
					}
#line 1586 "parser/src/parser.cpp"
    break;

  case 51:
#line 527 "parser/parser.yy"
                                        {
						if (yystack_[0].value.as < int > () < 0 || yystack_[0].value.as < int > () >= drv.m->directionsNum())
						{
							yy::parser::error(yystack_[0].location, "Unknown direction");
							YYERROR;
						}
						yystack_[2].value.as < std::vector<int> > ().push_back(yystack_[0].value.as < int > ());
						yylhs.value.as < std::vector<int> > () = yystack_[2].value.as < std::vector<int> > ();
					}
#line 1600 "parser/src/parser.cpp"
    break;

  case 52:
#line 537 "parser/parser.yy"
                                          { yylhs.value.as < std::vector<std::vector<int>> > () = yystack_[1].value.as < std::vector<std::vector<int>> > (); }
#line 1606 "parser/src/parser.cpp"
    break;

  case 53:
#line 538 "parser/parser.yy"
                                                          { yy::parser::error(yylhs.location, "Matrix cannot be empty"); YYERROR; }
#line 1612 "parser/src/parser.cpp"
    break;

  case 54:
#line 540 "parser/parser.yy"
                                    { yylhs.value.as < std::vector<std::vector<int>> > () = std::vector<std::vector<int>>{yystack_[0].value.as < std::vector<int> > ()}; }
#line 1618 "parser/src/parser.cpp"
    break;

  case 55:
#line 541 "parser/parser.yy"
                                                                        { yystack_[2].value.as < std::vector<std::vector<int>> > ().push_back(yystack_[0].value.as < std::vector<int> > ()); yylhs.value.as < std::vector<std::vector<int>> > () = yystack_[2].value.as < std::vector<std::vector<int>> > (); }
#line 1624 "parser/src/parser.cpp"
    break;

  case 56:
#line 543 "parser/parser.yy"
                                { yylhs.value.as < AbsSyn::Formula * > () = new AbsSyn::Formula(yystack_[0].value.as < AbsSyn::Expr * > ()->sub(yystack_[2].value.as < AbsSyn::Expr * > ())); }
#line 1630 "parser/src/parser.cpp"
    break;

  case 57:
#line 544 "parser/parser.yy"
                                                                         { yylhs.value.as < AbsSyn::Formula * > () = new AbsSyn::Formula(yystack_[0].value.as < AbsSyn::Expr * > ()->sub(yystack_[2].value.as < AbsSyn::Expr * > ())); }
#line 1636 "parser/src/parser.cpp"
    break;

  case 58:
#line 545 "parser/parser.yy"
                                                                        { yylhs.value.as < AbsSyn::Formula * > () = new AbsSyn::Formula(yystack_[2].value.as < AbsSyn::Expr * > ()->sub(yystack_[0].value.as < AbsSyn::Expr * > ())); }
#line 1642 "parser/src/parser.cpp"
    break;

  case 59:
#line 546 "parser/parser.yy"
                                                                         { yylhs.value.as < AbsSyn::Formula * > () = new AbsSyn::Formula(yystack_[2].value.as < AbsSyn::Expr * > ()->sub(yystack_[0].value.as < AbsSyn::Expr * > ())); }
#line 1648 "parser/src/parser.cpp"
    break;

  case 60:
#line 548 "parser/parser.yy"
                                                        {
								AbsSyn::Formula *f1 = new AbsSyn::Formula(yystack_[2].value.as < AbsSyn::Expr * > ()->sub(yystack_[0].value.as < AbsSyn::Expr * > ()));
								AbsSyn::Formula *f2 = new AbsSyn::Formula(yystack_[0].value.as < AbsSyn::Expr * > ()->sub(yystack_[2].value.as < AbsSyn::Expr * > ()));
								yylhs.value.as < AbsSyn::Formula * > () = f1->conj(f2);
							}
#line 1658 "parser/src/parser.cpp"
    break;

  case 61:
#line 553 "parser/parser.yy"
                                                                                                        { yylhs.value.as < AbsSyn::Formula * > () = yystack_[2].value.as < AbsSyn::Formula * > ()->conj(yystack_[0].value.as < AbsSyn::Formula * > ()); }
#line 1664 "parser/src/parser.cpp"
    break;

  case 62:
#line 554 "parser/parser.yy"
                                                                                                        { yylhs.value.as < AbsSyn::Formula * > () = yystack_[2].value.as < AbsSyn::Formula * > ()->disj(yystack_[0].value.as < AbsSyn::Formula * > ()); }
#line 1670 "parser/src/parser.cpp"
    break;

  case 63:
#line 555 "parser/parser.yy"
                                                                                                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[0].value.as < AbsSyn::Formula * > ()->neg(); }
#line 1676 "parser/src/parser.cpp"
    break;

  case 64:
#line 556 "parser/parser.yy"
                                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[1].value.as < AbsSyn::Formula * > (); }
#line 1682 "parser/src/parser.cpp"
    break;

  case 65:
#line 558 "parser/parser.yy"
                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[0].value.as < AbsSyn::Formula * > ()->always(yystack_[1].value.as < std::pair<int, int> > ()); }
#line 1688 "parser/src/parser.cpp"
    break;

  case 66:
#line 559 "parser/parser.yy"
                                                                                                        { yylhs.value.as < AbsSyn::Formula * > () = yystack_[0].value.as < AbsSyn::Formula * > ()->eventually(yystack_[1].value.as < std::pair<int, int> > ()); }
#line 1694 "parser/src/parser.cpp"
    break;

  case 67:
#line 560 "parser/parser.yy"
                                                                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[3].value.as < AbsSyn::Formula * > ()->until(yystack_[1].value.as < std::pair<int, int> > (), yystack_[0].value.as < AbsSyn::Formula * > ()); }
#line 1700 "parser/src/parser.cpp"
    break;

  case 68:
#line 561 "parser/parser.yy"
                                                                                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[2].value.as < AbsSyn::Formula * > ()->conj(yystack_[0].value.as < AbsSyn::Formula * > ()); }
#line 1706 "parser/src/parser.cpp"
    break;

  case 69:
#line 562 "parser/parser.yy"
                                                                                                                                { yylhs.value.as < AbsSyn::Formula * > () = yystack_[2].value.as < AbsSyn::Formula * > ()->disj(yystack_[0].value.as < AbsSyn::Formula * > ()); }
#line 1712 "parser/src/parser.cpp"
    break;

  case 70:
#line 563 "parser/parser.yy"
                                                                                                                                                                        { yylhs.value.as < AbsSyn::Formula * > () = yystack_[0].value.as < AbsSyn::Formula * > ()->neg(); }
#line 1718 "parser/src/parser.cpp"
    break;

  case 71:
#line 564 "parser/parser.yy"
                                                                               { yylhs.value.as < AbsSyn::Formula * > () = yystack_[1].value.as < AbsSyn::Formula * > (); }
#line 1724 "parser/src/parser.cpp"
    break;

  case 72:
#line 566 "parser/parser.yy"
                         {}
#line 1730 "parser/src/parser.cpp"
    break;

  case 73:
#line 567 "parser/parser.yy"
                                                                    {}
#line 1736 "parser/src/parser.cpp"
    break;

  case 74:
#line 570 "parser/parser.yy"
                                {
					if (drv.m->isTransModeDefined())
					{
						yy::parser::error(yylhs.location, "Transformation type already defined");
						YYERROR;
					}
					drv.m->setTransMode(yystack_[1].value.as < AbsSyn::transType > ());
				}
#line 1749 "parser/src/parser.cpp"
    break;

  case 75:
#line 579 "parser/parser.yy"
                                {
					if (drv.m->isDecompositionDefined())
					{
						yy::parser::error(yylhs.location, "Decomposition option already defined");
						YYERROR;
					}
					
					drv.m->setDecomposition();
				}
#line 1763 "parser/src/parser.cpp"
    break;

  case 76:
#line 589 "parser/parser.yy"
                                {
					if (drv.m->isAlphaDefined())
					{
						yy::parser::error(yylhs.location, "Alpha already defined");
						YYERROR;
					}
					
					if (yystack_[1].value.as < double > () > 1)
					{
						yy::parser::error(yystack_[1].location, "Alpha must be between 0 and 1");
						YYERROR;
					}
					
					drv.m->setAlpha(yystack_[1].value.as < double > ());
				}
#line 1783 "parser/src/parser.cpp"
    break;

  case 77:
#line 605 "parser/parser.yy"
                { yylhs.value.as < AbsSyn::transType > () = AbsSyn::transType::AFO; }
#line 1789 "parser/src/parser.cpp"
    break;

  case 78:
#line 606 "parser/parser.yy"
                                              { yylhs.value.as < AbsSyn::transType > () = AbsSyn::transType::OFO; }
#line 1795 "parser/src/parser.cpp"
    break;


#line 1799 "parser/src/parser.cpp"

            default:
              break;
            }
        }
#if YY_EXCEPTIONS
      catch (const syntax_error& yyexc)
        {
          YYCDEBUG << "Caught exception: " << yyexc.what() << '\n';
          error (yyexc);
          YYERROR;
        }
#endif // YY_EXCEPTIONS
      YY_SYMBOL_PRINT ("-> $$ =", yylhs);
      yypop_ (yylen);
      yylen = 0;
      YY_STACK_PRINT ();

      // Shift the result of the reduction.
      yypush_ (YY_NULLPTR, YY_MOVE (yylhs));
    }
    goto yynewstate;


  /*--------------------------------------.
  | yyerrlab -- here on detecting error.  |
  `--------------------------------------*/
  yyerrlab:
    // If not already recovering from an error, report this error.
    if (!yyerrstatus_)
      {
        ++yynerrs_;
        error (yyla.location, yysyntax_error_ (yystack_[0].state, yyla));
      }


    yyerror_range[1].location = yyla.location;
    if (yyerrstatus_ == 3)
      {
        /* If just tried and failed to reuse lookahead token after an
           error, discard it.  */

        // Return failure if at end of input.
        if (yyla.type_get () == yyeof_)
          YYABORT;
        else if (!yyla.empty ())
          {
            yy_destroy_ ("Error: discarding", yyla);
            yyla.clear ();
          }
      }

    // Else will try to reuse lookahead token after shifting the error token.
    goto yyerrlab1;


  /*---------------------------------------------------.
  | yyerrorlab -- error raised explicitly by YYERROR.  |
  `---------------------------------------------------*/
  yyerrorlab:
    /* Pacify compilers when the user code never invokes YYERROR and
       the label yyerrorlab therefore never appears in user code.  */
    if (false)
      YYERROR;

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYERROR.  */
    yypop_ (yylen);
    yylen = 0;
    goto yyerrlab1;


  /*-------------------------------------------------------------.
  | yyerrlab1 -- common code for both syntax error and YYERROR.  |
  `-------------------------------------------------------------*/
  yyerrlab1:
    yyerrstatus_ = 3;   // Each real token shifted decrements this.
    {
      stack_symbol_type error_token;
      for (;;)
        {
          yyn = yypact_[+yystack_[0].state];
          if (!yy_pact_value_is_default_ (yyn))
            {
              yyn += yy_error_token_;
              if (0 <= yyn && yyn <= yylast_ && yycheck_[yyn] == yy_error_token_)
                {
                  yyn = yytable_[yyn];
                  if (0 < yyn)
                    break;
                }
            }

          // Pop the current state because it cannot handle the error token.
          if (yystack_.size () == 1)
            YYABORT;

          yyerror_range[1].location = yystack_[0].location;
          yy_destroy_ ("Error: popping", yystack_[0]);
          yypop_ ();
          YY_STACK_PRINT ();
        }

      yyerror_range[2].location = yyla.location;
      YYLLOC_DEFAULT (error_token.location, yyerror_range, 2);

      // Shift the error token.
      error_token.state = state_type (yyn);
      yypush_ ("Shifting", YY_MOVE (error_token));
    }
    goto yynewstate;


  /*-------------------------------------.
  | yyacceptlab -- YYACCEPT comes here.  |
  `-------------------------------------*/
  yyacceptlab:
    yyresult = 0;
    goto yyreturn;


  /*-----------------------------------.
  | yyabortlab -- YYABORT comes here.  |
  `-----------------------------------*/
  yyabortlab:
    yyresult = 1;
    goto yyreturn;


  /*-----------------------------------------------------.
  | yyreturn -- parsing is finished, return the result.  |
  `-----------------------------------------------------*/
  yyreturn:
    if (!yyla.empty ())
      yy_destroy_ ("Cleanup: discarding lookahead", yyla);

    /* Do not reclaim the symbols of the rule whose action triggered
       this YYABORT or YYACCEPT.  */
    yypop_ (yylen);
    while (1 < yystack_.size ())
      {
        yy_destroy_ ("Cleanup: popping", yystack_[0]);
        yypop_ ();
      }

    return yyresult;
  }
#if YY_EXCEPTIONS
    catch (...)
      {
        YYCDEBUG << "Exception caught: cleaning lookahead and stack\n";
        // Do not try to display the values of the reclaimed symbols,
        // as their printers might throw an exception.
        if (!yyla.empty ())
          yy_destroy_ (YY_NULLPTR, yyla);

        while (1 < yystack_.size ())
          {
            yy_destroy_ (YY_NULLPTR, yystack_[0]);
            yypop_ ();
          }
        throw;
      }
#endif // YY_EXCEPTIONS
  }

  void
  parser::error (const syntax_error& yyexc)
  {
    error (yyexc.location, yyexc.what ());
  }

  // Generate an error message.
  std::string
  parser::yysyntax_error_ (state_type yystate, const symbol_type& yyla) const
  {
    // Number of reported tokens (one for the "unexpected", one per
    // "expected").
    std::ptrdiff_t yycount = 0;
    // Its maximum.
    enum { YYERROR_VERBOSE_ARGS_MAXIMUM = 5 };
    // Arguments of yyformat.
    char const *yyarg[YYERROR_VERBOSE_ARGS_MAXIMUM];

    /* There are many possibilities here to consider:
       - If this state is a consistent state with a default action, then
         the only way this function was invoked is if the default action
         is an error action.  In that case, don't check for expected
         tokens because there are none.
       - The only way there can be no lookahead present (in yyla) is
         if this state is a consistent state with a default action.
         Thus, detecting the absence of a lookahead is sufficient to
         determine that there is no unexpected or expected token to
         report.  In that case, just report a simple "syntax error".
       - Don't assume there isn't a lookahead just because this state is
         a consistent state with a default action.  There might have
         been a previous inconsistent state, consistent state with a
         non-default action, or user semantic action that manipulated
         yyla.  (However, yyla is currently not documented for users.)
       - Of course, the expected token list depends on states to have
         correct lookahead information, and it depends on the parser not
         to perform extra reductions after fetching a lookahead from the
         scanner and before detecting a syntax error.  Thus, state merging
         (from LALR or IELR) and default reductions corrupt the expected
         token list.  However, the list is correct for canonical LR with
         one exception: it will still contain any token that will not be
         accepted due to an error action in a later state.
    */
    if (!yyla.empty ())
      {
        symbol_number_type yytoken = yyla.type_get ();
        yyarg[yycount++] = yytname_[yytoken];

        int yyn = yypact_[+yystate];
        if (!yy_pact_value_is_default_ (yyn))
          {
            /* Start YYX at -YYN if negative to avoid negative indexes in
               YYCHECK.  In other words, skip the first -YYN actions for
               this state because they are default actions.  */
            int yyxbegin = yyn < 0 ? -yyn : 0;
            // Stay within bounds of both yycheck and yytname.
            int yychecklim = yylast_ - yyn + 1;
            int yyxend = yychecklim < yyntokens_ ? yychecklim : yyntokens_;
            for (int yyx = yyxbegin; yyx < yyxend; ++yyx)
              if (yycheck_[yyx + yyn] == yyx && yyx != yy_error_token_
                  && !yy_table_value_is_error_ (yytable_[yyx + yyn]))
                {
                  if (yycount == YYERROR_VERBOSE_ARGS_MAXIMUM)
                    {
                      yycount = 1;
                      break;
                    }
                  else
                    yyarg[yycount++] = yytname_[yyx];
                }
          }
      }

    char const* yyformat = YY_NULLPTR;
    switch (yycount)
      {
#define YYCASE_(N, S)                         \
        case N:                               \
          yyformat = S;                       \
        break
      default: // Avoid compiler warnings.
        YYCASE_ (0, YY_("syntax error"));
        YYCASE_ (1, YY_("syntax error, unexpected %s"));
        YYCASE_ (2, YY_("syntax error, unexpected %s, expecting %s"));
        YYCASE_ (3, YY_("syntax error, unexpected %s, expecting %s or %s"));
        YYCASE_ (4, YY_("syntax error, unexpected %s, expecting %s or %s or %s"));
        YYCASE_ (5, YY_("syntax error, unexpected %s, expecting %s or %s or %s or %s"));
#undef YYCASE_
      }

    std::string yyres;
    // Argument number.
    std::ptrdiff_t yyi = 0;
    for (char const* yyp = yyformat; *yyp; ++yyp)
      if (yyp[0] == '%' && yyp[1] == 's' && yyi < yycount)
        {
          yyres += yytnamerr_ (yyarg[yyi++]);
          ++yyp;
        }
      else
        yyres += *yyp;
    return yyres;
  }


  const signed char parser::yypact_ninf_ = -51;

  const signed char parser::yytable_ninf_ = -1;

  const short
  parser::yypact_[] =
  {
     102,   -51,   -38,   -36,   -32,    -9,    54,     7,   -51,   119,
     171,   171,    -4,   -51,   157,   -51,   -51,   -51,    37,   -51,
     -51,   -51,    56,    73,    82,    45,    45,    90,   110,    94,
     157,   -51,   -51,   -51,   -51,   -51,   -51,   -12,     3,   116,
     106,    39,   -51,   164,   130,   -51,   127,   130,   -51,   100,
     126,    39,   147,   147,   100,    39,   -51,   -51,   -51,   -51,
      97,   -23,    42,   144,   158,   146,   -51,   -51,   -51,   -51,
       8,   100,   145,   -51,   148,   100,   109,   160,   -51,   159,
     100,    47,    47,   -51,    80,    52,    -7,   100,   100,   100,
     100,   100,   100,   100,   100,    39,    39,   -51,    47,    47,
     147,   100,   152,   100,   -51,    81,   -51,    22,   -51,   -51,
     138,   -51,   100,    25,    47,    47,   -51,   -51,   -51,   -51,
     -51,   149,   149,   149,   161,   161,   -51,   149,   149,   169,
     -51,   170,   -51,    47,   149,   -15,    70,   -51,    38,   114,
     -51,   150,   100,   117,   100,   -51,   104,   -51,   186,   100,
     -25,   -51,   -51,   -14,   187,   -51,   -51,   155,   156,   122,
     -51,   125,   130,   149,   -51,   -51,    44,   -51,   162,   130,
     -51,   -51,   -51,   -51,   163,   -51,   165,   -51,   166,   -51,
     -51,   -51
  };

  const signed char
  parser::yydefact_[] =
  {
       0,     4,     0,     0,     0,     0,     0,     2,     5,     0,
       0,     0,     0,     1,     0,     6,    28,    29,     0,    30,
      31,    32,     0,     0,     0,     0,     0,     0,     0,     0,
      20,    11,     7,     8,     9,    10,    33,     0,     0,     0,
       0,     0,    12,    72,     0,    14,     0,     0,    16,     0,
       0,     0,     0,     0,     0,     0,    40,    38,    37,    39,
       0,     0,     0,     0,     0,     0,    21,    22,    23,    24,
       0,     0,     0,    34,     0,     0,     0,     0,    70,    63,
       0,     0,     0,    44,     0,     0,     0,     0,     0,     0,
       0,     0,     0,     0,     0,     0,     0,    19,     0,     0,
       0,     0,     0,     0,     3,     0,    73,     0,    13,    15,
       0,    17,     0,     0,     0,     0,    65,    66,    45,    71,
      64,    59,    57,    60,    43,    42,    41,    58,    56,    68,
      69,    61,    62,     0,    48,     0,     0,    26,     0,     0,
      75,     0,     0,     0,     0,    63,     0,    67,     0,     0,
       0,    53,    54,     0,     0,    77,    78,     0,     0,     0,
      18,     0,     0,    49,    47,    50,     0,    52,     0,     0,
      74,    76,    36,    35,     0,    46,     0,    55,     0,    25,
      51,    27
  };

  const short
  parser::yypgoto_[] =
  {
     -51,   -51,   195,   -47,   -50,   -51,   181,   -48,   105,    41,
     -51,   -51,   -51,   -39,   -33,   -51,   -51,   -51,   -51,   203,
     -51,   182,   -51,   -51,   -51,   -51,   -51,   -51,   -51
  };

  const short
  parser::yydefgoto_[] =
  {
      -1,    18,    22,    72,    81,    59,    37,    60,   135,   152,
     166,   137,   153,    61,    62,   157,     6,    14,     7,     8,
      30,    31,    43,    66,    67,    68,    69,    70,   106
  };

  const unsigned char
  parser::yytable_[] =
  {
      74,    76,    44,    82,    95,    96,    83,    84,   104,     9,
       2,    10,    78,     3,     4,    11,    85,    47,    79,   164,
      98,    99,    86,   107,     5,    97,   165,   110,   100,   105,
     167,   148,   113,   120,   149,   168,    45,    46,    12,   121,
     122,   123,   124,   125,   126,   127,   128,    24,   116,   117,
     133,    48,    46,   134,    13,   134,   129,   130,    90,    91,
      92,    90,    91,    92,   143,   131,   132,    84,    51,    98,
      99,   142,    52,    53,   144,    54,   114,   100,    55,    95,
      96,   145,   146,    54,   154,    32,   115,   149,   175,    56,
      57,    58,   119,   176,   159,    36,   161,    56,    57,    58,
     147,   163,     1,   139,    33,     2,   140,   141,     3,     4,
      87,    88,    89,   150,   151,   174,    90,    91,    92,     5,
     118,    34,   178,    16,    17,    93,    94,    87,    88,    89,
      35,    98,    99,    90,    91,    92,    54,   155,   156,    75,
      39,    41,    93,    94,   120,    90,    91,    92,    49,    40,
      56,    57,    58,    90,    91,    92,    50,   111,    90,    91,
      92,    90,    91,    92,   172,   160,    77,   173,    25,    26,
      27,    71,    28,    29,    90,    91,    92,    73,   118,    19,
      20,    21,    63,    64,    65,    90,    91,    92,    80,   101,
     102,   103,   112,   108,   100,   136,   109,    96,    99,    92,
     162,   169,   158,   170,   171,   150,    23,    38,   138,   177,
      15,   179,    42,     0,   181,     0,   180
  };

  const short
  parser::yycheck_[] =
  {
      47,    49,    14,    53,    27,    28,    54,    55,     0,    47,
       3,    47,    51,     6,     7,    47,    55,    14,    51,    44,
      27,    28,    55,    71,    17,    48,    51,    75,    35,    21,
      44,    46,    80,    40,    49,    49,    48,    49,    47,    87,
      88,    89,    90,    91,    92,    93,    94,    51,    81,    82,
     100,    48,    49,   101,     0,   103,    95,    96,    36,    37,
      38,    36,    37,    38,   112,    98,    99,   115,    29,    27,
      28,    49,    33,    34,    49,    36,    29,    35,    39,    27,
      28,   114,   115,    36,    46,    48,    39,    49,    44,    50,
      51,    52,    40,    49,   142,    50,   144,    50,    51,    52,
     133,   149,     0,    22,    48,     3,    25,    26,     6,     7,
      30,    31,    32,    43,    44,   162,    36,    37,    38,    17,
      40,    48,   169,     4,     5,    45,    46,    30,    31,    32,
      48,    27,    28,    36,    37,    38,    36,    23,    24,    39,
      50,    47,    45,    46,    40,    36,    37,    38,    32,    39,
      50,    51,    52,    36,    37,    38,    50,    48,    36,    37,
      38,    36,    37,    38,    42,    48,    40,    42,    11,    12,
      13,    41,    15,    16,    36,    37,    38,    50,    40,     8,
       9,    10,    18,    19,    20,    36,    37,    38,    41,    45,
      32,    45,    32,    48,    35,    43,    48,    28,    28,    38,
      14,    14,    52,    48,    48,    43,    11,    26,   103,   168,
       7,    48,    30,    -1,    48,    -1,    51
  };

  const signed char
  parser::yystos_[] =
  {
       0,     0,     3,     6,     7,    17,    70,    72,    73,    47,
      47,    47,    47,     0,    71,    73,     4,     5,    55,     8,
       9,    10,    56,    56,    51,    11,    12,    13,    15,    16,
      74,    75,    48,    48,    48,    48,    50,    60,    60,    50,
      39,    47,    75,    76,    14,    48,    49,    14,    48,    32,
      50,    29,    33,    34,    36,    39,    50,    51,    52,    59,
      61,    67,    68,    18,    19,    20,    77,    78,    79,    80,
      81,    41,    57,    50,    57,    39,    61,    40,    67,    68,
      41,    58,    58,    61,    61,    67,    68,    30,    31,    32,
      36,    37,    38,    45,    46,    27,    28,    48,    27,    28,
      35,    45,    32,    45,     0,    21,    82,    61,    48,    48,
      61,    48,    32,    61,    29,    39,    68,    68,    40,    40,
      40,    61,    61,    61,    61,    61,    61,    61,    61,    67,
      67,    68,    68,    58,    61,    62,    43,    65,    62,    22,
      25,    26,    49,    61,    49,    68,    68,    68,    46,    49,
      43,    44,    63,    66,    46,    23,    24,    69,    52,    61,
      48,    61,    14,    61,    44,    51,    64,    44,    49,    14,
      48,    48,    42,    42,    57,    44,    49,    63,    57,    48,
      51,    48
  };

  const signed char
  parser::yyr1_[] =
  {
       0,    54,    71,    70,    70,    72,    72,    73,    73,    73,
      73,    74,    74,    75,    75,    75,    75,    75,    75,    75,
      76,    76,    77,    77,    77,    78,    79,    80,    55,    55,
      56,    56,    56,    60,    60,    58,    57,    59,    59,    61,
      61,    61,    61,    61,    61,    61,    63,    63,    62,    62,
      64,    64,    65,    65,    66,    66,    68,    68,    68,    68,
      68,    68,    68,    68,    68,    67,    67,    67,    67,    67,
      67,    67,    81,    81,    82,    82,    82,    69,    69
  };

  const signed char
  parser::yyr2_[] =
  {
       0,     2,     0,     6,     1,     1,     2,     4,     4,     4,
       4,     1,     2,     5,     3,     5,     3,     5,     7,     4,
       0,     2,     1,     1,     1,     7,     3,     7,     1,     1,
       1,     1,     1,     1,     3,     5,     5,     1,     1,     1,
       1,     3,     3,     3,     2,     3,     3,     2,     1,     3,
       1,     3,     3,     2,     1,     3,     3,     3,     3,     3,
       3,     3,     3,     2,     3,     3,     3,     4,     3,     3,
       2,     3,     0,     2,     4,     2,     4,     1,     1
  };



  // YYTNAME[SYMBOL-NUM] -- String name of the symbol SYMBOL-NUM.
  // First, the terminals, then, starting at \a yyntokens_, nonterminals.
  const char*
  const parser::yytname_[] =
  {
  "END", "error", "$undefined", "PROB", "REACH", "SYNTH", "VARMODE",
  "PARAMMODE", "BOX", "PARAL", "POLY", "VAR", "PARAM", "CONST", "IN",
  "DYN", "SPEC", "ITER", "DIR", "TEMPL", "PDIR", "OPT", "TRANS", "AFO",
  "OFO", "DECOMP", "ALPHA", "\"&&\"", "\"||\"", "\"!\"", "\"<=\"",
  "\">=\"", "\"=\"", "\"G\"", "\"F\"", "\"U\"", "\"-\"", "\"+\"", "\"*\"",
  "\"(\"", "\")\"", "\"[\"", "\"]\"", "\"{\"", "\"}\"", "\"<\"", "\">\"",
  "\":\"", "\";\"", "\",\"", "IDENT", "INTEGER", "DOUBLE", "UMINUS",
  "$accept", "problemType", "modeType", "doubleInterval", "intInterval",
  "number", "identList", "expr", "numList", "matrixRow", "intList",
  "_matrix", "rowList", "path_formula", "state_formula", "transType", "s",
  "$@1", "headerList", "header", "symbolList", "symbol", "matricesList",
  "matrices", "direction", "template", "paramDir", "footerList", "footer", YY_NULLPTR
  };

#if YYDEBUG
  const short
  parser::yyrline_[] =
  {
       0,   115,   115,   114,   151,   153,   154,   156,   165,   174,
     183,   194,   195,   197,   216,   234,   253,   265,   279,   301,
     318,   319,   321,   335,   349,   358,   369,   398,   409,   410,
     412,   413,   414,   416,   417,   419,   445,   471,   472,   474,
     475,   488,   489,   490,   491,   492,   494,   495,   497,   506,
     517,   526,   537,   538,   540,   541,   543,   544,   545,   546,
     547,   553,   554,   555,   556,   558,   559,   560,   561,   562,
     563,   564,   566,   567,   569,   578,   588,   605,   606
  };

  // Print the state stack on the debug stream.
  void
  parser::yystack_print_ ()
  {
    *yycdebug_ << "Stack now";
    for (stack_type::const_iterator
           i = yystack_.begin (),
           i_end = yystack_.end ();
         i != i_end; ++i)
      *yycdebug_ << ' ' << int (i->state);
    *yycdebug_ << '\n';
  }

  // Report on the debug stream that the rule \a yyrule is going to be reduced.
  void
  parser::yy_reduce_print_ (int yyrule)
  {
    int yylno = yyrline_[yyrule];
    int yynrhs = yyr2_[yyrule];
    // Print the symbols being reduced, and their result.
    *yycdebug_ << "Reducing stack by rule " << yyrule - 1
               << " (line " << yylno << "):\n";
    // The symbols being reduced.
    for (int yyi = 0; yyi < yynrhs; yyi++)
      YY_SYMBOL_PRINT ("   $" << yyi + 1 << " =",
                       yystack_[(yynrhs) - (yyi + 1)]);
  }
#endif // YYDEBUG


} // yy
#line 2308 "parser/src/parser.cpp"

#line 608 "parser/parser.yy"


void
yy::parser::error (const location_type& l, const std::string& m)
{
	std::cerr << "Error at line " << l.end.line << ", column " << l.end.column << ": " << m << '\n';
}
