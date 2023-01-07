#include "driver.h"
#include "parser.h"

driver::driver():
    data(AbsSyn::InputData()), errors(false), trace_parsing(false),
    trace_scanning(false)
{
}

int driver::parse(const std::string &f)
{
  file = f;
  location.initialize(&file);
  scan_begin();
  yy::parser parse(*this);
  parse.set_debug_level(trace_parsing);
  int res = parse();
  scan_end();

  try {
    data.optimize_boundaries();
  } catch (std::domain_error &e) {
    // the parameter set is empty

    std::cerr << "The parameter set cannot be empty." << std::endl;

    return 1;
  }
  if (res == 0 && !errors && data.check())
    return 0;
  else
    return 1;
}

driver::~driver() {}

void driver::warning(const yy::location &l, const std::string &m,
                     const std::string filename)
{
  std::cerr << "\033[1;95mWarning\033[0m at line " << l.end.line << ", column "
            << l.end.column << ": " << m << '\n';
  printError(l, filename);
}

void driver::error(const yy::location &l, const std::string &m,
                   const std::string filename)
{
  std::cerr << "\033[1;31mError\033[0m at line " << l.end.line << ", column "
            << l.end.column << ": " << m << '\n';
  printError(l, filename);
  errors = true;
}

void driver::missingSemicolon(const yy::location &l,
                              const std::string filename)
{
  yy::location loc(l);
  loc.begin.line = l.end.line;
  loc.begin.column = l.end.column;
  loc.end.column = l.end.column + 1;
  error(loc, "Missing \";\"", filename);
}

void driver::printError(const yy::location &l, const std::string filename)
{
  using namespace std;

  // open file
  ifstream file;
  if (filename.empty() || filename == "-") {
    // TODO: check if we can read from stdin
    return;
  } else {
    file.open(filename);
  }

  // skip to correct line
  string line;
  for (int i = 0; i < l.begin.line - 1; i++) {
    getline(file, line);
  }

  // print useful lines
  int digits = floor(log10(l.begin.line)) + 2;
  for (int i = l.begin.line + 1; i <= l.end.line + 1; i++) {
    getline(file, line);
    std::string space_line = "";
    for (unsigned i = 0; i < line.size(); i++) {
      if (line[i] == '\t') {
        space_line += "    ";
      } else {
        space_line += line[i];
      }
    }
    int current_digits = floor(log10(i)) + 1;

    cerr << "  ";
    for (int j = 0; j < digits - current_digits; j++) {
      cerr << " ";
    }
    cerr << (i - 1) << " | " << space_line << endl;
  }

  if (l.begin.line == l.end.line) {
    cerr << "  ";
    for (int i = 0; i < digits; i++) {
      cerr << " ";
    }
    cerr << " |";
    for (int i = 0; i < l.begin.column; i++) {
      cerr << " ";
    }
    for (int i = l.begin.column; i < l.end.column; i++) {
      cerr << "^";
    }
    cerr << endl;
  }

  // close file
  file.close();

  return;
}
