#include "driver.h"
#include "parser.h"

driver::driver():
    data(AbsSyn::InputData()), ctx(AbsSyn::Context()),
    errors(false), trace_parsing(false), trace_scanning(false)
{
}

int driver::parse(const std::string &f)
{
  file = f;
  location.initialize(&file);
  scan_begin();
  yy::parser parse(*this);
//  parse.set_debug_level(trace_parsing);
  int res = parse();
  scan_end();
  if (res == 0 && !errors && data.check(ctx))
    return 0;
  else
    return 1;
}

driver::~driver() {}


void driver::warning(const yy::location &l, const std::string &m)
{
	std::cerr << "\033[1;95mWarning\033[0m at line " << l.end.line << ", column " << l.end.column << ": " << m << '\n';
}

void driver::error(const yy::location &l, const std::string &m)
{
	std::cerr << "\033[1;31mError\033[0m at line " << l.end.line << ", column " << l.end.column << ": " << m << '\n';
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
	for (int i = l.begin.line+1; i <= l.end.line + 1; i++) {
		getline(file, line);
		int current_digits = floor(log10(i)) + 1;
		
		cerr << "  ";
		for (int j = 0; j < digits - current_digits; j++) {
			cerr << " ";
		}
		cerr << i << " | " << line << endl;
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
