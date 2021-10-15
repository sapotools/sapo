#include "driver.h"
#include "parser.h"

driver::driver()
	: trace_parsing (false), trace_scanning (false)
{
	data = new AbsSyn::InputData();
}

int driver::parse (const std::string &f)
{
	file = f;
	location.initialize (&file);
	scan_begin();
	yy::parser parse(*this);
	parse.set_debug_level(trace_parsing);
	int res = parse();
	scan_end();
	if (res == 0 && data->check())
		return 0;
	else
		return 1;
}

driver::~driver()
{
	delete data;
}